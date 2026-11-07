use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};
use crate::db::index_item::{index_item_size, split_record};
use crate::db::with_path;
use crate::io::smiles::reader::molecule::parse_molecule;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};
use std::ops::RangeInclusive;
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct SearchResult {
    pub line: String,
    pub similarity: f32,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SearchMode {
    /// The `limit` most similar records, best first.
    Best,
    /// The first `limit` records in file order; stops reading early.
    First,
}

/// Returns records whose similarity to `query` is in `min_similarity..=max_similarity`.
pub fn similarity_search(
    smiles_file: &str,
    query: &str,
    min_similarity: f32,
    max_similarity: f32,
    limit: usize,
    mode: SearchMode,
) -> io::Result<Vec<SearchResult>> {
    let query = parse_molecule(query)
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("could not parse query {query}"),
            )
        })?
        .1;

    // smiles file
    let fis = File::open(smiles_file).map_err(with_path(smiles_file))?;
    let mut reader = BufReader::new(fis);

    // fingerprints file
    let fp_path = smiles_file.to_owned() + ".fp";
    let mut fif = File::open(&fp_path).map_err(with_path(&fp_path))?;

    // the query must use the index's fingerprint parameters
    let mut header = [0u8; INDEX_HEADER_SIZE];
    fif.read_exact(&mut header).map_err(with_path(&fp_path))?;
    let header = IndexHeader::from_bytes(&header)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{fp_path}: {e}")))?;
    // in the stored byte order, so records are scored without decoding them
    let query_words: Vec<u64> = query
        .ecfp(header.radius, header.bits)
        .to_bytes()
        .as_chunks::<8>()
        .0
        .iter()
        .map(|chunk| u64::from_le_bytes(*chunk))
        .collect();
    let item_size = index_item_size(header.bits);

    // stored offsets are only valid for the exact file that was indexed
    let stale = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg + ": re-index");
    // shortcut: compares lengths only, so an edit that keeps the size goes unnoticed; store a hash if that matters
    let smiles_len = reader
        .get_ref()
        .metadata()
        .map_err(with_path(smiles_file))?
        .len();
    if smiles_len != header.smiles_len {
        return Err(stale(format!("{smiles_file} changed since it was indexed")));
    }
    // shortcut: a program truncating the index during a search crashes it with SIGBUS; read
    // blocks into buffers instead (slower on Linux, faster on Windows) if that ever happens
    #[expect(
        unsafe_code,
        reason = "ringo-index never writes an index in place: it renames a finished file over \
                  it, so the mapped bytes stay fixed"
    )]
    // SAFETY: as above, nothing in ringo changes the mapped file while it is mapped
    let map = unsafe { Mmap::map(&fif) }.map_err(with_path(&fp_path))?;
    let expected_len = header
        .count
        .checked_mul(item_size as u64)
        .and_then(|records| records.checked_add(INDEX_HEADER_SIZE as u64));
    if expected_len != Some(map.len() as u64) {
        return Err(stale(format!(
            "{fp_path} does not hold the {} records its header lists",
            header.count
        )));
    }

    let mut hits = scan(
        &map[INDEX_HEADER_SIZE..],
        item_size,
        &query_words,
        min_similarity..=max_similarity,
        limit,
        mode,
    )?;
    if mode == SearchMode::Best {
        // best first; the sort is stable, so ties keep file order
        // shortcut: keeps every hit in memory, use a bounded BinaryHeap if low thresholds on huge files run out of memory
        hits.sort_by(|a, b| b.0.total_cmp(&a.0));
    }
    hits.truncate(limit);

    hits.into_iter()
        .map(|(similarity, position)| {
            reader.seek(io::SeekFrom::Start(position as u64))?;

            let mut line = String::new();
            reader
                .read_line(&mut line)
                .map_err(with_path(smiles_file))?;
            Ok(SearchResult { line, similarity })
        })
        .collect()
}

/// `(similarity, line offset)` pairs
type Hits = Vec<(f32, usize)>;

/// `(similarity, line offset)` of the `item_size`-byte `records` whose similarity is in
/// `range`, in file order. In `first` mode only the first `limit` hits are complete.
///
/// Threads take 1 MiB blocks of whole records in file order and score them in place, so the
/// finished blocks are always a prefix of the file: once they hold `limit` hits, the first
/// `limit` hits are among them.
fn scan(
    records: &[u8],
    item_size: usize,
    query_words: &[u64],
    range: RangeInclusive<f32>,
    limit: usize,
    mode: SearchMode,
) -> io::Result<Hits> {
    let query_ones = query_words.iter().map(|w| w.count_ones() as usize).sum();
    let block_len = ((1 << 20) / item_size).max(1) * item_size;
    let next_block = AtomicUsize::new(0);
    let found = AtomicUsize::new(0);
    let scan_blocks = || -> io::Result<Vec<(usize, Hits)>> {
        let mut scanned = Vec::new();
        while mode == SearchMode::Best || found.load(Ordering::Relaxed) < limit {
            let block = next_block.fetch_add(1, Ordering::Relaxed);
            let start = block * block_len;
            if start >= records.len() {
                break;
            }
            let block_bytes = &records[start..records.len().min(start + block_len)];

            let mut hits = Vec::new();
            for record in block_bytes.chunks_exact(item_size) {
                let (position, fingerprint) = split_record(record)?;
                let Some(similarity) = tanimoto(fingerprint, query_words, query_ones) else {
                    continue;
                };
                if range.contains(&similarity) {
                    hits.push((similarity, position));
                    if mode == SearchMode::First && hits.len() >= limit {
                        break;
                    }
                }
            }
            found.fetch_add(hits.len(), Ordering::Relaxed);
            scanned.push((block, hits));
        }
        Ok(scanned)
    };
    let scanned = (0..rayon::current_num_threads())
        .into_par_iter()
        .map(|_| scan_blocks())
        .collect::<io::Result<Vec<_>>>()?;
    let mut scanned: Vec<_> = scanned.into_iter().flatten().collect();
    scanned.sort_unstable_by_key(|(block, _)| *block);
    Ok(scanned.into_iter().flat_map(|(_, hits)| hits).collect())
}

/// Tanimoto of stored fingerprint bytes against the query, in one pass and without
/// allocating, the hot loop of a search. `None` for an empty record: a line the indexer
/// could not parse, whose similarity 0 would pass a range from 0.
#[expect(
    clippy::cast_precision_loss,
    reason = "counts above 2^24 bits only blur digits far below what a score shows"
)]
fn tanimoto(record: &[u8], query_words: &[u64], query_ones: usize) -> Option<f32> {
    let (mut ones, mut common) = (0, 0);
    for (chunk, query) in record.as_chunks::<8>().0.iter().zip(query_words) {
        let word = u64::from_le_bytes(*chunk);
        ones += word.count_ones() as usize;
        common += (word & query).count_ones() as usize;
    }
    (ones > 0).then(|| common as f32 / (ones + query_ones - common) as f32)
}

#[cfg(test)]
mod test {
    use crate::db::index::test::{copy_sample_molecules, TestDir};
    use crate::db::index::{index_file, DEFAULT_BATCH_LINES};
    use crate::db::index_item::index_item_size;
    use crate::db::search::SearchMode::{Best, First};
    use crate::db::search::{scan, similarity_search, tanimoto};
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;

    #[test]
    fn test_tanimoto_matches_bitset() {
        let fp = |s| parse_molecule(s).unwrap().1.ecfp(2, 1024);
        let smiles = ["C", "CCO", "c1ccccc1O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"];
        for a in smiles {
            for b in smiles {
                let (a, b) = (fp(a), fp(b));
                let words: Vec<u64> = b
                    .to_bytes()
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|chunk| u64::from_le_bytes(*chunk))
                    .collect();
                assert_eq!(
                    tanimoto(&a.to_bytes(), &words, b.0.count_ones(..)),
                    Some(tanimoto_bitset(&a.0, &b.0))
                );
            }
        }
        assert_eq!(tanimoto(&[0; 128], &[u64::MAX; 16], 1024), None);
    }

    #[test]
    fn test_similarity_search() {
        let dir = TestDir::new("test_similarity_search");
        let smiles_file = copy_sample_molecules(dir.path());
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(
            smiles_file,
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            0.7,
            1.0,
            100,
            Best,
        )
        .unwrap();
        assert_eq!(results.len(), 1);
        assert!(results[0].line.starts_with("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"));
        assert_eq!(results[0].similarity, 1.0);
        let results = similarity_search(
            smiles_file,
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            0.3,
            1.0,
            100,
            Best,
        )
        .unwrap();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_search_returns_best_matches_first() {
        let dir = TestDir::new("test_search_best_matches");
        let smiles_file = dir.path().join("m.smi");
        // the exact match is last, so a search that stops early never sees it
        std::fs::write(&smiles_file, "CCO\nCCN\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 1, Best).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].line.trim_end(), "CCC");
        assert_eq!(results[0].similarity, 1.0);

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).unwrap();
        assert_eq!(results.len(), 3);
        assert!(results
            .windows(2)
            .all(|w| w[0].similarity >= w[1].similarity));

        // max_similarity drops the exact match
        let results = similarity_search(smiles_file, "CCC", 0.0, 0.99, 10, Best).unwrap();
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.line.trim_end() != "CCC"));
    }

    #[test]
    fn test_search_first_mode_keeps_file_order() {
        let dir = TestDir::new("test_search_first_mode");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCN\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 2, First).unwrap();
        let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
        assert_eq!(lines, ["CCO", "CCN"]);
    }

    #[test]
    fn test_search_skips_unparsed_records() {
        let dir = TestDir::new("test_search_skips_unparsed");
        let smiles_file = dir.path().join("m.smi");
        // the unsupported allene line gets an empty fingerprint; it used to fill the one result
        std::fs::write(&smiles_file, "NC(Br)=[C@AL1]=C(O)C\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        for mode in [First, Best] {
            let results = similarity_search(smiles_file, "CCO", 0.0, 1.0, 10, mode).unwrap();
            let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
            assert_eq!(lines, ["CCC"], "{mode:?}");
        }
    }

    #[test]
    fn test_search_reads_across_blocks() {
        let dir = TestDir::new("test_search_across_blocks");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCN\nCCC\nCCCC\nCCCCO\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        // 400 KB records, so a 1 MiB block holds 2 and the last block holds 1
        index_file(smiles_file, 2, 3_200_000, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).unwrap();
        assert_eq!(results.len(), 5);
        assert_eq!(results[0].line.trim_end(), "CCC");
        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 3, First).unwrap();
        let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
        assert_eq!(lines, ["CCO", "CCN", "CCC"]);
    }

    #[test]
    fn test_scan_keeps_file_order_across_threads() {
        // 512-bit records, 14,563 to a 1 MiB block, so 7 blocks; every 1000th record matches
        // and the others are empty, which never match
        let count = 100_000;
        let mut records = Vec::new();
        for i in 0..count {
            records.extend((i as u64).to_le_bytes());
            records.extend([if i % 1000 == 0 { 0xff } else { 0 }; 64]);
        }
        let positions = |limit, mode| -> Vec<usize> {
            let hits = scan(
                &records,
                index_item_size(512),
                &[u64::MAX; 8],
                0.0..=1.0,
                limit,
                mode,
            );
            hits.unwrap()
                .into_iter()
                .map(|(_, position)| position)
                .collect()
        };

        let all: Vec<usize> = (0..count).step_by(1000).collect();
        assert_eq!(positions(10, Best), all);
        // limits within the first block and across several
        for limit in [1, 20, 50, 100] {
            assert_eq!(positions(limit, First)[..limit], all[..limit], "{limit}");
        }
        assert_eq!(positions(0, First), Vec::<usize>::new());
    }

    #[test]
    fn test_search_uses_index_fingerprint_params() {
        let dir = TestDir::new("test_search_index_params");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 3, 1024, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).unwrap();
        let fp = |s| parse_molecule(s).unwrap().1.ecfp(3, 1024);
        assert_eq!(results[0].similarity, 1.0);
        assert_eq!(
            results[1].similarity,
            tanimoto_bitset(&fp("CCO").0, &fp("CCC").0)
        );
    }

    #[test]
    fn test_search_rejects_index_without_header() {
        let dir = TestDir::new("test_search_old_index");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCC\n").unwrap();
        // a pre-header index: one 72-byte record starting with offset 0
        std::fs::write(dir.path().join("m.smi.fp"), [0u8; 72]).unwrap();

        let e = similarity_search(smiles_file.to_str().unwrap(), "CCC", 0.0, 1.0, 10, Best)
            .err()
            .unwrap();
        assert!(e.to_string().contains("re-index"), "{e}");
    }

    #[test]
    fn test_search_rejects_stale_index() {
        let dir = TestDir::new("test_search_stale_index");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        let fp_file = smiles_file.to_owned() + ".fp";
        let search = || similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).err();

        // SMILES file edited after indexing: stored offsets would point at the wrong lines
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();
        std::fs::write(smiles_file, "CCCO\nCCC\n").unwrap();
        let e = search().unwrap();
        assert!(
            e.to_string().contains("changed since it was indexed"),
            "{e}"
        );
        assert!(e.to_string().contains("re-index"), "{e}");

        // index missing its last byte, or with bytes appended
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();
        let fp = std::fs::read(&fp_file).unwrap();
        for bad in [&fp[..fp.len() - 1], &[fp.as_slice(), &[0]].concat()] {
            std::fs::write(&fp_file, bad).unwrap();
            let e = search().unwrap();
            assert!(e.to_string().contains("records its header lists"), "{e}");
        }

        // the untouched index still works
        std::fs::write(&fp_file, &fp).unwrap();
        assert!(search().is_none());
    }

    #[test]
    fn test_search_errors_instead_of_panicking() {
        let dir = TestDir::new("test_search_errors");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        let search = |query| similarity_search(smiles_file, query, 0.0, 1.0, 10, Best).err();

        // no index yet
        let e = search("CCC").unwrap();
        assert_eq!(e.kind(), std::io::ErrorKind::NotFound);
        assert!(e.to_string().contains("m.smi.fp"), "{e}");

        // index cut off inside its header
        std::fs::write(smiles_file.to_owned() + ".fp", [0u8; 4]).unwrap();
        assert_eq!(
            search("CCC").unwrap().kind(),
            std::io::ErrorKind::UnexpectedEof
        );

        // a query the reader cannot parse
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();
        assert_eq!(
            search("C1CC").unwrap().kind(),
            std::io::ErrorKind::InvalidInput
        );
    }

    #[test]
    fn test_search_returns_correct_lines_for_lf_and_crlf() {
        for (name, ending) in [("lf", "\n"), ("crlf", "\r\n")] {
            let dir = TestDir::new(&format!("test_search_line_endings_{name}"));
            let smiles_file = dir.path().join("m.smi");
            let content = ["CCO ethanol", "CCN ethylamine", "CCC propane"].join(ending) + ending;
            std::fs::write(&smiles_file, content).unwrap();
            let smiles_file = smiles_file.to_str().unwrap();
            index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

            let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).unwrap();
            let mut lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
            lines.sort_unstable();
            assert_eq!(
                lines,
                ["CCC propane", "CCN ethylamine", "CCO ethanol"],
                "{name}"
            );
        }
    }
}
