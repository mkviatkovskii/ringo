use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};
use crate::db::index_item::{index_item_size, split_record};
use crate::db::with_path;
use crate::io::smiles::reader::molecule::parse_molecule_with_validation;
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
    /// The first `limit` records in index order, by bit count and then by line; stops
    /// reading early.
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
    similarity_search_with_validation(
        smiles_file,
        query,
        min_similarity,
        max_similarity,
        limit,
        mode,
        false,
    )
}

/// As `similarity_search`, optionally validating the query's main-group valence.
/// Indexed records follow the policy chosen when building the index; search does not revalidate them.
pub fn similarity_search_with_validation(
    smiles_file: &str,
    query: &str,
    min_similarity: f32,
    max_similarity: f32,
    limit: usize,
    mode: SearchMode,
    validate: bool,
) -> io::Result<Vec<SearchResult>> {
    let query = parse_molecule_with_validation(query, validate)
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("could not parse query {query} (validation enabled: {validate})"),
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
    let records_start = header.records_start();
    let expected_len = header
        .count
        .checked_mul(item_size as u64)
        .and_then(|records| records.checked_add(records_start as u64));
    if expected_len != Some(map.len() as u64) {
        return Err(stale(format!(
            "{fp_path} does not hold the {} records its header lists",
            header.count
        )));
    }

    // the groups that can reach min_similarity: the bound rises up to the query's bit count
    // and falls after it, so they are one run of records
    let query_ones = query_words.iter().map(|w| w.count_ones() as usize).sum();
    let bad_table = || {
        stale(format!(
            "{fp_path}: group table does not add up to {} records",
            header.count
        ))
    };
    let (mut first, mut last, mut start) = (None, 0, 0usize);
    let table = &map[INDEX_HEADER_SIZE..records_start];
    for (ones, group_len) in (header.min_ones..).zip(table.as_chunks::<8>().0) {
        let end = usize::try_from(u64::from_le_bytes(*group_len))
            .ok()
            .and_then(|group_len| start.checked_add(group_len))
            .ok_or_else(bad_table)?;
        if can_reach(ones, query_ones, min_similarity) {
            first.get_or_insert(start);
            last = end;
        }
        start = end;
    }
    if start as u64 != header.count {
        return Err(bad_table());
    }
    let first = first.unwrap_or(last);

    let mut hits = scan(
        &map[records_start + first * item_size..records_start + last * item_size],
        item_size,
        &query_words,
        query_ones,
        min_similarity..=max_similarity,
        limit,
        mode,
    )?;
    if mode == SearchMode::Best {
        // best first, ties in line order
        // shortcut: keeps every hit in memory, use a bounded BinaryHeap if low thresholds on huge files run out of memory
        hits.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
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
/// `range`, in the order of `records`. In `first` mode only the first `limit` hits are complete.
///
/// Threads take 1 MiB blocks of whole records in order and score them in place, so the
/// finished blocks are always a prefix of `records`: once they hold `limit` hits, the first
/// `limit` hits are among them.
fn scan(
    records: &[u8],
    item_size: usize,
    query_words: &[u64],
    query_ones: usize,
    range: RangeInclusive<f32>,
    limit: usize,
    mode: SearchMode,
) -> io::Result<Hits> {
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

/// Whether a record with `ones` bits set can score `min_similarity` against a query with
/// `query_ones`: at best it shares every bit of the smaller one. Divides like `tanimoto`, and
/// rounding keeps the order of exact quotients, so it never rules out a record that matches.
#[expect(
    clippy::cast_precision_loss,
    reason = "counts above 2^24 bits round the same way in tanimoto"
)]
fn can_reach(ones: usize, query_ones: usize, min_similarity: f32) -> bool {
    ones > 0 && ones.min(query_ones) as f32 / ones.max(query_ones) as f32 >= min_similarity
}

#[cfg(test)]
mod test {
    use crate::db::index::test::{copy_sample_molecules, TestDir};
    use crate::db::index::{index_file, DEFAULT_BATCH_LINES};
    use crate::db::index_item::index_item_size;
    use crate::db::search::SearchMode::{Best, First};
    use crate::db::search::{can_reach, scan, similarity_search, tanimoto};
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
    fn test_can_reach_keeps_every_possible_score() {
        // every score a record can have, used as the threshold, which is where rounding
        // could drop a group
        for query_ones in 0..=100 {
            for ones in 1..=100 {
                for common in 0..=ones.min(query_ones) {
                    let mut record = vec![0u8; 128];
                    let mut query = [0u64; 16];
                    for bit in 0..ones {
                        record[bit / 8] |= 1 << (bit % 8);
                    }
                    // the query shares `common` bits with the record, the rest lie after them
                    for bit in (ones - common)..(ones - common + query_ones) {
                        query[bit / 64] |= 1 << (bit % 64);
                    }
                    let score = tanimoto(&record, &query, query_ones).unwrap();
                    assert!(
                        can_reach(ones, query_ones, score),
                        "{ones} bits, query {query_ones}, common {common}"
                    );
                }
            }
        }
        assert!(!can_reach(0, 10, 0.0));
        assert!(!can_reach(10, 30, 0.34));
        assert!(!can_reach(30, 10, 0.34));
    }

    #[test]
    fn test_search_matches_full_scan() {
        let dir = TestDir::new("test_search_full_scan");
        let smiles_file = dir.path().join("m.smi");
        let smiles = [
            "C",
            "CCO",
            "CC(=O)O",
            "c1ccccc1",
            "c1ccccc1O",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(C)NCC(O)COc1cccc2ccccc12",
            "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
            "C1CCCCC1",
            "C(",
            "OC(=O)CCCCCCCCCCCCCCCCC",
            "Nc1ccc(cc1)S(=O)(=O)N",
            "CC(=O)Nc1ccc(O)cc1",
        ];
        std::fs::write(&smiles_file, smiles.join("\n")).unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        let fp = |s| parse_molecule(s).ok().map(|(_, m)| m.ecfp(2, 512));
        for query in [
            "CCO",
            "c1ccccc1O",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(=O)Nc1ccc(O)cc1",
        ] {
            let query_fp = fp(query).unwrap();
            for min_similarity in [0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.95, 1.0] {
                // score every line, best first and ties in line order
                let mut expected: Vec<(f32, &str)> = smiles
                    .iter()
                    .filter_map(|line| {
                        let line_fp = fp(line).filter(|f| f.0.count_ones(..) > 0)?;
                        let similarity = tanimoto_bitset(&line_fp.0, &query_fp.0);
                        (similarity >= min_similarity).then_some((similarity, *line))
                    })
                    .collect();
                expected.sort_by(|a, b| b.0.total_cmp(&a.0));

                let results =
                    similarity_search(smiles_file, query, min_similarity, 1.0, 100, Best).unwrap();
                let results: Vec<(f32, &str)> = results
                    .iter()
                    .map(|r| (r.similarity, r.line.trim_end()))
                    .collect();
                assert_eq!(results, expected, "{query} at {min_similarity}");
            }
        }
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

    /// `lines` in index order: by bit count, then in line order
    fn index_order<'a>(lines: &[&'a str], bits: usize) -> Vec<&'a str> {
        let mut sorted = lines.to_vec();
        sorted.sort_by_key(|line| {
            parse_molecule(line)
                .unwrap()
                .1
                .ecfp(2, bits)
                .0
                .count_ones(..)
        });
        sorted
    }

    #[test]
    fn test_search_first_mode_keeps_index_order() {
        let dir = TestDir::new("test_search_first_mode");
        let smiles_file = dir.path().join("m.smi");
        let smiles = ["CCO", "CCN", "CCC", "c1ccccc1O", "CC"];
        std::fs::write(&smiles_file, smiles.join("\n")).unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 4, First).unwrap();
        let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
        assert_eq!(lines, index_order(&smiles, 512)[..4]);
    }

    #[test]
    fn test_search_skips_unparsed_records() {
        let dir = TestDir::new("test_search_skips_unparsed");
        let smiles_file = dir.path().join("m.smi");
        // the unsupported bond gets an empty fingerprint; it used to fill the one result
        std::fs::write(&smiles_file, "C~C\nCCC\n").unwrap();
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
        let smiles = ["CCO", "CCN", "CCC", "CCCC", "CCCCO"];
        std::fs::write(&smiles_file, smiles.join("\n")).unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        // 400 KB records, so a 1 MiB block holds 2 and the last block holds 1
        index_file(smiles_file, 2, 3_200_000, DEFAULT_BATCH_LINES).unwrap();

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best).unwrap();
        assert_eq!(results.len(), 5);
        assert_eq!(results[0].line.trim_end(), "CCC");
        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 3, First).unwrap();
        let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
        assert_eq!(lines, index_order(&smiles, 3_200_000)[..3]);
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
                512,
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
