use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};
use crate::db::index_item::{index_item_size, IndexItem};
use crate::io::smiles::reader::molecule::parse_molecule;
use crate::math::similarity::tanimoto::tanimoto_bitset;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek};

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
) -> Vec<SearchResult> {
    let query = parse_molecule(query).unwrap().1;

    // smiles file
    let fis = File::open(smiles_file).expect("Could not open file");
    let mut reader = BufReader::new(fis);

    //fingerprings file
    let fif = File::open(smiles_file.to_owned() + ".fp").expect("Could not open file");
    let file_len = fif.metadata().unwrap().len();
    let mut buf_reader = BufReader::new(fif);

    // the query must use the index's fingerprint parameters
    let mut header = [0u8; INDEX_HEADER_SIZE];
    buf_reader
        .read_exact(&mut header)
        .expect("index file is too short");
    let header =
        IndexHeader::from_bytes(&header).unwrap_or_else(|e| panic!("{smiles_file}.fp: {e}"));
    let query_fp = query.ecfp(header.radius, header.bits);
    let item_size = index_item_size(header.bits);
    let index_count = (file_len - INDEX_HEADER_SIZE as u64) / item_size as u64;

    let mut hits = Vec::new();
    let mut buf = vec![0u8; item_size];

    for _ in 0..index_count {
        // read index item from file
        buf_reader.read_exact(&mut buf).unwrap();

        // decode index item
        let index_item = IndexItem::from_bytes(&buf);

        // calculate similarity
        let similarity = tanimoto_bitset(&index_item.fingerprint.0, &query_fp.0);
        // print similarity if it is greater than min_similarity
        if (min_similarity..=max_similarity).contains(&similarity) {
            hits.push((similarity, index_item.position));
            if mode == SearchMode::First && hits.len() >= limit {
                break;
            }
        }
    }

    if mode == SearchMode::Best {
        // best first; the sort is stable, so ties keep file order
        // shortcut: keeps every hit in memory, use a bounded BinaryHeap if low thresholds on huge files run out of memory
        hits.sort_by(|a, b| b.0.total_cmp(&a.0));
    }
    hits.truncate(limit);

    hits.into_iter()
        .map(|(similarity, position)| {
            reader
                .seek(std::io::SeekFrom::Start(position as u64))
                .unwrap();

            let mut line = String::new();
            reader.read_line(&mut line).unwrap();
            // println!("{i} {similarity} {position} {line}");
            SearchResult { line, similarity }
        })
        .collect()
}

#[cfg(test)]
mod test {
    use crate::db::index::index_file;
    use crate::db::index::test::{copy_sample_molecules, TestDir};
    use crate::db::search::similarity_search;
    use crate::db::search::SearchMode::{Best, First};
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;

    #[test]
    fn test_similarity_search() {
        let dir = TestDir::new("test_similarity_search");
        let smiles_file = copy_sample_molecules(dir.path());
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512);

        let results = similarity_search(
            smiles_file,
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            0.7,
            1.0,
            100,
            Best,
        );
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
        );
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_search_returns_best_matches_first() {
        let dir = TestDir::new("test_search_best_matches");
        let smiles_file = dir.path().join("m.smi");
        // the exact match is last, so a search that stops early never sees it
        std::fs::write(&smiles_file, "CCO\nCCN\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512);

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 1, Best);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].line.trim_end(), "CCC");
        assert_eq!(results[0].similarity, 1.0);

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best);
        assert_eq!(results.len(), 3);
        assert!(results
            .windows(2)
            .all(|w| w[0].similarity >= w[1].similarity));

        // max_similarity drops the exact match
        let results = similarity_search(smiles_file, "CCC", 0.0, 0.99, 10, Best);
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.line.trim_end() != "CCC"));
    }

    #[test]
    fn test_search_first_mode_keeps_file_order() {
        let dir = TestDir::new("test_search_first_mode");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCN\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512);

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 2, First);
        let lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
        assert_eq!(lines, ["CCO", "CCN"]);
    }

    #[test]
    fn test_search_uses_index_fingerprint_params() {
        let dir = TestDir::new("test_search_index_params");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 3, 1024);

        let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best);
        let fp = |s| parse_molecule(s).unwrap().1.ecfp(3, 1024);
        assert_eq!(results[0].similarity, 1.0);
        assert_eq!(
            results[1].similarity,
            tanimoto_bitset(&fp("CCO").0, &fp("CCC").0)
        );
    }

    #[test]
    #[should_panic(expected = "re-index")]
    fn test_search_rejects_index_without_header() {
        let dir = TestDir::new("test_search_old_index");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCC\n").unwrap();
        // a pre-header index: one 72-byte record starting with offset 0
        std::fs::write(dir.path().join("m.smi.fp"), [0u8; 72]).unwrap();

        similarity_search(smiles_file.to_str().unwrap(), "CCC", 0.0, 1.0, 10, Best);
    }

    #[test]
    fn test_search_returns_correct_lines_for_lf_and_crlf() {
        for (name, ending) in [("lf", "\n"), ("crlf", "\r\n")] {
            let dir = TestDir::new(&format!("test_search_line_endings_{name}"));
            let smiles_file = dir.path().join("m.smi");
            let content = ["CCO ethanol", "CCN ethylamine", "CCC propane"].join(ending) + ending;
            std::fs::write(&smiles_file, content).unwrap();
            let smiles_file = smiles_file.to_str().unwrap();
            index_file(smiles_file, 2, 512);

            let results = similarity_search(smiles_file, "CCC", 0.0, 1.0, 10, Best);
            let mut lines: Vec<&str> = results.iter().map(|r| r.line.trim_end()).collect();
            lines.sort();
            assert_eq!(
                lines,
                ["CCC propane", "CCN ethylamine", "CCO ethanol"],
                "{name}"
            );
        }
    }
}
