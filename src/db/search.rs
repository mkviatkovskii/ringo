use crate::math::similarity::tanimoto::tanimoto_bitset;
use crate::molecule::smiles::reader::molecule::parse_molecule;
use crate::db::index_item::IndexItem;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek};

pub struct SearchResult {
    pub line: String,
    pub similarity: f32,
}

pub fn similarity_search(
    smiles_file: &str,
    query: &str,
    min_similarity: f32,
    limit: usize,
) -> Vec<SearchResult> {
    let query = parse_molecule(query).unwrap().1;
    let query_fp = query.ecfp(2, 512);

    // smiles file
    let fis = File::open(&smiles_file).expect("Could not open file");
    let mut reader = BufReader::new(fis);

    //fingerprings file
    let fif = File::open(smiles_file.to_owned() + ".fp").expect("Could not open file");
    let file_len = fif.metadata().unwrap().len();
    let index_item_size = 72u8;
    let index_count = file_len / index_item_size as u64;
    let mut buf_reader = BufReader::new(fif);

    let mut results = Vec::new();

    for _ in 0..index_count {
        // read index item from file
        let mut buf = vec![0u8; index_item_size as usize];
        buf_reader.read_exact(&mut buf).unwrap();

        // decode index item
        let index_item: IndexItem = bincode::decode_from_slice(&buf, bincode::config::standard())
            .unwrap()
            .0;

        // calculate similarity
        let similarity = tanimoto_bitset(&index_item.fingerprint.0, &query_fp.0);
        // print similarity if it is greater than min_similarity
        if similarity >= min_similarity {
            let position = index_item.position;
            reader
                .seek(std::io::SeekFrom::Start(position as u64))
                .unwrap();

            let mut line = String::new();
            reader.read_line(&mut line).unwrap();
            // println!("{i} {similarity} {position} {line}");
            results.push(SearchResult {
                line: line,
                similarity: similarity,
            });

            if results.len() >= limit {
                break;
            }
        }
    }

    results
}

fn main() {
    println!("db-search");
}

#[cfg(test)]
mod test {
    use crate::db::index::index_file;
    use crate::db::search::similarity_search;

    #[test]
    fn test_similarity_search() {
        index_file("molecles.smi");
        let results = similarity_search("molecules.smi", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 0.7, 100);
        assert_eq!(results.len(), 1);
        assert!(results[0].line.starts_with("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"));
        assert_eq!(results[0].similarity, 1.0);
        let results = similarity_search("molecules.smi", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 0.5, 100);
        assert_eq!(results.len(), 2);
    }
}
