use std::fs::File;
use std::io::{BufRead};
use crate::ringo::molecule::smiles::reader::molecule::parse_molecule;
use crate::ringo::ringo::index_item::IndexItem;

fn index(smiles_file: &str) -> Vec<IndexItem> {
    // open file for reading
    let fi = File::open(smiles_file).expect("Could not open file");

    let mut result = Vec::new();
    // open binary file for index
    let mut offset = 0;
    // let mut fo = File::create(smiles_file.to_owned() + ".fp");
    for line in std::io::BufReader::new(fi).lines() {
        let line = line.unwrap();
        let molecule = parse_molecule(&line).unwrap().1;
        // let ecfp = molecule.ecfp(2, 512);
        // write ecfp and offset to binary file
        result.push(IndexItem::new(offset, molecule.ecfp(2, 512)));

        offset += line.len() + 1;
    }
    return result;
}

#[test]
fn test_index() {
    let result = index("molecules.smi");
    assert_eq!(result.len(), 2);
}
