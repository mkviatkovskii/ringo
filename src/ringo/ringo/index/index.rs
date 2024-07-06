use crate::ringo::fingerprint::fingerprint::FINGERPRINT_SIZE;
use crate::ringo::molecule::smiles::reader::molecule::parse_molecule;
use crate::ringo::ringo::index::index_item::IndexItem;
use bincode::encode_into_slice;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

#[cfg(windows)]
const LINE_ENDING_LENGTH: usize = 2;
#[cfg(not(windows))]
const LINE_ENDING_LENGTH: usize = 1;

pub(crate) fn index(smiles_file: &str) {
    // open file for reading
    let fi = File::open(smiles_file).expect("Could not open file");

    // open binary file for index
    let mut offset = 0;
    let fo = File::create(smiles_file.to_owned() + ".fp");
    let mut writer = BufWriter::new(fo.unwrap());

    for line in std::io::BufReader::new(fi).lines() {
        let line = line.unwrap();
        let molecule = parse_molecule(&line).unwrap().1;
        let index_item = IndexItem {
            position: offset,
            fingerprint: molecule.ecfp(2, 512),
        };
        offset += line.len() + LINE_ENDING_LENGTH;

        let mut buf = vec![0u8; FINGERPRINT_SIZE / 8 + 8];

        encode_into_slice(&index_item, buf.as_mut_slice(), bincode::config::standard()).unwrap();
        writer.write(&buf).unwrap();
    }
}

#[cfg(test)]
mod test {
    use crate::ringo::ringo::index::index::index;

    #[test]
    fn test_index() {
        index("molecules.smi");
    }
}
