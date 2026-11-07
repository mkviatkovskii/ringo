use crate::db::index_header::IndexHeader;
use crate::db::index_item::IndexItem;
use crate::io::smiles::reader::molecule::parse_molecule;
use crate::model::fingerprint::Fingerprint;
use fixedbitset::FixedBitSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Writes `<smiles_file>.fp` with ECFP fingerprints of the given radius and length
/// (a positive multiple of 64), and returns how many lines could not be parsed.
/// Those lines get an empty fingerprint, so records stay aligned with lines.
pub fn index_file(smiles_file: &str, radius: usize, bits: usize) -> usize {
    let header = IndexHeader::new(radius, bits).unwrap_or_else(|e| panic!("{e}"));

    // open file for reading
    let fi = File::open(smiles_file).expect("Could not open file");

    // write to a temp file and rename at the end, so a failed run keeps the old index
    let fp_path = smiles_file.to_owned() + ".fp";
    let tmp_path = fp_path.clone() + ".tmp";
    let mut offset = 0;
    let mut skipped = 0;
    let mut writer = BufWriter::new(File::create(&tmp_path).unwrap());
    writer.write_all(&header.to_bytes()).unwrap();

    let mut reader = BufReader::new(fi);
    let mut line = String::new();
    loop {
        line.clear();
        // read_line returns the bytes consumed, including "\n" or "\r\n"
        let read = reader.read_line(&mut line).unwrap();
        if read == 0 {
            break;
        }
        let fingerprint = match parse_molecule(line.trim_end_matches(['\n', '\r'])) {
            Ok((_, molecule)) => molecule.ecfp(radius, bits),
            Err(_) => {
                skipped += 1;
                Fingerprint(FixedBitSet::with_capacity(bits))
            }
        };
        let index_item = IndexItem {
            position: offset,
            fingerprint,
        };
        offset += read;

        writer.write_all(&index_item.to_bytes()).unwrap();
    }
    // flush, sync and close before rename (an open file cannot be renamed on Windows)
    writer.into_inner().unwrap().sync_all().unwrap();
    std::fs::rename(&tmp_path, &fp_path).unwrap();
    skipped
}

#[cfg(test)]
pub(crate) mod test {
    use crate::db::index::index_file;
    use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};
    use crate::db::index_item::index_item_size;
    use std::path::{Path, PathBuf};

    /// Per-test scratch directory under the system temp dir, removed on drop
    /// (also when the test panics).
    pub(crate) struct TestDir(PathBuf);

    impl TestDir {
        pub(crate) fn new(name: &str) -> TestDir {
            let path = std::env::temp_dir().join(format!("ringo-{}-{}", name, std::process::id()));
            let _ = std::fs::remove_dir_all(&path);
            std::fs::create_dir_all(&path).unwrap();
            TestDir(path)
        }

        pub(crate) fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Copies the sample `molecules.smi` into `dir` and returns the copy's path,
    /// so tests never write index files into the repository.
    pub(crate) fn copy_sample_molecules(dir: &Path) -> PathBuf {
        let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("molecules.smi");
        let target = dir.join("molecules.smi");
        std::fs::copy(source, &target).unwrap();
        target
    }

    #[test]
    fn test_index() {
        let dir = TestDir::new("test_index");
        let smiles_file = copy_sample_molecules(dir.path());
        index_file(smiles_file.to_str().unwrap(), 3, 1024);

        let fp = std::fs::read(dir.path().join("molecules.smi.fp")).unwrap();
        // header, then 2 molecules, one fixed-size record each
        assert_eq!(fp.len(), INDEX_HEADER_SIZE + 2 * index_item_size(1024));
        let header = IndexHeader::from_bytes(fp[..INDEX_HEADER_SIZE].try_into().unwrap());
        assert_eq!(header, IndexHeader::new(3, 1024));
    }

    #[test]
    fn test_index_skips_unparseable_lines() {
        let dir = TestDir::new("test_index_skips_unparseable_lines");
        let smiles_file = dir.path().join("m.smi");
        // the chiral line is not supported yet and must not abort indexing
        std::fs::write(&smiles_file, "CCO\nN[C@H](C)C(=O)O\nCCC\n").unwrap();

        let skipped = index_file(smiles_file.to_str().unwrap(), 2, 512);

        assert_eq!(skipped, 1);
        let fp_file = dir.path().join("m.smi.fp");
        assert_eq!(
            std::fs::metadata(fp_file).unwrap().len(),
            (INDEX_HEADER_SIZE + 3 * index_item_size(512)) as u64
        );
        assert!(!dir.path().join("m.smi.fp.tmp").exists());
    }
}
