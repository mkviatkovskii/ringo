use crate::db::index_item::IndexItem;
use crate::io::smiles::reader::molecule::parse_molecule;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};

#[cfg(windows)]
const LINE_ENDING_LENGTH: usize = 2;
#[cfg(not(windows))]
const LINE_ENDING_LENGTH: usize = 1;

pub fn index_file(smiles_file: &str) {
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

        writer.write_all(&index_item.to_bytes()).unwrap();
    }
    writer.flush().unwrap();
}

#[cfg(test)]
pub(crate) mod test {
    use crate::db::index::index_file;
    use crate::db::index_item::INDEX_ITEM_SIZE;
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
        index_file(smiles_file.to_str().unwrap());

        let fp_file = dir.path().join("molecules.smi.fp");
        // 2 molecules, one fixed-size record each
        assert_eq!(
            std::fs::metadata(fp_file).unwrap().len(),
            2 * INDEX_ITEM_SIZE as u64
        );
    }
}
