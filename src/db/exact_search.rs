//! `<smiles_file>.key`: structure keys for exact search (#58).
//!
//! Magic (8 bytes), format version (`u16`), 6 zero bytes, SMILES file length (`u64`),
//! record count (`u64`), all little-endian. Then one record per parsed line: its
//! `Molecule::exact_key` and its line offset (`u64` each), sorted by key, then offset.
use crate::db::with_path;
use crate::io::smiles::reader::molecule::parse_molecule;
use crate::model::exact_match::ExactMatchError;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Seek, Write};

const MAGIC: [u8; 8] = *b"RINGOKEY";
// Bump once per release when the key or the layout changes.
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 32;
const RECORD_SIZE: usize = 16;

/// Writes `path` with the `(key, line offset)` pairs of a SMILES file `smiles_len` bytes long,
/// through a temporary file renamed over it at the end.
pub(crate) fn write_key_index(
    path: &str,
    smiles_len: u64,
    mut records: Vec<(u64, u64)>,
) -> io::Result<()> {
    records.sort_unstable();
    let tmp_path = path.to_owned() + ".tmp";
    let mut out = io::BufWriter::new(File::create(&tmp_path).map_err(with_path(&tmp_path))?);
    out.write_all(&MAGIC)?;
    out.write_all(&VERSION.to_le_bytes())?;
    out.write_all(&[0; 6])?;
    out.write_all(&smiles_len.to_le_bytes())?;
    out.write_all(&(records.len() as u64).to_le_bytes())?;
    for (key, offset) in records {
        out.write_all(&key.to_le_bytes())?;
        out.write_all(&offset.to_le_bytes())?;
    }
    let out = out.into_inner().map_err(io::IntoInnerError::into_error)?;
    // sync before the rename
    out.sync_all()?;
    std::fs::rename(&tmp_path, path).map_err(with_path(path))
}

pub struct ExactHits {
    /// The lines that hold the query, in file order.
    pub lines: Vec<String>,
    /// Lines with the query's key that could not be compared, and why.
    pub undecided: Vec<(String, ExactMatchError)>,
}

/// Finds the lines of `smiles_file` that hold the same molecule as `query`, by its key in
/// `<smiles_file>.key` and then `Molecule::exact_match` on each line with that key.
pub fn exact_search(smiles_file: &str, query: &str) -> io::Result<ExactHits> {
    let query = parse_molecule(query)
        .map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("could not parse query {query}"),
            )
        })?
        .1;
    let mut reader = BufReader::new(File::open(smiles_file).map_err(with_path(smiles_file))?);
    let key_path = smiles_file.to_owned() + ".key";
    let file = File::open(&key_path).map_err(with_path(&key_path))?;
    #[expect(
        unsafe_code,
        reason = "ringo-index never writes an index in place: it renames a finished file over \
                  it, so the mapped bytes stay fixed"
    )]
    // SAFETY: as above, nothing in ringo changes the mapped file while it is mapped
    let map = unsafe { Mmap::map(&file) }.map_err(with_path(&key_path))?;

    let stale = |msg: String| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{key_path}: {msg}: re-index"),
        )
    };
    let u64_at = |i: usize| {
        let mut word = [0u8; 8];
        word.copy_from_slice(&map[i..i + 8]);
        u64::from_le_bytes(word)
    };
    if map.len() < HEADER_SIZE || map[..8] != MAGIC {
        return Err(stale("not a ringo key index".into()));
    }
    let version = u16::from_le_bytes([map[8], map[9]]);
    if version != VERSION {
        return Err(stale(format!(
            "format version {version} is not supported, expected {VERSION}"
        )));
    }
    // shortcut: compares lengths only, like the similarity index
    let smiles_len = reader
        .get_ref()
        .metadata()
        .map_err(with_path(smiles_file))?
        .len();
    if u64_at(16) != smiles_len {
        return Err(stale(format!("{smiles_file} changed since it was indexed")));
    }
    let count = u64_at(24);
    if count.checked_mul(RECORD_SIZE as u64) != Some((map.len() - HEADER_SIZE) as u64) {
        return Err(stale(format!(
            "does not hold the {count} records its header lists"
        )));
    }

    // (key, offset) word pairs
    let records = map[HEADER_SIZE..].as_chunks::<8>().0.as_chunks::<2>().0;
    let key = query.exact_key();
    let start = records.partition_point(|[k, _]| u64::from_le_bytes(*k) < key);
    let mut hits = ExactHits {
        lines: Vec::new(),
        undecided: Vec::new(),
    };
    for [_, offset] in records[start..]
        .iter()
        .take_while(|[k, _]| u64::from_le_bytes(*k) == key)
    {
        reader.seek(io::SeekFrom::Start(u64::from_le_bytes(*offset)))?;
        let mut line = String::new();
        reader
            .read_line(&mut line)
            .map_err(with_path(smiles_file))?;
        let (_, molecule) = parse_molecule(line.trim_end_matches(['\n', '\r']))
            .map_err(|_| stale(format!("{smiles_file} changed since it was indexed")))?;
        match query.exact_match(&molecule) {
            Ok(true) => hits.lines.push(line),
            Ok(false) => {}
            Err(error) => hits.undecided.push((line, error)),
        }
    }
    Ok(hits)
}

#[cfg(test)]
mod test {
    use crate::db::exact_search::exact_search;
    use crate::db::index::test::TestDir;
    use crate::db::index::{index_file, DEFAULT_BATCH_LINES};
    use crate::model::exact_match::ExactMatchError;

    fn index(dir: &TestDir, lines: &[&str]) -> String {
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, lines.join("\n") + "\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap().to_owned();
        index_file(&smiles_file, 2, 512, DEFAULT_BATCH_LINES, false).unwrap();
        smiles_file
    }

    fn found(smiles_file: &str, query: &str) -> Vec<String> {
        let hits = exact_search(smiles_file, query).unwrap();
        assert_eq!(hits.undecided, []);
        hits.lines
            .iter()
            .map(|line| line.trim_end().to_owned())
            .collect()
    }

    #[test]
    fn test_exact_search_finds_every_spelling_in_line_order() {
        let dir = TestDir::new("test_exact_search");
        let smiles_file = index(
            &dir,
            &[
                "OCC ethanol",
                "C[C@@H](C(=O)O)N L-alanine",
                "COC dimethyl ether",
                "C(",
                "C[C@H](C(=O)O)N D-alanine",
                "CC(N)C(=O)O alanine",
                "CCO ethanol again",
                "c1ccccc1 benzene",
            ],
        );
        assert_eq!(
            found(&smiles_file, "CCO"),
            ["OCC ethanol", "CCO ethanol again"]
        );
        assert_eq!(
            found(&smiles_file, "N[C@@H](C)C(=O)O"),
            ["C[C@@H](C(=O)O)N L-alanine"]
        );
        assert_eq!(found(&smiles_file, "NC(C)C(=O)O"), ["CC(N)C(=O)O alanine"]);
        assert_eq!(found(&smiles_file, "C1=CC=CC=C1"), ["c1ccccc1 benzene"]);
        assert_eq!(found(&smiles_file, "CCN"), Vec::<String>::new());
    }

    #[test]
    fn test_exact_search_reports_lines_it_cannot_compare() {
        let dir = TestDir::new("test_exact_search_undecided");
        let smiles_file = index(&dir, &["C[Pt@SP1](F)(Cl)N", "CCO"]);
        let hits = exact_search(&smiles_file, "C[Pt@SP2](F)(Cl)N").unwrap();
        assert_eq!(hits.lines, Vec::<String>::new());
        let undecided: Vec<_> = hits
            .undecided
            .iter()
            .map(|(line, error)| (line.trim_end(), *error))
            .collect();
        assert_eq!(
            undecided,
            [("C[Pt@SP1](F)(Cl)N", ExactMatchError::UnsupportedStereo)]
        );
    }

    #[test]
    fn test_exact_search_rejects_stale_or_missing_index() {
        let dir = TestDir::new("test_exact_search_stale");
        let smiles_file = index(&dir, &["CCO", "CCC"]);
        let key_file = smiles_file.clone() + ".key";
        let key = std::fs::read(&key_file).unwrap();

        std::fs::write(&smiles_file, "CCCO\nCCC\n").unwrap();
        let e = exact_search(&smiles_file, "CCC").err().unwrap();
        assert!(
            e.to_string().contains("changed since it was indexed"),
            "{e}"
        );

        std::fs::write(&smiles_file, "CCO\nCCC\n").unwrap();
        for bad in [&key[..key.len() - 1], &key[..4]] {
            std::fs::write(&key_file, bad).unwrap();
            let e = exact_search(&smiles_file, "CCC").err().unwrap();
            assert!(e.to_string().contains("re-index"), "{e}");
        }
        std::fs::write(&key_file, &key).unwrap();
        assert_eq!(found(&smiles_file, "CCC"), ["CCC"]);

        std::fs::remove_file(&key_file).unwrap();
        let e = exact_search(&smiles_file, "CCC").err().unwrap();
        assert_eq!(e.kind(), std::io::ErrorKind::NotFound);
        assert_eq!(
            exact_search(&smiles_file, "C1CC").err().unwrap().kind(),
            std::io::ErrorKind::InvalidInput
        );
    }
}
