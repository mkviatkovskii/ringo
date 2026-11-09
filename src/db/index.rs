use crate::db::index_header::IndexHeader;
use crate::db::index_item::{index_item_size, split_record, IndexItem};
use crate::db::with_path;
use crate::io::smiles::reader::molecule::parse_molecule;
use crate::model::fingerprint::Fingerprint;
use fixedbitset::FixedBitSet;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, Write};

/// Lines fingerprinted in parallel per batch: enough to keep every core busy, a few MB of memory
pub const DEFAULT_BATCH_LINES: usize = 65_536;

#[derive(Clone, Copy, Debug)]
pub struct IndexOptions {
    pub radius: usize,
    pub bits: usize,
    pub batch_lines: usize,
    /// Run `Molecule::validate()` after loading each record, before fingerprinting.
    pub validate: bool,
}

/// Writes `<smiles_file>.fp` with ECFP fingerprints of the given radius and length
/// (a positive multiple of 64), and returns how many lines could not be parsed.
/// Those lines get an empty fingerprint, so records stay aligned with lines.
/// Lines are fingerprinted in parallel `batch_lines` at a time. Records are grouped by bit
/// count, so search can skip the groups that cannot reach its threshold.
pub fn index_file(
    smiles_file: &str,
    radius: usize,
    bits: usize,
    batch_lines: usize,
) -> io::Result<usize> {
    index_file_with_options(
        smiles_file,
        IndexOptions {
            radius,
            bits,
            batch_lines,
            validate: false,
        },
    )
}

/// Builds an index with explicit options. Loading and optional validation failures
/// are counted as skipped records and stored as empty fingerprints.
pub fn index_file_with_options(smiles_file: &str, options: IndexOptions) -> io::Result<usize> {
    let IndexOptions {
        radius,
        bits,
        batch_lines,
        validate,
    } = options;
    let header = IndexHeader::new(radius, bits)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    if batch_lines == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "batch_lines must be positive",
        ));
    }

    // open file for reading
    let fi = File::open(smiles_file).map_err(with_path(smiles_file))?;

    // write to temp files and rename at the end, so a failed run keeps the old index; records
    // go to the first one in line order, since grouping them needs the size of every group
    let fp_path = smiles_file.to_owned() + ".fp";
    let unsorted_path = fp_path.clone() + ".unsorted.tmp";
    let tmp_path = fp_path.clone() + ".tmp";
    let mut offset = 0;
    let mut skipped = 0;
    let mut count = 0;
    // records per bit count
    let mut groups = vec![0u64; bits + 1];
    let mut writer =
        BufWriter::new(File::create(&unsorted_path).map_err(with_path(&unsorted_path))?);

    let mut reader = BufReader::new(fi);
    let mut lines =
        read_batch(&mut reader, &mut offset, batch_lines).map_err(with_path(smiles_file))?;
    let mut fingerprinted: Vec<(usize, Option<Fingerprint>)> = Vec::new();
    // one thread writes the previous batch and reads the next while the others fingerprint this
    // one (about 10% faster on ChEMBL than taking turns); the loop ends after writing the last batch
    loop {
        let (next, fingerprints) = rayon::join(
            || -> io::Result<_> {
                for (position, fingerprint) in fingerprinted.drain(..) {
                    let fingerprint = fingerprint.unwrap_or_else(|| {
                        skipped += 1;
                        Fingerprint(FixedBitSet::with_capacity(bits))
                    });
                    groups[fingerprint.0.count_ones(..)] += 1;
                    writer.write_all(
                        &IndexItem {
                            position,
                            fingerprint,
                        }
                        .to_bytes(),
                    )?;
                    count += 1;
                }
                read_batch(&mut reader, &mut offset, batch_lines).map_err(with_path(smiles_file))
            },
            || {
                lines
                    .par_iter()
                    .map(|(_, line)| {
                        parse_molecule(line.trim_end_matches(['\n', '\r']))
                            .ok()
                            .filter(|(_, molecule)| !validate || molecule.validate().is_ok())
                            .map(|(_, molecule)| molecule.ecfp(radius, bits))
                    })
                    .collect::<Vec<_>>()
            },
        );
        let next = next?;
        if lines.is_empty() {
            break;
        }
        fingerprinted = lines
            .into_iter()
            .map(|(position, _)| position)
            .zip(fingerprints)
            .collect();
        lines = next;
    }
    // flush and close
    writer
        .into_inner()
        .map_err(io::IntoInnerError::into_error)?;
    let min_ones = groups.iter().position(|&n| n > 0).unwrap_or(0);
    let max_ones = groups.iter().rposition(|&n| n > 0).unwrap_or(0);
    let header = IndexHeader {
        smiles_len: offset as u64,
        count,
        min_ones,
        max_ones,
        ..header
    };
    write_grouped(
        &unsorted_path,
        &tmp_path,
        &header,
        &groups[min_ones..=max_ones],
    )?;
    std::fs::rename(&tmp_path, &fp_path).map_err(with_path(&fp_path))?;
    std::fs::remove_file(&unsorted_path).map_err(with_path(&unsorted_path))?;
    Ok(skipped)
}

/// Bytes of a group buffered before they are written at the group's place in the index
const GROUP_BUFFER: usize = 1 << 16;

/// Writes `header`, the group table and the records of `unsorted_path` to `path`, grouped
/// by bit count and in the same order within a group. `groups` holds the record count of
/// each bit count from `header.min_ones` to `header.max_ones`.
fn write_grouped(
    unsorted_path: &str,
    path: &str,
    header: &IndexHeader,
    groups: &[u64],
) -> io::Result<()> {
    let item_size = index_item_size(header.bits);
    let mut out = File::create(path).map_err(with_path(path))?;
    out.write_all(&header.to_bytes())?;
    out.write_all(
        &groups
            .iter()
            .flat_map(|n| n.to_le_bytes())
            .collect::<Vec<_>>(),
    )?;
    // where the next record of each group goes
    let mut next: Vec<u64> = groups
        .iter()
        .scan(header.records_start() as u64, |start, n| {
            let group_start = *start;
            *start += n * item_size as u64;
            Some(group_start)
        })
        .collect();
    let mut flush = |group: usize, buffer: &mut Vec<u8>| -> io::Result<()> {
        out.seek(io::SeekFrom::Start(next[group]))?;
        out.write_all(buffer)?;
        next[group] += buffer.len() as u64;
        buffer.clear();
        Ok(())
    };

    let mut buffers = vec![Vec::new(); groups.len()];
    let mut input = BufReader::with_capacity(
        1 << 20,
        File::open(unsorted_path).map_err(with_path(unsorted_path))?,
    );
    let mut record = vec![0; item_size];
    for _ in 0..header.count {
        input
            .read_exact(&mut record)
            .map_err(with_path(unsorted_path))?;
        let ones: usize = split_record(&record)?
            .1
            .as_chunks::<8>()
            .0
            .iter()
            .map(|word| u64::from_le_bytes(*word).count_ones() as usize)
            .sum();
        let group = ones.wrapping_sub(header.min_ones);
        let buffer = buffers.get_mut(group).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{unsorted_path} changed while indexing"),
            )
        })?;
        buffer.extend_from_slice(&record);
        if buffer.len() >= GROUP_BUFFER {
            flush(group, buffer)?;
        }
    }
    for (group, buffer) in buffers.iter_mut().enumerate() {
        flush(group, buffer)?;
    }
    // sync before the caller renames it
    out.sync_all()
}

/// Reads up to `batch_lines` lines with their byte offsets, moving `offset` past them.
fn read_batch(
    reader: &mut impl BufRead,
    offset: &mut usize,
    batch_lines: usize,
) -> io::Result<Vec<(usize, String)>> {
    let mut lines = Vec::with_capacity(batch_lines);
    while lines.len() < batch_lines {
        let mut line = String::new();
        // read_line returns the bytes consumed, including "\n" or "\r\n"
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        lines.push((*offset, line));
        *offset += read;
    }
    Ok(lines)
}

#[cfg(test)]
pub(crate) mod test {
    use crate::db::index::{index_file, DEFAULT_BATCH_LINES};
    use crate::db::index_header::{IndexHeader, INDEX_HEADER_SIZE};
    use crate::db::index_item::{index_item_size, split_record};
    use crate::io::smiles::reader::molecule::parse_molecule;
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
        let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/molecules.smi");
        let target = dir.join("molecules.smi");
        std::fs::copy(source, &target).unwrap();
        target
    }

    #[test]
    fn test_index() {
        let dir = TestDir::new("test_index");
        let smiles_file = copy_sample_molecules(dir.path());
        index_file(smiles_file.to_str().unwrap(), 3, 1024, DEFAULT_BATCH_LINES).unwrap();

        let fp = std::fs::read(dir.path().join("molecules.smi.fp")).unwrap();
        let header = IndexHeader::from_bytes(fp[..INDEX_HEADER_SIZE].try_into().unwrap()).unwrap();
        let smiles = std::fs::read_to_string(&smiles_file).unwrap();
        let ones: Vec<usize> = smiles
            .lines()
            .map(|line| {
                parse_molecule(line)
                    .unwrap()
                    .1
                    .ecfp(3, 1024)
                    .0
                    .count_ones(..)
            })
            .collect();
        assert_eq!(
            header,
            IndexHeader {
                smiles_len: smiles.len() as u64,
                count: 2,
                min_ones: *ones.iter().min().unwrap(),
                max_ones: *ones.iter().max().unwrap(),
                ..IndexHeader::new(3, 1024).unwrap()
            }
        );
        // header, group table, then 2 molecules, one fixed-size record each
        assert_eq!(fp.len(), header.records_start() + 2 * index_item_size(1024));
    }

    #[test]
    fn test_index_skips_unparseable_lines() {
        let dir = TestDir::new("test_index_skips_unparseable_lines");
        let smiles_file = dir.path().join("m.smi");
        // Unsupported syntax and bond-order overflow must not abort indexing.
        let overflow = format!("C{}", "(C)".repeat(300));
        std::fs::write(&smiles_file, format!("CCO\nC~C\n{overflow}\nCCC\n")).unwrap();

        let skipped =
            index_file(smiles_file.to_str().unwrap(), 2, 512, DEFAULT_BATCH_LINES).unwrap();

        assert_eq!(skipped, 2);
        let fp = std::fs::read(dir.path().join("m.smi.fp")).unwrap();
        let header = IndexHeader::from_bytes(fp[..INDEX_HEADER_SIZE].try_into().unwrap()).unwrap();
        // the unparsed lines have no bits set
        assert_eq!(header.min_ones, 0);
        assert_eq!(fp.len(), header.records_start() + 4 * index_item_size(512));
        assert!(!dir.path().join("m.smi.fp.tmp").exists());
        assert!(!dir.path().join("m.smi.fp.unsorted.tmp").exists());
    }

    #[test]
    fn test_index_groups_records_by_bit_count() {
        let dir = TestDir::new("test_index_groups_records");
        let smiles_file = dir.path().join("m.smi");
        let cycle = ["CCO", "c1ccccc1 benzene", "C(", "CCN"];
        // batches of 999 split the cycle at different points, and each group holds more
        // than one buffer of records
        let lines: Vec<&str> = cycle.iter().copied().cycle().take(4001).collect();
        std::fs::write(&smiles_file, lines.join("\n")).unwrap();

        let skipped = index_file(smiles_file.to_str().unwrap(), 2, 512, 999).unwrap();

        assert_eq!(skipped, lines.iter().filter(|line| **line == "C(").count());
        // line offsets and fingerprints, the empty one for the line that does not parse,
        // then sorted by bit count; the sort is stable, so line order stays within a group
        let mut offset = 0;
        let mut expected: Vec<(usize, Vec<u8>)> = Vec::new();
        for line in &lines {
            let fingerprint = match parse_molecule(line) {
                Ok((_, molecule)) => molecule.ecfp(2, 512).to_bytes(),
                Err(_) => vec![0; 64],
            };
            expected.push((offset, fingerprint));
            offset += line.len() + 1;
        }
        let ones = |fingerprint: &[u8]| -> usize {
            fingerprint
                .iter()
                .map(|byte| byte.count_ones() as usize)
                .sum()
        };
        expected.sort_by_key(|(_, fingerprint)| ones(fingerprint));

        let fp = std::fs::read(dir.path().join("m.smi.fp")).unwrap();
        let header = IndexHeader::from_bytes(fp[..INDEX_HEADER_SIZE].try_into().unwrap()).unwrap();
        let records: Vec<(usize, Vec<u8>)> = fp[header.records_start()..]
            .chunks_exact(index_item_size(512))
            .map(|record| {
                let (position, fingerprint) = split_record(record).unwrap();
                (position, fingerprint.to_vec())
            })
            .collect();
        assert_eq!(records, expected);

        assert_eq!(header.min_ones, ones(&expected[0].1));
        assert_eq!(header.max_ones, ones(&expected[expected.len() - 1].1));
        let table: Vec<u64> = fp[INDEX_HEADER_SIZE..header.records_start()]
            .as_chunks::<8>()
            .0
            .iter()
            .map(|count| u64::from_le_bytes(*count))
            .collect();
        let mut counts = vec![0; header.groups()];
        for (_, fingerprint) in &expected {
            counts[ones(fingerprint) - header.min_ones] += 1;
        }
        assert_eq!(table, counts);
    }

    #[test]
    fn test_index_errors_instead_of_panicking() {
        let dir = TestDir::new("test_index_errors");
        let missing = dir.path().join("missing.smi");
        let e = index_file(missing.to_str().unwrap(), 2, 512, DEFAULT_BATCH_LINES).unwrap_err();
        assert_eq!(e.kind(), std::io::ErrorKind::NotFound);
        assert!(e.to_string().contains("missing.smi"), "{e}");

        let smiles_file = copy_sample_molecules(dir.path());
        let e = index_file(smiles_file.to_str().unwrap(), 2, 100, DEFAULT_BATCH_LINES).unwrap_err();
        assert_eq!(e.kind(), std::io::ErrorKind::InvalidInput);
        // an empty batch would end indexing at once and write an empty index
        let e = index_file(smiles_file.to_str().unwrap(), 2, 512, 0).unwrap_err();
        assert_eq!(e.kind(), std::io::ErrorKind::InvalidInput);
    }

    #[test]
    fn test_failed_rebuild_keeps_old_index() {
        let dir = TestDir::new("test_failed_rebuild");
        let smiles_file = dir.path().join("m.smi");
        std::fs::write(&smiles_file, "CCO\nCCC\n").unwrap();
        let smiles_file = smiles_file.to_str().unwrap();
        index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap();
        let fp_file = dir.path().join("m.smi.fp");
        let old = std::fs::read(&fp_file).unwrap();

        // the second line is not UTF-8, so the rebuild fails after it has started writing
        std::fs::write(smiles_file, b"CCO\n\xff\xfe\n").unwrap();
        let e = index_file(smiles_file, 2, 512, DEFAULT_BATCH_LINES).unwrap_err();
        assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
        assert!(dir.path().join("m.smi.fp.unsorted.tmp").exists());
        assert_eq!(std::fs::read(&fp_file).unwrap(), old);
    }
}
