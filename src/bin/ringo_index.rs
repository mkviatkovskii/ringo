use clap::{CommandFactory, Parser};
use ringo::db::index::{index_file_with_options, IndexOptions, DEFAULT_BATCH_LINES};
use ringo::db::index_header::IndexHeader;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "ringo-index",
    version,
    about = "Build an ECFP index from a SMILES file"
)]
struct Cli {
    /// File containing one SMILES record per line
    smiles_file: String,
    /// Fingerprint radius
    #[arg(default_value_t = 2)]
    radius: usize,
    /// Fingerprint length, a positive multiple of 64
    #[arg(default_value_t = 512)]
    bits: usize,
    /// Records processed per batch
    #[arg(default_value_t = DEFAULT_BATCH_LINES,
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..))]
    batch_lines: usize,
    /// Validate each loaded molecule before fingerprinting; skip failed records
    #[arg(long)]
    validate: bool,
}

fn main() -> ExitCode {
    let args = Cli::parse();
    if let Err(error) = IndexHeader::new(args.radius, args.bits) {
        Cli::command()
            .error(clap::error::ErrorKind::ValueValidation, error)
            .exit();
    }
    match index_file_with_options(
        &args.smiles_file,
        IndexOptions {
            radius: args.radius,
            bits: args.bits,
            batch_lines: args.batch_lines,
            validate: args.validate,
        },
    ) {
        Ok(skipped) => {
            if skipped > 0 {
                eprintln!(
                    "skipped {skipped} line(s) that could not be parsed or failed validation"
                );
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("ringo-index: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};

    #[test]
    fn command_arguments() {
        Cli::command().debug_assert();
        let defaults = Cli::try_parse_from(["ringo-index", "m.smi"]).unwrap();
        assert_eq!(
            (
                defaults.radius,
                defaults.bits,
                defaults.batch_lines,
                defaults.validate
            ),
            (2, 512, 65_536, false)
        );
        let full = Cli::try_parse_from(["ringo-index", "--validate", "m.smi", "3", "1024", "1000"])
            .unwrap();
        assert_eq!(
            (full.radius, full.bits, full.batch_lines, full.validate),
            (3, 1024, 1000, true)
        );
        assert!(
            Cli::try_parse_from(["ringo-index", "m.smi", "--validate"])
                .unwrap()
                .validate
        );
        for argv in [
            vec!["ringo-index"],
            vec!["ringo-index", "m.smi", "two"],
            vec!["ringo-index", "m.smi", "2", "512", "0"],
            vec!["ringo-index", "m.smi", "2", "512", "1000", "extra"],
        ] {
            assert!(Cli::try_parse_from(argv).is_err());
        }
    }
}
