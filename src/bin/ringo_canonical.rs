use clap::Parser;
use ringo::io::smiles::canonical::{canonical_smiles_with_limit, DEFAULT_SEARCH_LIMIT};
use ringo::io::smiles::reader::molecule::parse_molecule;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::process::ExitCode;

#[derive(Parser)]
#[command(about = "Write Ringo canonical SMILES, one record per line")]
struct Cli {
    /// Input SMILES file, or - for standard input
    smiles_file: String,
    /// Maximum search nodes per molecule
    #[arg(long, default_value_t = DEFAULT_SEARCH_LIMIT)]
    search_limit: usize,
}

fn run(args: &Cli) -> Result<(), String> {
    let input: Box<dyn BufRead> = if args.smiles_file == "-" {
        Box::new(io::stdin().lock())
    } else {
        Box::new(BufReader::new(
            File::open(&args.smiles_file).map_err(|e| e.to_string())?,
        ))
    };
    let mut output = io::BufWriter::new(io::stdout().lock());
    let mut failed = false;
    for (index, line) in input.lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim();
        let (token, name) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        let result = if token.is_empty() {
            Ok(String::new())
        } else {
            parse_molecule(token)
                .map_err(|_| "invalid or unsupported SMILES".to_owned())
                .and_then(|(_, molecule)| {
                    canonical_smiles_with_limit(&molecule, args.search_limit)
                        .map_err(|e| format!("{e:?}"))
                })
        };
        let smiles = match result {
            Ok(smiles) => smiles,
            Err(error) => {
                eprintln!("ringo-canonical: line {}: {error}", index + 1);
                failed = true;
                String::new()
            }
        };
        write!(output, "{smiles}").map_err(|e| e.to_string())?;
        if !name.trim().is_empty() {
            write!(output, "\t{}", name.trim()).map_err(|e| e.to_string())?;
        }
        writeln!(output).map_err(|e| e.to_string())?;
    }
    output.flush().map_err(|e| e.to_string())?;
    if failed {
        Err("one or more records failed".to_owned())
    } else {
        Ok(())
    }
}

fn main() -> ExitCode {
    match run(&Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ringo-canonical: {error}");
            ExitCode::FAILURE
        }
    }
}
