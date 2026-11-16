use clap::Parser;
use ringo::io::smiles::canonical::{canonical_smiles_with, DEFAULT_SEARCH_LIMIT};
use ringo::io::smiles::reader::molecule::{parse_molecule_with, DirectionalStereo};
use ringo::io::smiles::writer::DoubleBondStereoPolicy;
use ringo::io::smiles::TolerateStereo;
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
    /// How to handle invalid, incomplete, or unpreserved `/` and `\` directional marks
    #[arg(long, value_enum, default_value_t = DirectionalStereo::Reject)]
    directional_stereo: DirectionalStereo,
    /// How to handle double-bond stereo configurations that `/` and `\` cannot express
    #[arg(long, value_enum, default_value_t = DoubleBondStereoPolicy::Reject)]
    double_bond_stereo: DoubleBondStereoPolicy,
    /// Tolerate invalid directional marks on input and unrepresentable double-bond stereo on output
    #[arg(long)]
    tolerate_stereo: bool,
}

fn run(args: &Cli) -> Result<(), String> {
    let policy = TolerateStereo::from_flags(
        args.tolerate_stereo,
        args.directional_stereo,
        args.double_bond_stereo,
    );
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
            parse_molecule_with(token, policy.directional)
                .map_err(|_| "invalid or unsupported SMILES".to_owned())
                .and_then(|(_, molecule)| {
                    canonical_smiles_with(&molecule, args.search_limit, policy.double_bond)
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
