use clap::Parser;
use ringo::db::search::{similarity_search_with_options, SearchMode, SearchOptions};
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "ringo-search",
    version,
    about = "Search a SMILES fingerprint index"
)]
struct Cli {
    /// Indexed SMILES file
    smiles_file: String,
    /// Query molecule in SMILES notation
    query: String,
    /// Minimum similarity, inclusive
    min_similarity: f32,
    /// Maximum number of results
    limit: usize,
    /// Maximum similarity, inclusive
    #[arg(default_value_t = 1.0)]
    max_similarity: f32,
    /// Return the best matches or the first matches in index order
    #[arg(value_enum, default_value = "best")]
    mode: SearchMode,
    /// Validate the loaded query before fingerprinting
    #[arg(long)]
    validate: bool,
}

fn main() -> ExitCode {
    let args = Cli::parse();
    let results = match similarity_search_with_options(
        &args.smiles_file,
        &args.query,
        SearchOptions {
            min_similarity: args.min_similarity,
            max_similarity: args.max_similarity,
            limit: args.limit,
            mode: args.mode,
            validate: args.validate,
        },
    ) {
        Ok(results) => results,
        Err(error) => {
            eprintln!("ringo-search: {error}");
            return ExitCode::FAILURE;
        }
    };
    for result in results {
        println!("{:.3}\t{}", result.similarity, result.line.trim_end());
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};
    use ringo::db::search::SearchMode;

    #[test]
    fn command_arguments() {
        Cli::command().debug_assert();
        let defaults = Cli::try_parse_from(["ringo-search", "m.smi", "CCO", "0.5", "10"]).unwrap();
        assert_eq!(
            (
                defaults.min_similarity,
                defaults.limit,
                defaults.max_similarity,
                defaults.mode,
                defaults.validate
            ),
            (0.5, 10, 1.0, SearchMode::Best, false)
        );
        let full = Cli::try_parse_from([
            "ringo-search",
            "m.smi",
            "CCO",
            "--validate",
            "0.5",
            "10",
            "0.9",
            "first",
        ])
        .unwrap();
        assert_eq!(
            (full.max_similarity, full.mode, full.validate),
            (0.9, SearchMode::First, true)
        );
        assert!(
            Cli::try_parse_from(["ringo-search", "m.smi", "CCO", "0.5", "10", "--validate"])
                .unwrap()
                .validate
        );
        for argv in [
            vec!["ringo-search", "m.smi"],
            vec!["ringo-search", "m.smi", "CCO", "high", "10"],
            vec!["ringo-search", "m.smi", "CCO", "0.5", "10", "1.0", "worst"],
        ] {
            assert!(Cli::try_parse_from(argv).is_err());
        }
    }
}
