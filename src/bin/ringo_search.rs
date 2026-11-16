use clap::Parser;
use ringo::db::exact_search::exact_search_with;
use ringo::db::search::{similarity_search_with, SearchMode};
use ringo::io::smiles::reader::molecule::DirectionalStereo;
use ringo::io::smiles::writer::DoubleBondStereoPolicy;
use ringo::io::smiles::TolerateStereo;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "ringo-search",
    version,
    about = "Search a SMILES fingerprint index"
)]
struct Cli {
    /// Find the records with exactly the query structure instead of similar ones
    #[arg(long, conflicts_with_all = ["min_similarity", "limit", "validate"])]
    exact: bool,
    /// Indexed SMILES file
    smiles_file: String,
    /// Query molecule in SMILES notation
    query: String,
    /// Minimum similarity, inclusive
    #[arg(required_unless_present = "exact")]
    min_similarity: Option<f32>,
    /// Maximum number of results
    #[arg(required_unless_present = "exact")]
    limit: Option<usize>,
    /// Maximum similarity, inclusive
    #[arg(default_value_t = 1.0)]
    max_similarity: f32,
    /// Return the best matches or the first matches in index order
    #[arg(value_enum, default_value = "best")]
    mode: SearchMode,
    /// Validate the loaded query before fingerprinting
    #[arg(long)]
    validate: bool,
    /// How to handle invalid, incomplete, or unpreserved `/` and `\` directional marks
    #[arg(long, value_enum, default_value_t = DirectionalStereo::Reject)]
    directional_stereo: DirectionalStereo,
    /// Tolerate stereo issues on input (`--directional-stereo tolerate`)
    #[arg(long)]
    tolerate_stereo: bool,
}

fn main() -> ExitCode {
    let args = Cli::parse();
    let policy = TolerateStereo::from_flags(
        args.tolerate_stereo,
        args.directional_stereo,
        DoubleBondStereoPolicy::Reject,
    );
    let (Some(min_similarity), Some(limit)) = (args.min_similarity, args.limit) else {
        // clap requires both unless --exact is given
        return search_exact(&args.smiles_file, &args.query, policy.directional);
    };
    let results = match similarity_search_with(
        &args.smiles_file,
        &args.query,
        min_similarity,
        args.max_similarity,
        limit,
        args.mode,
        args.validate,
        policy.directional,
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

/// Prints the lines that hold `query`, and on stderr the lines it could not be compared with.
fn search_exact(smiles_file: &str, query: &str, directional_stereo: DirectionalStereo) -> ExitCode {
    let hits = match exact_search_with(smiles_file, query, directional_stereo) {
        Ok(hits) => hits,
        Err(e) => {
            eprintln!("ringo-search: {e}");
            return ExitCode::FAILURE;
        }
    };
    for line in hits.lines {
        println!("{}", line.trim_end());
    }
    for (line, error) in hits.undecided {
        eprintln!(
            "ringo-search: could not compare with {}: {error:?}",
            line.trim_end()
        );
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{CommandFactory, Parser};
    use ringo::db::search::SearchMode;
    use ringo::io::smiles::reader::molecule::DirectionalStereo;

    #[test]
    fn command_arguments() {
        Cli::command().debug_assert();
        let defaults = Cli::try_parse_from(["ringo-search", "m.smi", "CCO", "0.5", "10"]).unwrap();
        assert_eq!(
            (
                defaults.exact,
                defaults.min_similarity,
                defaults.limit,
                defaults.max_similarity,
                defaults.mode,
                defaults.validate,
                defaults.directional_stereo,
                defaults.tolerate_stereo,
            ),
            (
                false,
                Some(0.5),
                Some(10),
                1.0,
                SearchMode::Best,
                false,
                DirectionalStereo::Reject,
                false,
            )
        );
        let full = Cli::try_parse_from([
            "ringo-search",
            "m.smi",
            "CCO",
            "--validate",
            "--directional-stereo",
            "tolerate",
            "--tolerate-stereo",
            "0.5",
            "10",
            "0.9",
            "first",
        ])
        .unwrap();
        assert_eq!(
            (
                full.max_similarity,
                full.mode,
                full.validate,
                full.directional_stereo,
                full.tolerate_stereo,
            ),
            (
                0.9,
                SearchMode::First,
                true,
                DirectionalStereo::Tolerate,
                true,
            )
        );
        assert!(
            Cli::try_parse_from(["ringo-search", "m.smi", "CCO", "0.5", "10", "--validate"])
                .unwrap()
                .validate
        );
        let exact = Cli::try_parse_from(["ringo-search", "--exact", "m.smi", "CCO"]).unwrap();
        assert_eq!(
            (exact.exact, exact.min_similarity, exact.limit),
            (true, None, None)
        );
        for argv in [
            vec!["ringo-search", "m.smi"],
            vec!["ringo-search", "m.smi", "CCO", "0.5"],
            vec!["ringo-search", "m.smi", "CCO", "high", "10"],
            vec!["ringo-search", "m.smi", "CCO", "0.5", "10", "1.0", "worst"],
            vec!["ringo-search", "--exact", "m.smi"],
            vec!["ringo-search", "--exact", "m.smi", "CCO", "0.5", "10"],
            vec!["ringo-search", "--exact", "--validate", "m.smi", "CCO"],
        ] {
            assert!(Cli::try_parse_from(argv).is_err());
        }
    }
}
