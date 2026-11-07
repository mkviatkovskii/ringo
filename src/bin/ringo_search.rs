use ringo::db::search::{similarity_search, SearchMode};
use std::process::ExitCode;

const USAGE: &str = "usage: ringo-search <smiles_file> <query> <min_similarity> <limit> \
                     [max_similarity=1.0] [best|first]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(args) = parse_args(&args) else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };
    let results = match similarity_search(
        args.smiles_file,
        args.query,
        args.min_similarity,
        args.max_similarity,
        args.limit,
        args.mode,
    ) {
        Ok(results) => results,
        Err(e) => {
            eprintln!("ringo-search: {e}");
            return ExitCode::FAILURE;
        }
    };
    for result in results {
        println!("{:.3}\t{}", result.similarity, result.line.trim_end());
    }
    ExitCode::SUCCESS
}

#[derive(Debug, PartialEq)]
struct Args<'a> {
    smiles_file: &'a str,
    query: &'a str,
    min_similarity: f32,
    limit: usize,
    max_similarity: f32,
    mode: SearchMode,
}

fn parse_args(args: &[String]) -> Option<Args<'_>> {
    if !(4..=6).contains(&args.len()) {
        return None;
    }
    Some(Args {
        smiles_file: &args[0],
        query: &args[1],
        min_similarity: args[2].parse().ok()?,
        limit: args[3].parse().ok()?,
        max_similarity: args.get(4).map_or(Some(1.0), |s| s.parse().ok())?,
        mode: match args.get(5).map(String::as_str) {
            None | Some("best") => SearchMode::Best,
            Some("first") => SearchMode::First,
            Some(_) => return None,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::{parse_args, Args};
    use ringo::db::search::SearchMode;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(std::string::ToString::to_string).collect()
    }

    #[test]
    fn test_parse_args() {
        assert_eq!(
            parse_args(&args(&["m.smi", "CCO", "0.5", "10"])),
            Some(Args {
                smiles_file: "m.smi",
                query: "CCO",
                min_similarity: 0.5,
                limit: 10,
                max_similarity: 1.0,
                mode: SearchMode::Best,
            })
        );
        let full = args(&["m.smi", "CCO", "0.5", "10", "0.9", "first"]);
        let parsed = parse_args(&full).unwrap();
        assert_eq!(
            (parsed.max_similarity, parsed.mode),
            (0.9, SearchMode::First)
        );

        // the panics from issue #20: too few arguments, a bad number, a bad mode
        assert_eq!(parse_args(&args(&["m.smi"])), None);
        assert_eq!(parse_args(&args(&["m.smi", "CCO", "high", "10"])), None);
        assert_eq!(
            parse_args(&args(&["m.smi", "CCO", "0.5", "10", "1.0", "worst"])),
            None
        );
    }
}
