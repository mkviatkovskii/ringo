extern crate ringo;

use ringo::db::index::index_file;
use ringo::db::index_header::IndexHeader;
use std::process::ExitCode;

const USAGE: &str = "usage: ringo-index <smiles_file> [radius=2] [bits=512]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((smiles_file, radius, bits)) = parse_args(&args) else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };
    // index_file panics on invalid parameters, so check them here first
    if let Err(e) = IndexHeader::new(radius, bits) {
        eprintln!("{e}\n{USAGE}");
        return ExitCode::from(2);
    }
    let skipped = index_file(smiles_file, radius, bits);
    if skipped > 0 {
        eprintln!("skipped {skipped} line(s) that could not be parsed");
    }
    ExitCode::SUCCESS
}

fn parse_args(args: &[String]) -> Option<(&str, usize, usize)> {
    if args.is_empty() || args.len() > 3 {
        return None;
    }
    let radius = args.get(1).map_or(Some(2), |s| s.parse().ok())?;
    let bits = args.get(2).map_or(Some(512), |s| s.parse().ok())?;
    Some((&args[0], radius, bits))
}

#[cfg(test)]
mod tests {
    use super::parse_args;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_parse_args() {
        assert_eq!(parse_args(&args(&["m.smi"])), Some(("m.smi", 2, 512)));
        assert_eq!(
            parse_args(&args(&["m.smi", "3", "1024"])),
            Some(("m.smi", 3, 1024))
        );
        assert_eq!(parse_args(&args(&[])), None);
        assert_eq!(parse_args(&args(&["m.smi", "two"])), None);
        assert_eq!(parse_args(&args(&["m.smi", "2", "512", "extra"])), None);
    }
}
