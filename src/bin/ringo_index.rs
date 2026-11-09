use ringo::db::index::{index_file, DEFAULT_BATCH_LINES};
use ringo::db::index_header::IndexHeader;
use std::process::ExitCode;

const USAGE: &str = "usage: ringo-index <smiles_file> [radius=2] [bits=512] [batch_lines=65536]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((smiles_file, radius, bits, batch_lines)) = parse_args(&args) else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };
    // bad parameters are a usage error (exit 2), not a failed run
    if let Err(e) = IndexHeader::new(radius, bits) {
        eprintln!("{e}\n{USAGE}");
        return ExitCode::from(2);
    }
    match index_file(smiles_file, radius, bits, batch_lines) {
        Ok(skipped) => {
            if skipped > 0 {
                eprintln!("skipped {skipped} line(s) that could not be parsed");
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("ringo-index: {e}");
            ExitCode::FAILURE
        }
    }
}

fn parse_args(args: &[String]) -> Option<(&str, usize, usize, usize)> {
    if args.is_empty() || args.len() > 4 {
        return None;
    }
    let radius = args.get(1).map_or(Some(2), |s| s.parse().ok())?;
    let bits = args.get(2).map_or(Some(512), |s| s.parse().ok())?;
    let batch_lines = args
        .get(3)
        .map_or(Some(DEFAULT_BATCH_LINES), |s| s.parse().ok())
        .filter(|&n| n > 0)?;
    Some((&args[0], radius, bits, batch_lines))
}

#[cfg(test)]
mod tests {
    use super::parse_args;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(std::string::ToString::to_string).collect()
    }

    #[test]
    fn test_parse_args() {
        assert_eq!(
            parse_args(&args(&["m.smi"])),
            Some(("m.smi", 2, 512, 65_536))
        );
        assert_eq!(
            parse_args(&args(&["m.smi", "3", "1024"])),
            Some(("m.smi", 3, 1024, 65_536))
        );
        assert_eq!(
            parse_args(&args(&["m.smi", "3", "1024", "1000"])),
            Some(("m.smi", 3, 1024, 1000))
        );
        assert_eq!(parse_args(&args(&[])), None);
        assert_eq!(parse_args(&args(&["m.smi", "two"])), None);
        assert_eq!(parse_args(&args(&["m.smi", "2", "512", "0"])), None);
        assert_eq!(
            parse_args(&args(&["m.smi", "2", "512", "1000", "extra"])),
            None
        );
    }
}
