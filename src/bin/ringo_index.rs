extern crate ringo;

use ringo::db::index::index_file;

// usage: ringo_index <smiles_file> [radius] [bits]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let radius = args.get(2).map_or(2, |s| s.parse().unwrap());
    let bits = args.get(3).map_or(512, |s| s.parse().unwrap());
    let skipped = index_file(&args[1], radius, bits);
    if skipped > 0 {
        eprintln!("skipped {skipped} line(s) that could not be parsed");
    }
}
