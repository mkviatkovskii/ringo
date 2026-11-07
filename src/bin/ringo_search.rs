extern crate ringo;

use ringo::db::search::{similarity_search, SearchMode};

// usage: ringo_search <smiles_file> <query> <min_similarity> <limit> [max_similarity] [best|first]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let max_similarity = args.get(5).map_or(1.0, |s| s.parse().unwrap());
    let mode = match args.get(6).map(String::as_str) {
        None | Some("best") => SearchMode::Best,
        Some("first") => SearchMode::First,
        Some(other) => panic!("unknown mode {other:?}, expected best or first"),
    };
    let results = similarity_search(
        &args[1],
        &args[2],
        args[3].parse().unwrap(),
        max_similarity,
        args[4].parse().unwrap(),
        mode,
    );
    for result in results {
        println!("{:?} {:?}", result.line, result.similarity);
    }
}
