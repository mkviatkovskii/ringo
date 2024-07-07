extern crate ringo;

use ringo::db::search::similarity_search;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let results = similarity_search(
        &args[1],
        &args[2],
        args[3].parse().unwrap(),
        args[4].parse().unwrap(),
    );
    for result in results {
        println!("{:?} {:?}", result.line, result.similarity);
    }
}
