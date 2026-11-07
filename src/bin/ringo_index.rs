extern crate ringo;

use ringo::db::index::index_file;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let skipped = index_file(&args[1]);
    if skipped > 0 {
        eprintln!("skipped {skipped} line(s) that could not be parsed");
    }
}
