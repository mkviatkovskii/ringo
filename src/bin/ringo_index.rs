extern crate ringo;

use ringo::db::index::index_file;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    index_file(&args[1]);
}
