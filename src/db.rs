pub mod exact_search;
pub mod index;
pub mod index_header;
pub mod index_item;
pub mod search;

/// Prefixes an I/O error with the file it is about, since `std::io::Error` does not say.
fn with_path(path: &str) -> impl FnOnce(std::io::Error) -> std::io::Error + '_ {
    move |e| std::io::Error::new(e.kind(), format!("{path}: {e}"))
}
