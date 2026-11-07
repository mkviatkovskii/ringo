#[expect(
    clippy::unwrap_used,
    reason = "integration helper requires valid fixtures and successful file operations"
)]
pub fn assert_search<'a>(name: &str, records: &str, queries: impl IntoIterator<Item = &'a str>) {
    let path = std::env::temp_dir().join(format!("ringo-{name}-{}.smi", std::process::id()));
    std::fs::write(&path, records).unwrap();
    let filename = path.to_str().unwrap();
    let skipped =
        ringo::db::index::index_file(filename, 2, 2048, ringo::db::index::DEFAULT_BATCH_LINES)
            .unwrap();
    assert_eq!(skipped, 0);
    for query in queries {
        let results = ringo::db::search::similarity_search(
            filename,
            query,
            1.0,
            1.0,
            10,
            ringo::db::search::SearchMode::Best,
        )
        .unwrap();
        assert_eq!(results.len(), 1, "{query}");
        assert_eq!(results[0].similarity, 1.0, "{query}");
    }
    std::fs::remove_file(filename.to_owned() + ".fp").unwrap();
    std::fs::remove_file(path).unwrap();
}
