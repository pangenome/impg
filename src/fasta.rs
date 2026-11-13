//! One FASTA parser for plain and compressed (niffler) files.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

fn open_niffler_reader(path: &Path) -> io::Result<Box<dyn Read + Send>> {
    let file = File::open(path)?;
    let (reader, _format) = niffler::send::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!(
            "Failed to open reader for '{}': {e}",
            path.display()
        ))
    })?;
    Ok(reader)
}

/// Open a plain or compressed file for line reading.
pub fn open_reader(path: &Path) -> io::Result<BufReader<Box<dyn Read + Send>>> {
    Ok(BufReader::new(open_niffler_reader(path)?))
}

/// Open a plain or compressed file with an explicit `BufReader` buffer capacity.
pub fn open_reader_with_capacity(
    path: &Path,
    capacity: usize,
) -> io::Result<BufReader<Box<dyn Read + Send>>> {
    Ok(BufReader::with_capacity(
        capacity,
        open_niffler_reader(path)?,
    ))
}

/// Write a single FASTA record wrapped to 80 bases per line.
pub fn write_wrapped_fasta<W: Write + ?Sized>(
    out: &mut W,
    name: &str,
    seq: &[u8],
) -> io::Result<()> {
    writeln!(out, ">{name}")?;
    for chunk in seq.chunks(80) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// Pass every record to `f`: the name is the first token after '>', the
/// sequence the trimmed lines that follow. Lines before the first header
/// are dropped. Empty names and empty sequences are passed on.
pub fn for_each_record(
    path: &Path,
    mut f: impl FnMut(String, Vec<u8>) -> io::Result<()>,
) -> io::Result<()> {
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    for line in open_reader(path)?.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            let next = header.split_whitespace().next().unwrap_or("").to_string();
            match name.replace(next) {
                Some(done) => f(done, std::mem::take(&mut seq))?,
                None => seq.clear(),
            }
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(done) = name {
        f(done, seq)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn records_keep_first_token_and_drop_preamble() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.fa");
        std::fs::write(&path, "junk\n>a desc\nAC\n GT \n>\n>b\n\n>c\nT\n").unwrap();
        let mut got = Vec::new();
        super::for_each_record(&path, |n, s| {
            got.push((n, String::from_utf8(s).unwrap()));
            Ok(())
        })
        .unwrap();
        let want = [("a", "ACGT"), ("", ""), ("b", ""), ("c", "T")];
        assert_eq!(got, want.map(|(n, s)| (n.to_string(), s.to_string())));
    }

    #[test]
    fn write_wrapped_fasta_wraps_at_80_and_handles_empty() {
        let mut out = Vec::new();
        super::write_wrapped_fasta(&mut out, "empty", b"").unwrap();
        super::write_wrapped_fasta(&mut out, "long", &[b'A'; 85]).unwrap();
        let expected = format!(">empty\n>long\n{}\nAAAAA\n", "A".repeat(80));
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }
}
