//! Bounded JSON framing. Preflight validates the whole syntax/EOF without an index.
use super::*;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
pub struct Json {
    reader: BufReader<File>,
    pub position: u64,
    cap: usize,
}
impl Json {
    pub fn open(path: &Path, cap: usize) -> io::Result<Self> {
        Ok(Self {
            reader: BufReader::with_capacity(65536, File::open(path)?),
            position: 0,
            cap,
        })
    }
    pub fn seek(&mut self, position: u64) -> io::Result<()> {
        self.reader.seek(SeekFrom::Start(position))?;
        self.position = position;
        Ok(())
    }
    fn peek(&mut self) -> io::Result<Option<u8>> {
        Ok(self.reader.fill_buf()?.first().copied())
    }
    fn byte(&mut self) -> io::Result<u8> {
        let b = self.peek()?.ok_or_else(|| invalid("premature JSON EOF"))?;
        self.reader.consume(1);
        self.position += 1;
        Ok(b)
    }
    pub fn ws(&mut self) -> io::Result<()> {
        while self.peek()?.is_some_and(|b| b" \r\n\t".contains(&b)) {
            self.byte()?;
        }
        Ok(())
    }
    pub fn expect(&mut self, b: u8) -> io::Result<()> {
        self.ws()?;
        ensure(self.byte()? == b, "JSON delimiter mismatch")
    }
    pub fn end(&mut self) -> io::Result<()> {
        self.ws()?;
        ensure(self.peek()?.is_none(), "trailing JSON bytes")
    }
    pub fn line_end(&mut self) -> io::Result<()> {
        while self.peek()?.is_some_and(|b| b" \t\r".contains(&b)) {
            self.byte()?;
        }
        match self.peek()? {
            Some(b'\n') => {
                self.byte()?;
                Ok(())
            }
            None => Ok(()),
            _ => Err(invalid("multiple values in ledger line")),
        }
    }
    fn push(&self, v: &mut Vec<u8>, b: u8) -> io::Result<()> {
        ensure(v.len() < self.cap, "input-record-byte-budget-exhausted")?;
        v.push(b);
        Ok(())
    }
    fn string(&mut self) -> io::Result<Vec<u8>> {
        self.expect(b'"')?;
        let mut v = vec![b'"'];
        let mut escaped = false;
        loop {
            let b = self.byte()?;
            self.push(&mut v, b)?;
            if !escaped && b == b'"' {
                break;
            }
            if !escaped && b == b'\\' {
                escaped = true;
            } else {
                escaped = false;
            }
        }
        // This also rejects invalid escapes, controls, and invalid UTF-8.
        let _: String = serde_json::from_slice(&v).map_err(io::Error::other)?;
        Ok(v)
    }
    pub fn key(&mut self) -> io::Result<String> {
        let v = self.string()?;
        self.expect(b':')?;
        serde_json::from_slice(&v).map_err(io::Error::other)
    }
    /// No retained value when capture=false, even for multi-GB container fields.
    fn value_into(&mut self, depth: usize, mut out: Option<&mut Vec<u8>>) -> io::Result<()> {
        ensure(depth <= 128, "JSON nesting budget exhausted")?;
        self.ws()?;
        let b = self.peek()?.ok_or_else(|| invalid("missing JSON value"))?;
        if b == b'{' || b == b'[' {
            self.byte()?;
            if let Some(v) = out.as_deref_mut() {
                self.push(v, b)?;
            }
            let end = if b == b'{' { b'}' } else { b']' };
            self.ws()?;
            if self.peek()? != Some(end) {
                loop {
                    if b == b'{' {
                        let key = self.string()?;
                        self.expect(b':')?;
                        if let Some(v) = out.as_deref_mut() {
                            for c in key {
                                self.push(v, c)?;
                            }
                            self.push(v, b':')?;
                        }
                    }
                    self.value_into(depth + 1, out.as_deref_mut())?;
                    self.ws()?;
                    if self.peek()? == Some(end) {
                        break;
                    }
                    self.expect(b',')?;
                    if let Some(v) = out.as_deref_mut() {
                        self.push(v, b',')?;
                    }
                }
            }
            self.expect(end)?;
            if let Some(v) = out {
                self.push(v, end)?;
            }
        } else {
            let v = if b == b'"' {
                self.string()?
            } else {
                let mut v = Vec::new();
                while self.peek()?.is_some_and(|c| !b" \r\n\t,]}".contains(&c)) {
                    let c = self.byte()?;
                    self.push(&mut v, c)?;
                }
                let _: Value = serde_json::from_slice(&v).map_err(io::Error::other)?;
                v
            };
            if let Some(out) = out {
                for b in v {
                    self.push(out, b)?;
                }
            }
        }
        Ok(())
    }
    pub fn skip(&mut self) -> io::Result<()> {
        self.value_into(0, None)
    }
    pub fn raw(&mut self) -> io::Result<Vec<u8>> {
        let mut v = Vec::new();
        self.value_into(0, Some(&mut v))?;
        Ok(v)
    }
    pub fn value<T: serde::de::DeserializeOwned>(&mut self) -> io::Result<T> {
        serde_json::from_slice(&self.raw()?).map_err(io::Error::other)
    }
    /// Called before the first element, then after each value. Rejects trailing commas.
    pub fn next(&mut self, first: &mut bool, end: u8) -> io::Result<bool> {
        self.ws()?;
        if self.peek()? == Some(end) {
            self.byte()?;
            return Ok(false);
        }
        if !*first {
            self.expect(b',')?;
        }
        *first = false;
        Ok(true)
    }
}
pub fn hash_file(path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut hash = 0xcbf29ce484222325u64;
    let mut bytes = 0u64;
    let mut buf = [0; 65536];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        bytes += n as u64;
        for b in &buf[..n] {
            hash = (hash ^ *b as u64).wrapping_mul(0x100000001b3);
        }
    }
    Ok(format!("fnv1a64:{bytes}:{hash:016x}"))
}
#[cfg(test)]
mod tests {
    use super::*;
    fn parse(bytes: &[u8], cap: usize) -> io::Result<()> {
        let mut f = tempfile::NamedTempFile::new()?;
        f.write_all(bytes)?;
        let mut j = Json::open(f.path(), cap)?;
        j.skip()?;
        j.end()
    }
    #[test]
    fn bounded_stream_syntax_eof_and_large_skipped_containers() {
        for bad in [
            b"{\"a\":[1,]}".as_slice(),
            b"{}x",
            b"{",
            b"[01]",
            b"\"\\x\"",
            b"[true false]",
        ] {
            assert!(parse(bad, 64).is_err(), "{bad:?}");
        }
        assert!(parse(b"{\"a\":[1,true,null,\"x\\\"y\"]}", 64).is_ok());
        let large = format!("[{}0]", "0,".repeat(100000));
        assert!(parse(large.as_bytes(), 16).is_ok());
        assert!(parse(b"\"abcdefghijklmnop\"", 16).is_err());
    }
}
