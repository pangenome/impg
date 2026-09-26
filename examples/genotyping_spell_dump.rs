//! Assessment-side helper for the genotyping scoreboard (EVALUATION ONLY —
//! no model/router/objective code is touched; reads the panel AGC through the
//! repo's own AgcIndex and dumps the FULL spelled sequence of the requested
//! panel sources, so the scoreboard can spell every chain's per-window call
//! with the same sequences the production machinery fetches).
//!
//! Usage: genotyping_spell_dump <source.agc> <out.fa> id:name [id:name ...]
//!
//! The FASTA header is ">[id]:[name]" and the sequence is the source's whole
//! [0,length) span (the strand orientation of a chain segment is applied by
//! the caller, which reverse-complements as needed).

use impg::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use std::io::{BufWriter, Write};

fn fail(msg: String) -> std::io::Error {
    std::io::Error::other(msg)
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        return Err(fail(
            "usage: genotyping_spell_dump <source.agc> <out.fa> id:name [id:name ...]".into(),
        ));
    }
    let agc = args[1].clone();
    let out_path = args[2].clone();
    let mut pairs = Vec::new();
    for spec in &args[3..] {
        let (id, name) = spec
            .split_once(':')
            .ok_or_else(|| fail(format!("bad id:name spec {spec}")))?;
        pairs.push((id.to_string(), name.to_string()));
    }

    let index = UnifiedSequenceIndex::from_files(&[agc])?;

    let out = std::fs::File::create(&out_path)?;
    let mut writer = BufWriter::new(out);
    for (id, name) in &pairs {
        let length = index.get_sequence_length(name)?;
        let dna = index.fetch_sequence(name, 0, length as i32)?;
        if dna.len() != length {
            return Err(fail(format!(
                "fetch length mismatch for {name}: {length} vs {}",
                dna.len()
            )));
        }
        writeln!(writer, ">{}:{}", id, name)?;
        for chunk in dna.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}
