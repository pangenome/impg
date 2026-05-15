//! `impg subgraph` — extract a local syng index + GFA for a query region.
//!
//! Pipeline: full-syng `query_region_ext` → pad/merge target hits →
//! pull DNA slices from FASTA/AGC (named `contig:start-end`) → build a
//! local `SyngIndex` over the slices → dump to GFA via `syng2gfa`.

use std::fs::File;
use std::io::{self, BufWriter, Write};

use log::{info, warn};

use crate::commands::partition::parse_target_range;
use crate::commands::syng2gfa::{self, GfaVersion};
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::syng::{SyncmerParams, SyngIndex};

pub struct SubgraphOpts {
    pub index_prefix: String,
    pub sequence_files: Vec<String>,
    pub region: String,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub pad: u64,
    pub merge_distance: Option<u64>,
    pub output_prefix: String,
    pub gfa_version: GfaVersion,
    pub position_sample_rate: u32,
    pub blunt: bool,
}

struct Range {
    genome: String,
    start: u64,
    end: u64,
}

pub fn run(opts: SubgraphOpts) -> io::Result<()> {
    // ---- load full index + sequence index ----
    info!("[subgraph] loading full syng index from '{}'", opts.index_prefix);
    let full = SyngIndex::load(&opts.index_prefix, SyncmerParams::default())?;
    let syncmer_len = full.syncmer_length_bp() as u64;

    if opts.sequence_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "subgraph requires --sequence-files / --sequence-list to extract region DNA",
        ));
    }
    info!(
        "[subgraph] opening {} sequence file(s)",
        opts.sequence_files.len()
    );
    let seq_idx = UnifiedSequenceIndex::from_files(&opts.sequence_files)?;

    // ---- parse query region ----
    let (q_genome, (q_start_i, q_end_i), q_name) = parse_target_range(&opts.region)?;
    let q_start = q_start_i.max(0) as u64;
    let q_end = q_end_i.max(0) as u64;
    info!("[subgraph] query region: {}", q_name);

    // ---- query full index for homologous target ranges ----
    let hits = full.query_region_ext(
        &q_genome,
        q_start,
        q_end,
        opts.syng_padding,
        opts.syng_extension,
    )?;
    info!(
        "[subgraph] full-index query returned {} raw target interval(s)",
        hits.len()
    );

    let mut ranges: Vec<Range> = hits
        .into_iter()
        .map(|h| Range {
            genome: h.genome,
            start: h.start,
            end: h.end,
        })
        .collect();
    ranges.push(Range {
        genome: q_genome,
        start: q_start,
        end: q_end,
    });

    // ---- pad (clamped to contig length) + sort + sweep-merge ----
    let merge_dist = opts.merge_distance.unwrap_or(opts.pad);
    let contig_len = |name: &str| -> u64 {
        full.name_map
            .name_to_path
            .get(name)
            .and_then(|idx| full.name_map.path_to_length.get(*idx as usize))
            .copied()
            .unwrap_or(u64::MAX)
    };
    for r in &mut ranges {
        let len = contig_len(&r.genome);
        r.start = r.start.saturating_sub(opts.pad);
        r.end = r.end.saturating_add(opts.pad).min(len);
    }
    ranges.sort_by(|a, b| (a.genome.as_str(), a.start).cmp(&(b.genome.as_str(), b.start)));

    let mut merged: Vec<Range> = Vec::with_capacity(ranges.len());
    for r in ranges {
        match merged.last_mut() {
            Some(cur) if cur.genome == r.genome && r.start <= cur.end.saturating_add(merge_dist) => {
                cur.end = cur.end.max(r.end);
            }
            _ => merged.push(r),
        }
    }
    info!(
        "[subgraph] {} merged region(s) after pad={} merge_distance={}",
        merged.len(),
        opts.pad,
        merge_dist
    );

    // ---- write <o>.regions.bed ----
    let bed_path = format!("{}.regions.bed", opts.output_prefix);
    {
        let mut w = BufWriter::new(File::create(&bed_path)?);
        for r in &merged {
            writeln!(w, "{}\t{}\t{}", r.genome, r.start, r.end)?;
        }
        w.flush()?;
    }
    info!("[subgraph] wrote {}", bed_path);

    // ---- extract DNA slices, write <o>.local.fa ----
    let fa_path = format!("{}.local.fa", opts.output_prefix);
    let mut extracted: Vec<(String, Vec<u8>)> = Vec::with_capacity(merged.len());
    {
        let mut w = BufWriter::new(File::create(&fa_path)?);
        let mut total_bp: u64 = 0;
        for r in &merged {
            if r.end > i32::MAX as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "region end {} on '{}' exceeds i32 range (fetch_sequence limit)",
                        r.end, r.genome
                    ),
                ));
            }
            let seq = seq_idx.fetch_sequence(&r.genome, r.start as i32, r.end as i32)?;
            if (seq.len() as u64) < syncmer_len {
                warn!(
                    "[subgraph] skipping {}:{}-{} ({} bp < syncmer length {})",
                    r.genome,
                    r.start,
                    r.end,
                    seq.len(),
                    syncmer_len
                );
                continue;
            }
            let name = format!("{}:{}-{}", r.genome, r.start, r.end);
            writeln!(w, ">{}", name)?;
            for chunk in seq.chunks(60) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
            total_bp += seq.len() as u64;
            extracted.push((name, seq));
        }
        w.flush()?;
        info!(
            "[subgraph] wrote {} ({} sequences, {} bp)",
            fa_path,
            extracted.len(),
            total_bp
        );
    }
    if extracted.is_empty() {
        return Err(io::Error::other(
            "no extractable regions (all shorter than syncmer length)",
        ));
    }

    // ---- build local syng index over extracted slices ----
    info!(
        "[subgraph] building local syng index ({} paths, sample_rate={})",
        extracted.len(),
        opts.position_sample_rate
    );
    let mut local = SyngIndex::new(full.params);
    local.enable_online_sampled_positions(opts.position_sample_rate)?;
    for (name, seq) in extracted {
        local.add_sequence(name, seq);
    }
    local.save(&opts.output_prefix)?;
    info!(
        "[subgraph] saved local syng index: {}.{{1gbwt,1khash,syng.*}} ({} syncmer nodes)",
        opts.output_prefix,
        local.num_syncmer_nodes()
    );

    // ---- dump local index to GFA (gap-fill from local.fa) ----
    let gfa_path = format!("{}.gfa", opts.output_prefix);
    let local_seq = UnifiedSequenceIndex::from_files(&[fa_path.clone()])?;
    let mut buf: Vec<u8> = Vec::new();
    let (s, l, p, skipped, gaps, gap_bp) =
        syng2gfa::write_gfa(&local, &mut buf, opts.gfa_version, Some(&local_seq))?;
    std::fs::write(&gfa_path, &buf)?;
    info!(
        "[subgraph] wrote {}: {} S ({} syncmer + {} gap, {} gap bp), {} L, {} P/W ({} skipped)",
        gfa_path,
        s,
        s - gaps,
        gaps,
        gap_bp,
        l,
        p,
        skipped
    );

    if opts.blunt {
        if opts.gfa_version != GfaVersion::V1_0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--blunt requires --gfa-version 1.0 (bluntg parses P lines only)",
            ));
        }
        let g = bluntg::parse_gfa(&buf);
        let blunt = bluntg::bluntify_auto(g);
        let blunt_path = format!("{}.blunt.gfa", opts.output_prefix);
        let mut w = BufWriter::new(File::create(&blunt_path)?);
        bluntg::write_gfa(&mut w, &blunt)?;
        w.flush()?;
        info!(
            "[subgraph] wrote {} ({} S, {} L, {} P, all 0M)",
            blunt_path,
            blunt.segments.len(),
            blunt.links.len(),
            blunt.paths.len()
        );
    }

    Ok(())
}
