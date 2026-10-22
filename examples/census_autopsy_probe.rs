//! The census-profile autopsy probe (receipt-side, no behavior change
//! anywhere): the L2 locus-17 exception's expressibility receipt.
//!
//! The A/B's one clean global-native win (L2 locus 17: profile B calls
//! the truth pair, the restricted profile A does not) traced to census
//! records whose ONLY anchor is a k-mer the LOCALITY substrate's
//! dictionary never selected, while the GLOBAL panel's dictionary
//! resolves it (membership earned by selection elsewhere in the panel —
//! docs/syng-gmem-bwt/global-panel-call-mode.md §6 risk #1's
//! "conserved-elsewhere k-mer" class). This probe receipts that
//! difference at the exact band, against any panel:
//!
//!   - the canonical anchor walk (the expressibility map's own rule —
//!     `mem_records::canonical_anchor_walk`) over a stated slice, and
//!   - the raw dictionary membership of every k-mer in a stated band
//!     (both orientations, with the node id the dictionary returns).
//!
//! Usage: census_autopsy_probe <panel.syng> <slice.fa> <lo> <hi>
//!
//! The committed receipt run (docs §10.6.17): the L2 window-17 band
//! [1230, 1340) of the HG002#1 and HG002#2 axis slices, against BOTH
//! the regenerated L2 locality panel (the three-step partition's
//! locus-pos64.syng) and the ONE global panel — the array-adjacent
//! k-mers resolve globally and are absent from the locality
//! dictionary, which is why the locality cut's node map (and with it
//! the committed kept-read set — the restricted profile's basis)
//! dropped the single-anchor records that carried the truth allele.

use impg::genome_inference::mem_records;
use impg::syng::{SyncmerParams, SyngIndex};
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "usage: census_autopsy_probe <panel.syng> <slice.fa> <lo> <hi>"
        );
        std::process::exit(2);
    }
    let panel = SyngIndex::load(&args[1], SyncmerParams::default())?;
    let k = panel.syncmer_length_bp();
    let lo: u64 = args[3].parse().map_err(|_| {
        std::io::Error::other(format!("bad lo {}", args[3]))
    })?;
    let hi: u64 = args[4].parse().map_err(|_| {
        std::io::Error::other(format!("bad hi {}", args[4]))
    })?;

    // The slice fasta's FIRST record (the band's own content, any frame).
    let file = std::fs::File::open(&args[2])?;
    let mut lines = BufReader::new(file).lines();
    let _header = lines
        .next()
        .transpose()?
        .ok_or_else(|| std::io::Error::other("empty slice fasta"))?;
    let mut seq: Vec<u8> = Vec::new();
    for line in lines {
        let line = line?;
        if line.starts_with('>') {
            break;
        }
        seq.extend_from_slice(line.trim().as_bytes());
    }
    let seq = seq
        .iter()
        .map(|&b| b.to_ascii_uppercase())
        .collect::<Vec<u8>>();
    eprintln!(
        "[probe] panel {} k={} slice {} bytes band [{lo},{hi})",
        args[1],
        k,
        seq.len()
    );

    // (a) The canonical anchor walk over the slice — the expressibility
    // map's own extraction rule, verbatim.
    let walk = mem_records::canonical_anchor_walk(&panel, &seq)?;
    let band: Vec<(i32, u64)> = walk
        .iter()
        .filter(|&(_, pos)| *pos >= lo && *pos < hi)
        .copied()
        .collect();
    println!("## canonical_anchor_walk anchors in [{lo},{hi}): {}", band.len());
    for (node, pos) in band {
        println!("walk\t{pos}\t{node}");
    }

    // (b) The raw dictionary membership of every k-mer in the band,
    // both orientations (the dictionary's own node id when found).
    let kh = panel.kmer_hash_ptr();
    let find = |kmer: &[u8]| -> Option<i64> {
        let mut buf = kmer.to_vec();
        buf.push(0);
        let mut index: i64 = 0;
        let found = unsafe {
            impg::syng_ffi::kmerHashFind(
                kh,
                buf.as_mut_ptr() as *mut std::os::raw::c_char,
                &mut index,
            )
        };
        if found {
            Some(index)
        } else {
            None
        }
    };
    let revcomp = |seq: &[u8]| -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    };
    println!("## dictionary membership per position in [{lo},{hi}) (fwd-node/rc-node, '-' = absent)");
    for pos in lo..hi {
        let start = pos as usize;
        let Some(kmer) = seq.get(start..start + k) else {
            continue;
        };
        if kmer.contains(&b'N') {
            continue;
        }
        let fwd = find(kmer);
        let rc = find(&revcomp(kmer));
        if fwd.is_some() || rc.is_some() {
            println!(
                "dict\t{pos}\t{}\t{}",
                fwd.map(|i| i.to_string()).unwrap_or_else(|| "-".into()),
                rc.map(|i| i.to_string()).unwrap_or_else(|| "-".into())
            );
        }
    }
    Ok(())
}
