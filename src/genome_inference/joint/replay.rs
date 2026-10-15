//! Exact start-domain partitioning; no original-source incidence totals reused.
use super::*;
use crate::genome_inference::{observations::profile, sample::collect_raw_views};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Contribution {
    pub tokens: [u64; 3],
    /// Read length order in Compiled; input orientations averaged, never summed.
    pub totals: Vec<[u64; 2]>,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub contributions: Vec<Contribution>,
    pub admitted_starts: Vec<u64>,
    pub event_runs: u64,
    pub max_crop_bp: u64,
}

pub(super) fn compile_walk(
    panel: &SyngIndex,
    molecule_length: u64,
    lengths: &[u64],
    core_bp: u64,
    mut fetch: impl FnMut(u64, u64) -> io::Result<Vec<u8>>,
) -> io::Result<Profile> {
    require(core_bp > 0, "empty replay core")?;
    let k = panel.syncmer_length_bp() as u64;
    let mut rows: BTreeMap<[u64; 3], Vec<[u64; 2]>> = BTreeMap::new();
    let mut profile = Profile::default();
    for (li, &length) in lengths.iter().enumerate() {
        let starts = if length <= molecule_length {
            molecule_length - length + 1
        } else {
            0
        };
        profile.admitted_starts.push(starts);
        let mut lo = 0;
        while lo < starts {
            let end = lo.saturating_add(core_bp).min(starts);
            // Owned starts [lo,end), read-only right halo. All syncmers in every
            // admitted read are fully contained in this crop, including novel joins.
            let dna = fetch(lo, end - 1 + length)?;
            profile.max_crop_bp = profile.max_crop_bp.max(dna.len() as u64);
            let views = profile::raw_views(panel, &dna)?;
            let mut events = vec![0, end - lo];
            for view in &views {
                for &(_, p) in view {
                    for event in [(p + k).saturating_sub(length), p + 1] {
                        if event > 0 && event < end - lo {
                            events.push(event);
                        }
                    }
                }
            }
            events.sort_unstable();
            events.dedup();
            for range in events.windows(2) {
                let start = range[0];
                let f = profile::restrict(&views[0], start, length, k);
                let r = profile::restrict(&views[1], start, length, k);
                // Native selection/pruning is unchanged. Raw observers are not
                // physical copies. RC swaps the two native passes, so q_rc=q_f.
                for record in collect_raw_views(panel, &f, &r, length)? {
                    for pair in record.windows(2) {
                        let tokens = canonical(&encode_walk(pair)?).try_into().unwrap();
                        let row = rows
                            .entry(tokens)
                            .or_insert_with(|| vec![[0, 0]; lengths.len()]);
                        for q in &mut row[li] {
                            profile::add(q, range[1] - start)?;
                        }
                    }
                }
                profile::add(&mut profile.event_runs, 1)?;
            }
            lo = end;
        }
    }
    profile.contributions = rows
        .into_iter()
        .map(|(tokens, totals)| Contribution { tokens, totals })
        .collect();
    Ok(profile)
}
