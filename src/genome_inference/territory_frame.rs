//! THE TERRITORY FRAME (the global-panel call mode's only new
//! machinery): the committed locality substrate is stated in the
//! territory-LOCAL frame — every member row's path name is a locus.fa
//! record header `path:start-end` (the same convention the cut's own
//! provenance parse and the region-as-interface resolution consume)
//! and every coordinate is slice-relative (0-based within the record).
//! The record header IS the frame statement: the same row in the
//! GLOBAL panel's frame is `(path, start-interval + start, ...)` — a
//! pure coordinate rename, the re-anchor principle (k-independent
//! panel coordinates; conversion-not-materialization).
//!
//! Consumers: `realign_score` (the `--global-frame` engine switch —
//! maps translated on load, product rows translated back on emission)
//! and `sample_projection cut --frame global` (the universe rows).

use crate::sample_mem_bwt::invalid;
use std::collections::HashMap;
use std::io;

/// Split a territory-local record name `path:start-end` into the
/// global panel path name and the record's interval origin. The final
/// `:` is the separator (the committed convention: contig names never
/// carry one; `parse_locus_header` splits the same way).
pub fn split_local_name(name: &str) -> io::Result<(&str, u64)> {
    let (path, interval) = name
        .rsplit_once(':')
        .ok_or_else(|| invalid(&format!("the territory-local name carries no interval: {name:?}")))?;
    let (start, end): (u64, u64) = interval
        .split_once('-')
        .and_then(|(s, e)| Some((s.parse().ok()?, e.parse().ok()?)))
        .ok_or_else(|| invalid(&format!("the territory-local name carries no interval range: {name:?}")))?;
    if start >= end {
        return Err(invalid(&format!(
            "the territory-local record interval is empty: {name:?}"
        )));
    }
    Ok((path, start))
}

/// One translated member row: the global path name and the
/// slice-relative extent shifted into the global frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GlobalRow {
    pub path_name: String,
    pub start: u64,
    pub end: u64,
}

/// Translate one territory-local row `("path:s-e", start, end)` into
/// the global frame. Fail-closed on malformed names.
pub fn translate_row(local_name: &str, start: u64, end: u64) -> io::Result<GlobalRow> {
    let (path, origin) = split_local_name(local_name)?;
    if start >= end {
        return Err(invalid(&format!(
            "the territory-local row is empty: {local_name:?} {start}-{end}"
        )));
    }
    let (global_start, global_end) = (origin.checked_add(start), origin.checked_add(end));
    let (global_start, global_end) = match (global_start, global_end) {
        (Some(start), Some(end)) => (start, end),
        _ => {
            return Err(invalid(&format!(
                "the translated row overflows the coordinate range: {local_name:?} {start}-{end}"
            )))
        }
    };
    Ok(GlobalRow {
        path_name: path.to_string(),
        start: global_start,
        end: global_end,
    })
}

/// The back-translation table: the territory's own record statements
/// (the locality maps' member names), answering "which local record
/// contains this global row" — the product emission's inverse map.
#[derive(Clone, Debug, Default)]
pub struct TerritoryFrame {
    /// global path name -> the record intervals in global coordinates
    /// (the slice-relative local interval shifted by the header
    /// origin), in map order.
    records: HashMap<String, Vec<(u64, u64)>>,
}

impl TerritoryFrame {
    /// Build the table from the territory's local member-row names.
    /// Duplicate identical records (the same name in several windows)
    /// collapse; two DIFFERENT records of one path that overlap are a
    /// named substrate error (the back-translation would be
    /// ambiguous).
    pub fn from_local_names<'a>(names: impl IntoIterator<Item = &'a str>) -> io::Result<Self> {
        let mut frame = TerritoryFrame::default();
        for name in names {
            let (path, origin) = split_local_name(name)?;
            let (path, origin) = (path.to_string(), origin);
            let end = match name.rsplit_once(':').and_then(|(_, interval)| {
                interval
                    .split_once('-')
                    .and_then(|(_, e)| e.parse::<u64>().ok())
            }) {
                Some(end) => end,
                None => {
                    return Err(invalid(&format!(
                        "the territory-local name carries no interval: {name:?}"
                    )))
                }
            };
            let intervals = frame.records.entry(path).or_default();
            if !intervals.iter().any(|&(lo, hi)| lo == origin && hi == end) {
                for &(lo, hi) in intervals.iter() {
                    if lo < end && origin < hi {
                        return Err(invalid(&format!(
                            "the territory carries overlapping records of one path: {name:?}"
                        )));
                    }
                }
                intervals.push((origin, end));
            }
        }
        Ok(frame)
    }

    /// The number of distinct territory records in the table (the
    /// frame's own receipt counter).
    pub fn records_len(&self) -> usize {
        self.records.values().map(|intervals| intervals.len()).sum()
    }

    /// The global panel path name of a local record name.
    pub fn global_name(&self, local_name: &str) -> io::Result<String> {
        Ok(split_local_name(local_name)?.0.to_string())
    }

    /// Back-translate a global row `(path, g_start, g_end)` into the
    /// territory's own frame: the containing record's local name and
    /// the slice-relative extent. Fail-closed when no single record
    /// interval contains the row (a named structural error — the
    /// merged-fold hulls stay inside their member rows by
    /// construction, so a miss is a substrate or machinery divergence
    /// that must be receipted, not guessed).
    pub fn local_row(&self, path_name: &str, start: u64, end: u64) -> io::Result<(String, u64, u64)> {
        let intervals = self
            .records
            .get(path_name)
            .ok_or_else(|| invalid(&format!("the global path {path_name:?} carries no record in the territory frame")))?;
        let mut found = None;
        for &(lo, hi) in intervals {
            if lo <= start && end <= hi {
                if found.is_some() {
                    return Err(invalid(&format!(
                        "the global row {path_name:?}:{start}-{end} is contained in more than one record"
                    )));
                }
                found = Some((format!("{path_name}:{lo}-{hi}"), start - lo, end - lo));
            }
        }
        found.ok_or_else(|| {
            invalid(&format!(
                "the global row {path_name:?}:{start}-{end} lies outside every record of the territory frame"
            ))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_local_name_split_is_the_committed_convention() {
        let (path, origin) = split_local_name("HG002#1#chr1:4632437-4883593").unwrap();
        assert_eq!(path, "HG002#1#chr1");
        assert_eq!(origin, 4632437);
        // contig names with dots carry no colon: a GLOBAL name is
        // refused loudly (it is not a territory-local record name).
        assert!(split_local_name("HG002#1#CM094075.1").is_err());
        assert!(split_local_name("HG002#1#chr1:4632437").is_err());
        assert!(split_local_name("HG002#1#chr1:4883593-4632437").is_err());
    }

    #[test]
    fn the_row_translation_is_the_header_offset() {
        let row = translate_row("HG002#1#chr1:4632437-4883593", 0, 10036).unwrap();
        assert_eq!(row.path_name, "HG002#1#chr1");
        assert_eq!(row.start, 4632437);
        assert_eq!(row.end, 4632437 + 10036);
        assert!(translate_row("HG002#1#chr1:4632437-4883593", 5, 5).is_err());
    }

    #[test]
    fn the_back_translation_round_trips_and_fails_closed() {
        let frame = TerritoryFrame::from_local_names([
            "HG002#1#chr1:4632437-4883593",
            "HG00097#2#CM094075.1:4700806-4951892",
            "HG002#1#chr1:4632437-4883593", // the duplicate collapses
        ])
        .unwrap();
        assert_eq!(
            frame.records.get("HG002#1#chr1").unwrap(),
            &vec![(4632437u64, 4883593u64)]
        );
        let (name, start, end) = frame.local_row("HG002#1#chr1", 4632437, 4642473).unwrap();
        assert_eq!(name, "HG002#1#chr1:4632437-4883593");
        assert_eq!((start, end), (0, 10036));
        // outside every record: refused, never clamped
        assert!(frame.local_row("HG002#1#chr1", 4883000, 4884000).is_err());
        // a path with no record in the territory: refused
        assert!(frame.local_row("CHM13#0#chr1", 0, 100).is_err());
        // empty territory: everything refused
        let empty = TerritoryFrame::default();
        assert!(empty.local_row("HG002#1#chr1", 0, 100).is_err());
    }

    #[test]
    fn overlapping_records_of_one_path_are_a_named_substrate_error() {
        let frame = TerritoryFrame::from_local_names([
            "HG002#1#chr1:100-200",
            "HG002#1#chr1:150-250",
        ]);
        assert!(frame.is_err());
        // disjoint records of one path are fine
        let frame = TerritoryFrame::from_local_names([
            "HG002#1#chr1:100-200",
            "HG002#1#chr1:300-400",
        ])
        .unwrap();
        assert!(frame.local_row("HG002#1#chr1", 350, 360).is_ok());
        // a row spanning the gap between two records is refused (it
        // belongs to no single record)
        assert!(frame.local_row("HG002#1#chr1", 150, 350).is_err());
    }

    #[test]
    fn the_global_name_helper_translates_the_component_lane() {
        let frame = TerritoryFrame::from_local_names(["HG002#1#chr1:4632437-4883593"]).unwrap();
        assert_eq!(
            frame.global_name("HG002#1#chr1:4632437-4883593").unwrap(),
            "HG002#1#chr1"
        );
    }
}
