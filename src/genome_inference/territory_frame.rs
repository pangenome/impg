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
//!
//! THE RECORD-ORIGIN CONVENTION (the `--record-origin <d>`
//! parameter): the committed substrate's locus.fa content need not
//! sit at its header's own interval — the L1-L4 substrates' content
//! sits ONE BASE LATER (d=1, the uniform delta=1 AGC-identity
//! receipts), the pilot's at the header interval exactly (d=0).
//! The global frame, which deletes locus.fa by design, must be TOLD
//! the convention: it shifts every record's content origin to
//! `header_start + d` (default 0 = the header's own interval). The
//! committed local frame derives it per record against the AGC
//! (fail-closed on neither convention matching) and needs no
//! parameter. The back-translation cancels d exactly (the same
//! shift lands in the forward table and the inverse), so the product
//! emission is origin-independent.

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
/// the global frame, at the substrate's stated record origin: the
/// record's content sits at panel interval `[s + d, e + d)` (d = the
/// `--record-origin` convention: 0 = the header's own interval, 1 =
/// the committed L1-L4 substrates' one-base-later content). Fail-closed
/// on malformed names.
pub fn translate_row(local_name: &str, start: u64, end: u64, record_origin: u64) -> io::Result<GlobalRow> {
    let (path, origin) = split_local_name(local_name)?;
    if start >= end {
        return Err(invalid(&format!(
            "the territory-local row is empty: {local_name:?} {start}-{end}"
        )));
    }
    let origin = origin
        .checked_add(record_origin)
        .ok_or_else(|| invalid(&format!(
            "the record origin overflows the coordinate range: {local_name:?} + {record_origin}"
        )))?;
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
    /// The substrate's record-origin convention (the `--record-origin`
    /// d): every record's content sits at `header_start + d` — the
    /// forward translation and this table's containment share it, so
    /// the back-translation cancels it exactly.
    record_origin: u64,
    /// global path name -> the record statements in map order, as the
    /// HEADER states them: `(start, end)` of `path:start-end` — the
    /// record NAME for the back-translation; the containment interval
    /// is this shifted by `record_origin`.
    records: HashMap<String, Vec<(u64, u64)>>,
}

impl TerritoryFrame {
    /// Build the table from the territory's local member-row names,
    /// at the substrate's stated record-origin convention `d`
    /// (0 = the header's own interval; 1 = the committed L1-L4
    /// substrates' one-base-later content). Duplicate identical records
    /// (the same name in several windows) collapse; two DIFFERENT
    /// records of one path that overlap are a named substrate error
    /// (the back-translation would be ambiguous). The overlap check
    /// runs on the CONTENT intervals (`start + d`), the frame the
    /// machinery actually occupies.
    pub fn from_local_names<'a>(
        names: impl IntoIterator<Item = &'a str>,
        record_origin: u64,
    ) -> io::Result<Self> {
        let mut frame = TerritoryFrame {
            record_origin,
            ..TerritoryFrame::default()
        };
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
                // The overlap test is origin-independent (the uniform d
                // shifts both records alike, so it cancels).
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
    /// territory's own frame: the containing record's local name (the
    /// HEADER's own interval — the committed output shape, whatever
    /// the record-origin convention) and the slice-relative extent,
    /// measured from the record's CONTENT origin `start + d` (the same
    /// shift the forward translation applied, so the round trip is
    /// the identity and the product emission is origin-independent).
    /// Fail-closed when no single record interval contains the row
    /// (a named structural error — the merged-fold hulls stay inside
    /// their member rows by construction, so a miss is a substrate or
    /// machinery divergence that must be receipted, not guessed).
    pub fn local_row(&self, path_name: &str, start: u64, end: u64) -> io::Result<(String, u64, u64)> {
        let intervals = self
            .records
            .get(path_name)
            .ok_or_else(|| invalid(&format!("the global path {path_name:?} carries no record in the territory frame")))?;
        let mut found = None;
        for &(lo, hi) in intervals {
            let content_lo = lo + self.record_origin;
            let content_hi = hi + self.record_origin;
            if content_lo <= start && end <= content_hi {
                if found.is_some() {
                    return Err(invalid(&format!(
                        "the global row {path_name:?}:{start}-{end} is contained in more than one record"
                    )));
                }
                found = Some((format!("{path_name}:{lo}-{hi}"), start - content_lo, end - content_lo));
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
        let row = translate_row("HG002#1#chr1:4632437-4883593", 0, 10036, 0).unwrap();
        assert_eq!(row.path_name, "HG002#1#chr1");
        assert_eq!(row.start, 4632437);
        assert_eq!(row.end, 4632437 + 10036);
        assert!(translate_row("HG002#1#chr1:4632437-4883593", 5, 5, 0).is_err());
    }

    #[test]
    fn the_record_origin_shifts_the_content_to_the_substrate_convention() {
        // d=1: the committed L1-L4 substrates' content sits ONE BASE
        // LATER than the header interval (the uniform delta=1
        // AGC-identity receipts).
        let row = translate_row("HG002#1#chr1:4632437-4883593", 0, 10036, 1).unwrap();
        assert_eq!(row.start, 4632438);
        assert_eq!(row.end, 4632438 + 10036);
        // overflow of the shifted origin: fail-closed, never wrapped
        assert!(translate_row(
            "HG002#1#chr1:18446744073709551614-18446744073709551615",
            0,
            1,
            3,
        )
        .is_err());
    }

    #[test]
    fn the_back_translation_round_trips_and_fails_closed() {
        let frame = TerritoryFrame::from_local_names(
            [
                "HG002#1#chr1:4632437-4883593",
                "HG00097#2#CM094075.1:4700806-4951892",
                "HG002#1#chr1:4632437-4883593", // the duplicate collapses
            ],
            0,
        )
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
    fn the_record_origin_cancels_in_the_back_translation() {
        // THE ORIGIN-INDEPENDENCE PROPERTY: a row translated forward at
        // d back-translates to EXACTLY its own local extent (the same
        // shift lands in both directions), and the record NAME keeps
        // the header's own interval — the committed output shape.
        let name = "HG002#1#chr1:4632437-4883593";
        for d in [0u64, 1] {
            let frame = TerritoryFrame::from_local_names([name], d).unwrap();
            let translated = translate_row(name, 100, 200, d).unwrap();
            let (back, start, end) = frame
                .local_row(&translated.path_name, translated.start, translated.end)
                .unwrap();
            assert_eq!(back, name);
            assert_eq!((start, end), (100, 200));
        }
        // The d=1 content interval: containment is measured from
        // s+1, so the header's own first base [s, s+1) is OUTSIDE the
        // record (fail-closed), while [s+1, s+2) is inside.
        let frame = TerritoryFrame::from_local_names([name], 1).unwrap();
        assert!(frame.local_row("HG002#1#chr1", 4632437, 4632438).is_err());
        let (_, start, end) = frame.local_row("HG002#1#chr1", 4632438, 4632439).unwrap();
        assert_eq!((start, end), (0, 1));
    }

    #[test]
    fn overlapping_records_of_one_path_are_a_named_substrate_error() {
        // the overlap test is origin-independent: the same substrate
        // error at both conventions
        for d in [0u64, 1] {
            let frame = TerritoryFrame::from_local_names(
                ["HG002#1#chr1:100-200", "HG002#1#chr1:150-250"],
                d,
            );
            assert!(frame.is_err());
            // disjoint records of one path are fine
            let frame = TerritoryFrame::from_local_names(
                ["HG002#1#chr1:100-200", "HG002#1#chr1:300-400"],
                d,
            )
            .unwrap();
            assert!(frame.local_row("HG002#1#chr1", 350, 360).is_ok());
            // a row spanning the gap between two records is refused (it
            // belongs to no single record)
            assert!(frame.local_row("HG002#1#chr1", 150, 350).is_err());
        }
    }

    #[test]
    fn the_global_name_helper_translates_the_component_lane() {
        let frame = TerritoryFrame::from_local_names(["HG002#1#chr1:4632437-4883593"], 0).unwrap();
        assert_eq!(
            frame.global_name("HG002#1#chr1:4632437-4883593").unwrap(),
            "HG002#1#chr1"
        );
    }
}
