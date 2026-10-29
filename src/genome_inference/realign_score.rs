//! THE SCORING LAYER (slice B of the local-realignment evidence model;
//! owner direction 2026-11-06, slice A's anchor-projection layer stands
//! as the placement machinery of record). Assessment-side only.
//!
//! THE MODEL, derived and stated (no thresholds, no tuning constants):
//!
//! 1. THE SCORING RATE comes from the reads' own base qualities: the
//!    sample's FASTQ carries one uniform quality byte ('I' = Phred 40,
//!    asserted by streaming the whole file), so the per-base error
//!    probability is eps = 10^(-40/10) = 1e-4 and the substitution
//!    model gives P(match) = 1 - eps, P(mismatch) = eps/3 (uniform over
//!    the three non-read bases). The per-base log weights are
//!    A = ln(1 - eps) and B = ln(eps/3). The measured error profile is
//!    stated beside them: slice A measured zero unmatched placement bp
//!    over 1,063,522 verified occurrences (every base of every chained
//!    placement span lies inside an exactly-matched anchor window), so
//!    the quality byte is the sample's own statement of its error rate.
//!
//! 2. PER (RECORD, CANDIDATE ROW) SKIPPED-BASE VOTES, computed
//!    in-process from slice A's committed machinery (the anchor
//!    correspondence, the serving-anchor projection, the mirror
//!    equations — the same functions the slice A receipts and checker
//!    validated; the tens-of-GB full-fidelity per-path context is
//!    computed here, never materialized). A candidate row's
//!    log-likelihood for one read instance is
//!
//!    ```text
//!    LL(read | row) = -ln(2 * (len_row - 149))
//!                     + logsumexp over anchored placements of
//!                       [ backbone term + skipped-base evidence ]
//!    ```
//!
//!    * THE BACKBONE TERM is the MEM-anchored shared material: the
//!      merged bp of the PINNED anchor windows (the anchors the row's
//!      own walk shares with the required sign, placed by the
//!      unique-monotone correspondence), each a match by node identity
//!      and orientation sign. It is identical across near-identical
//!      paths (rows sharing the anchor chain share the backbone) and
//!      cancels in their comparisons; it is computed once per
//!      (record, placement), never per base.
//!    * THE SKIPPED-BASE EVIDENCE is the read's bases outside its
//!      merged anchor windows (the flanks; interior hull gaps measured
//!      zero on this sample, emitted honestly per mirror state when
//!      present): each skipped base is projected onto the row through
//!      its serving anchor (the bracketing pinned anchor, slice A's
//!      committed rule) and votes A on a match, B on a mismatch
//!      against the row's spelled sequence (the read base compared as
//!      its complement at physically-reverse placements). Bases
//!      projecting outside the row's extent carry no vote (the row
//!      spells nothing there) — and any placement whose voted bases do
//!      not all project inside the extent is not a placement of the
//!      model (the read must lie fully inside the row, the uniform
//!      placement prior's own support).
//!    * THE PLACEMENT PRIOR is uniform over the row's full-read
//!      placements, both strands: 2 * (len_row - 149) of them — the
//!      derived count, stated; the anchored placements are the
//!      dominant terms of the sum (the exactness receipt measures
//!      the dominance directly).
//!    * Bases inside anchor windows the row does NOT pin (absent or
//!      ambiguous anchors) abstain — the named remainder from slice A
//!      ("the within-anchor variant positions"), counted and reported
//!      per locus, not silently scored.
//!
//! 3. NO SHARE-SPLITTING: each candidate row is an independent
//!    hypothesis; a record votes ONCE per row (its instances'
//!    likelihoods under that row are computed independently of every
//!    other row and of how many occurrences discovered the placement —
//!    multiple occurrences of one record pin the SAME placement on a
//!    row and count once). The record's instances are its evidence:
//!    every distinct read variant (mirror state, hull start, sequence)
//!    votes with its own instance count.
//!
//! 4. THE LOCALITY RULE (the ruling for occurrences outside built
//!    axis-partition rows, measured before implementation): a locus's
//!    evidence set is the records with at least one occurrence
//!    TOUCHING the locus's window (the routing's own territory-touch
//!    rule — the same convention that defines the current instrument's
//!    per-locus record sets), and the occurrences that vote at the
//!    locus are exactly those touching occurrences, whether or not
//!    the occurrence's own placement is contained in a built
//!    axis-partition row. The measured outside class is row-edge
//!    overhangs (the read's anchor hull crosses its own path's
//!    partition-membership row boundary) plus territory-extension
//!    material — both window-seam artifacts of the alignment-induced
//!    membership, and the likelihood a read supports a row is a
//!    property of (read, row sequence), not of a row boundary; the
//!    placement-validity rule already guarantees such a read testifies
//!    only where the row actually spells its bases.
//!
//! 5. CANDIDATES fold by identical sequence AND walk (coalesced copies
//!    are ONE hypothesis — the identical-through-graph pair
//!    BTE#3/#4 block28_contig1 folds by construction; the fold is
//!    verified against the committed GFAs by the checker). Classes are
//!    all unordered fold pairs i <= j; class LL = sum over the locus's
//!    evidence units (record, variant) of count * ln(0.5 e^{LL_i} +
//!    0.5 e^{LL_j}) — the equal 1/2 mixture is the balanced 15x/15x
//!    diploid prior. QUAL is the existing cluster machinery verbatim
//!    (spectrum_knee / cluster_form_qual / qual_p / qual_from_p and
//!    the material-distance spectrum, mirrored pure functions from
//!    panel_route_spine/cosine_probe.rs) over the likelihood ratios.
//!
//! 5b. THE ALIGNMENT-CONNECTED FOLD RELAXATION (the segdup follow-up's
//!    stage-1 rule; env IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS,
//!    beside the exact rule — env-unset is the exact code path
//!    verbatim). The exact fold rule coalesces rows only on identical
//!    spelled sequence AND identical contained walk — exact node+edge
//!    multiset equality at identical offsets — which is too strict at
//!    near-identical segdup copies: the truth's second haplotype
//!    exists as graph paths but the variant pockets (the
//!    node-multiset deltas between the near-identical copies) break
//!    the exact multiset, so the copies never coalesce and the exact
//!    correspondence refuses. Under the relaxation rows fold by their
//!    CONNECTION in the locality's alignment structure — the pggb
//!    bubble: the connected alignment component containing the axis
//!    correspondence (folds are alignment-adjacent iff their contained
//!    node multisets share a node) — SUBJECT TO the safety property
//!    that preserves the fold's guarantee, A FOLD STILL MEANS NO
//!    CANDIDATE CAN DISTINGUISH THE MEMBERS: the members must be
//!    indistinguishable AT READ LEVEL — identical on the
//!    observed-material projection, i.e. every evidence unit's
//!    per-fold log-likelihood column is BIT-IDENTICAL across the
//!    members (the operative exactly-verified predicate, computed
//!    from the scoring matrix itself), whose data-level form is the
//!    derived condition: every node and edge in the symmetric
//!    difference of the members' contained multisets — the variant
//!    pockets — carries ZERO observed census mass at this locus, OR
//!    THE FOLD DOES NOT HAPPEN (the merge is refused and the exact
//!    folds stand; an SNV pocket with even one observed read refuses
//!    the merge — the property that prevents folding
//!    genuinely-distinguishable material). The same alignment
//!    connection carries the truth-second correspondence: the fold
//!    carrying the second haplotype's rows with the MAXIMUM
//!    shared-node mass with the axis fold (argmax, not a threshold;
//!    ties and disconnected candidates refuse — no fabrication),
//!    replacing the exact rule's single-window/exactly-one-fold
//!    attribution refusal. The merged fold carries the representative's
//!    sequence/walk/skeleton (any member is valid under the asserted
//!    invariant) and the additive multiset merge; the report names
//!    every merge group with its pocket masses and refusal counters,
//!    and the truth-QV artifact carries the correspondence receipt
//!    (checker inputs). Merging LL-identical columns preserves every
//!    class LL value, so winner and truth likelihoods are unchanged;
//!    fold indices, class counts, called strain lists and QUAL spectra
//!    can change where merges land — every change attributable to a
//!    named merge group.
//!
//! 6. THE EXACTNESS PROOF: for EVERY scored (unit, fold, placement) in
//!    the pilot domain the factorized score is asserted equal to the
//!    direct per-base comparison (the backbone verified base-by-base
//!    against the row's sequence, not trusted from node identity) —
//!    integer (match, mismatch) equality, loudly failing. A measured
//!    evenly-spaced sample of (read, candidate row) pairs — plus the
//!    named strata: the truth folds, the identical-through-graph fold,
//!    and pairs with mismatch votes (the variant-pocket class) — is
//!    emitted with vote-level detail for the independent checker,
//!    which re-derives the scores from the committed GFAs, runs the
//!    bounded-offset dominance scan and the bounded edit-distance
//!    alignment on the sample.
//!
//! 7. THE GENERATED-HERE/GENERATED-ELSEWHERE MARGINALIZATION (slice C,
//!    the owner's rho-squared principle: conserved material cancels;
//!    locally-variable material discriminates). A locus's record set
//!    is defined by the routing's territory touch, so conserved
//!    pockets bring foreign-row-generated and cross-window reads into
//!    the locale's evidence; a candidate that spells such material
//!    would otherwise collect it at full match score while a candidate
//!    that does not pays the all-mismatch floor — even though the
//!    read's own genome generated it somewhere the candidate models
//!    only as "outside this locale". Each read's likelihood under a
//!    candidate therefore MARGINALIZES the two generation branches:
//!
//!    ```text
//!    LL(read | candidate) = logsumexp( local-spell branch ,
//!                                      E(read) )
//!    ```
//!
//!    where the local-spell branch is the per-(record, row) likelihood
//!    of rule 2 (the uniform placement prior + the anchored
//!    placements' backbone + skipped-base votes) and E(read) is the
//!    CANDIDATE-INDEPENDENT elsewhere branch — the read explained by
//!    the genome outside this locale at the derived background rate
//!    (rule: the read matches its origin at the MEASURED per-base rate
//!    — the census's zero-unmatched-bp measurement over 1,063,522
//!    verified occurrences — under the max-entropy uniform origin
//!    over the donor's diploid genome, the same max-entropy spirit as
//!    the old Poisson model's beta = M/|U|: the unit of match
//!    likelihood spread uniformly over the universe of placements;
//!    the genome size measured from the panel, the instrument's own
//!    truth-free genome universe). THE FLOOR for unexplained reads is
//!    E(read), not the all-mismatch 150*B (which survives only as the
//!    per-placement lower bound inside the local branch's logsumexp,
//!    a bit-exact no-op wherever an anchor pins). A read whose local
//!    placements score no better than E has ~zero discriminating
//!    power between candidates — it does not vote beyond E; the
//!    logsumexp is MONOTONE in the local branch, so a candidate that
//!    spells a read strictly better than every rival under the local
//!    branch still wins that read after the marginalization (the
//!    L4-control preservation property, unit-proven).
//!
//!    THE POCKET-READ ANATOMY (measured first, at chrI L7 — the
//!    form-deciding gate): the winner's pocket placements are 100%
//!    collinear (one continuous path through the pocket row's spell)
//!    with m mean 148.7/150 — the CONTINUOUS world, so no separate
//!    placement-continuity rule is imposed; the marginalization
//!    carries the whole weight. The anatomy's second finding is
//!    recorded in the receipts and the docs note: the pocket reads'
//!    FULL-READ donor occurrences sit on the sample's own S288C/SK1
//!    chrI paths — 66% INSIDE the window's own truth rows, 27% in the
//!    neighboring seam — i.e. the reads are (mostly) locally
//!    generated, and the truth's failure to place them is the pin
//!    skeleton's stored-walk frame blindness (slice A's binding
//!    lesson: rc-frame-qualified anchors are absent from the stored
//!    path walk), NOT an inability of the truth's rows to spell them.
//!    The marginalization bounds that artifact's cost at E; the frame
//!    repair itself is named as the next lever.
//!
//! 8. THE EXTENT NORMALIZATION (stage 2 of the owner's approved
//!    four; env IMPG_REALIGN_TERRITORY_NORMALIZED, beside the
//!    committed convention). THE MEASURED DEFECT: the candidates of
//!    one locus are rows of the same partition whose OWN extents
//!    differ by thousands of bases (the census's window territory is
//!    the PANGENOME territory — every path's full member row — so
//!    every candidate sits wholly "in territory" while a longer row
//!    can place reads a shorter one cannot: the likelihood rewards
//!    whoever happens to SPAN more, not who explains the locus
//!    better; the causal autopsy measured this unplaced-differential
//!    term winner-favored at 582/584 non-rank-1 loci against
//!    truth-favoring both-placed evidence). THE DERIVED RULE, no
//!    tuning constants: candidates are compared on EQUAL extents —
//!    the locus's TERRITORY, its own coordinate, the window's axis
//!    row. Per fold, the territory IMAGE is the interval spanned by
//!    the fold's walk steps shared with the axis row's walk (the
//!    alignment-induced correspondence — the same node-identity
//!    structure the pin machinery already trusts; interior diverged
//!    pockets lie between shared anchors and are INCLUDED). The
//!    local-spell branch is restricted to the image: a placement is
//!    valid iff the read's pinned windows and voted bases project
//!    inside it, and the uniform placement prior counts the image's
//!    placements (2 * (image_len - 149), both strands). Mass outside
//!    the image is handled by the derived symmetric rule: it enters
//!    as the E-branch background for BOTH candidates equally — the
//!    elsewhere branch is candidate-independent, so out-of-image mass
//!    contributes ZERO differential and cannot vote for the spanner.
//!    THE HONEST CHECK: the territory is the locus's own coordinate,
//!    not the truth's and not any candidate's; a true
//!    deletion/insertion haplotype keeps every bit of its IN-WINDOW
//!    evidence (the image is the interval between the fold's outermost
//!    shared anchors — all in-window material, diverged pockets
//!    included), and its material beyond the window is E for every
//!    candidate equally, so it loses no RELATIVE evidence and no rival
//!    can harvest it. A fold sharing no anchor with the axis row has
//!    an EMPTY image: it spells no territory material of this locus
//!    and every unit falls to E under it. Under the default (env
//!    unset) every fold's image is its full extent and the rule
//!    degenerates EXACTLY to the committed comparison — the identity
//!    gate.
//!
//!    THE EDGE REFINEMENT (measured at the two boundary cases the
//!    first fleet pass caught — chrV L24, the one regression, named
//!    and diagnosed): the interval between the outermost shared
//!    anchors EXCLUDES variant pockets at the territory's EDGE — a
//!    true edge variant (no shared anchor on its outer side) cannot
//!    be bounded by an anchor and fell to E, silencing the truth's
//!    own edge material. THE DERIVED WITNESS that separates the two
//!    edge populations (measured): the fold's edge-pocket steps are
//!    checked against EVERY OTHER WINDOW's own axis-row node set
//!    (the component's own walks, the panel's committed structure):
//!    edge material CLAIMED by another window's axis row is that
//!    window's territory (excluded — chrMT L2's winner row head:
//!    131 of 158 steps shared with the NEIGHBOR window's axis walk,
//!    the orphaned-mass harvest); edge material claimed by NO
//!    window's walk is this locus's own edge variant, in-territory
//!    BY ELIMINATION (included, bounded by the row's own extent —
//!    chrV L24's SK1 ortholog head: 12 steps shared with no window's
//!    walk, a true variant). The image therefore extends to the
//!    fold's own row edge on each side whose pocket no other window
//!    claims; a fold with no shared anchor at all is own-material
//!    iff no other window claims any of its steps. The window's own
//!    axis fold IS the territory by definition: its image is its
//!    full extent. No tuning constants — the attribution is exact
//!    set membership over the committed axis walks.
//!
//! 9. THE ALIGNMENT-INDUCED LOCALITY DOMAIN (the owner's go; env
//!    IMPG_REALIGN_LOCALITY_DOMAIN, beside the committed
//!    conventions; implies rule 8). THE DOMAIN CHANGE: the maps'
//!    member lists are the LOCALITY FOLDS (built receipt-side by
//!    locality-domain-maps.py): the universe = the connected
//!    component of the window's axis partition pggb graph that
//!    contains the axis row (partition membership over-joins
//!    disconnected material — those rows leave the locus's
//!    universe); the fold material on each universe path = the
//!    path-interval hull of the MAXIMAL MONOTONE AXIS-
//!    CORRESPONDENCE — order-preserving (path position, axis
//!    position) pairs pooled from the pggb shared segments (the
//!    within-partition real alignment, colinear LIS core) and the
//!    interned syng nodes shared with the axis row's walk, chained
//!    along the path's own tiling across partition boundaries only
//!    within the component's own axis partitions and admitted only
//!    at the correspondence's monotone ends (the seam rejoin: the
//!    partition boundary severed one continuous material and the
//!    shared-node coverage names exactly the severed piece). The
//!    records re-derive per locality: an occurrence votes iff its
//!    interval overlaps a locality fold row on the occurrence's own
//!    path (the census's own territory-touch convention on the new
//!    extents). The territory images become exact under this fold
//!    construction (the fold IS the correspondence core), so the
//!    locality mode implies the territory rule; all other machinery
//!    is unchanged.

use clap::Args;
use crate::genome_inference::{panel_routes as routes, read_json, territory_frame, PanelIdentity};
use crate::syng::{SyncmerParams, SyngIndex};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// The reads' exact length: the realignment model's placement
/// arithmetic, the per-read mismatch log-probability base count, and
/// the derive cache's read records are keyed to it. The yeast
/// validation's committed reads are 150 (the default, behavior
/// identical); the HG002 locus pilot reads are 148 (the GIAB BAM's
/// uniformly trimmed length). Derived from the derive cache's own
/// uniform read length at load (the cache is a pure function of the
/// FASTQ and the panel).
/// THE FRESH-PROCESS READ LENGTH (the committed initial value of
/// the read length in every standalone run): every committed receipt
/// computes the elsewhere branch E BEFORE the derive cache states the
/// cache's own declared read length, so E is measured at THIS value
/// in every standalone run.
/// THE CONCURRENT-CALLS AMENDMENT: the read length is now PER-CALL
/// state (a local in `run_injected`), not a process-global — under
/// the combined process's CONCURRENT calls a process-global would
/// race across calls (one call's derive-cache store landing before
/// another call's E read, the exact divergence class the serial
/// combined process's first gate receipted at the stale 148). The
/// per-call local reproduces the fresh-process trajectory BY
/// CONSTRUCTION at any call concurrency: initialized here, read at E
/// before the derive load stores the declared length, and shared
/// immutably with the seam threads afterwards.
const READ_LENGTH_FRESH: usize = 150;
#[inline]
fn read_len() -> usize {
    READ_LENGTH_FRESH
}

// ---------------------------------------------------------------------------
// Options and small utilities (the slice A runner conventions).
// ---------------------------------------------------------------------------

/// The product model identity (the CLI output manifest of the
/// per-locus diplotype calls).
pub const MODEL: &str = "anchor-realign-diplotype-calls-v1";

/// The selection convention. `Env` (the instrument default) defers to
/// the committed environment switches verbatim; the other values pin
/// the convention explicitly (the product CLI default is `Locality` —
/// the alignment-induced locality domain, the current best measured
/// selection: realignment likelihood + territory-image extents +
/// locality folds).
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Selection {
    /// The committed environment switches decide (IMPG_REALIGN_*)
    Env,
    /// The committed default convention (no territory/locality rules)
    Default,
    /// The territory-image extent normalization (rule 8)
    Territory,
    /// The alignment-induced locality domain (rule 9; implies rule 8)
    Locality,
}

#[derive(Args, Debug)]
pub struct Options {
    /// The panel syng prefix.
    #[arg(long)]
    pub panel: String,
    /// The route graph directory (graph.json + sources).
    #[arg(long)]
    pub routes: PathBuf,
    /// The partition-embedded graph directory (partition<N>.gfa.map.json).
    #[arg(long)]
    pub partition_graphs: PathBuf,
    /// The component lane (e.g. S288C#0#chrI) — the axis path whose
    /// member rows define the window-to-partition mapping.
    #[arg(long)]
    pub component: String,
    /// The committed multi-matching census receipt (records, placements).
    #[arg(long)]
    pub census: PathBuf,
    /// The sample FASTQ (read identity + the quality-byte assertion).
    #[arg(long)]
    pub reads: PathBuf,
    /// The record-derivation cache (the slice A binary; a pure
    /// function of the FASTQ and the panel).
    #[arg(long)]
    pub derive_cache: PathBuf,
    /// The scored loci (window ids of the component, comma-separated,
    /// e.g. 2,4,7) — the pilot set.
    #[arg(long)]
    pub loci: String,
    /// The likelihood receipt path (sidecars are derived from it).
    #[arg(long)]
    pub out: PathBuf,
    /// The exactness-sample receipt path (vote-level detail for the
    /// independent checker).
    #[arg(long)]
    pub exactness_out: PathBuf,
    /// The evenly-spaced exactness sample count per locus (a stated
    /// measurement sample; with IMPG_SCORE_EXACTNESS_ASSERT=1 the
    /// in-process assertion covers the whole pilot domain as well —
    /// the receipt/gate form of the removed always-on check).
    #[arg(long, default_value_t = 100)]
    pub exactness_sample: u64,
    /// The pocket-read anatomy sidecar (slice C: per unit, classified
    /// by whether the winner and truth pairs place it, the placement
    /// anatomy on the winner and truth folds plus the read's donor-
    /// path full-read verification). Pure emission; no scoring change.
    #[arg(long)]
    pub anatomy_out: Option<PathBuf>,
    /// The pin-skeleton sidecar (slice D: per fold, the canonical-
    /// scheme steps with per-step frame tags beside the stored path
    /// walk's steps and the position diff — the frame-repair audit).
    #[arg(long)]
    pub skeleton_out: Option<PathBuf>,
    /// The frame-audit-only mode (slice D's diagnosis): emit the
    /// per-row skeleton diff (stored walk vs canonical scheme) for
    /// every axis-partition row of the given loci and exit — no
    /// reads, no census, no scoring. The output path is --out.
    #[arg(long, default_value_t = false)]
    pub frame_audit_only: bool,
    /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard).
    #[arg(long, default_value_t = 64.0)]
    pub rss_budget_gib: f64,
    /// The selection convention. The instrument default is Env (the
    /// committed environment switches decide); the product CLI pins
    /// Locality (the alignment-induced locality domain).
    #[arg(long, value_enum, default_value_t = Selection::Env)]
    pub selection: Selection,
    /// THE PRODUCT EMISSION (the owner's product swap): per-locus
    /// diplotype calls — the two called paths (the winner fold pair,
    /// spelled as member rows), the class identity, the model-internal
    /// QUAL (the cluster machinery's Phred form over the likelihood
    /// separation) and the E/alternative-gap fields — truth-free, one
    /// JSON record per locus, stream-written beside the receipt. The
    /// instrument receipt (--out) is unchanged when this is unset.
    #[arg(long)]
    pub product_calls: Option<PathBuf>,
    /// THE DIPLOTYPE SEQUENCE EMISSION (the write-up codification,
    /// the owner's ruling: the diplotype sequence is a PRODUCT OUTPUT
    /// at every called locus): per-locus hap1/hap2 records — the
    /// chosen pair's walks from the local graph spelled into the
    /// reconstructed diploid sequence, each record carrying the
    /// model-internal QUAL, the called-class identity and the tied
    /// classes' coalesced called set (every called fold's sequence,
    /// the owner's ambiguity ruling). Requires --product-calls; emits
    /// <product>.sequences.jsonl beside the product calls. Truth-free
    /// (the called material is the product run's own; the called-vs-
    /// truth sequence QV lives only behind --truth-qv).
    #[arg(long, requires = "product_calls", default_value_t = false)]
    pub product_sequences: bool,
    /// THE MASKING OPERATOR (the held-out experiment's stranger filter,
    /// folded onto the call path; the region-as-interface flag set):
    /// candidate rows whose member path names match any of these
    /// fnmatch-style globs (e.g. 'HG002#*') leave the candidate universe
    /// - loudly counted per locus in the product record and the receipt.
    /// THE AXIS ROW (the --component lane) is never excluded: it is the
    /// coordinate frame. Default (empty): no row is excluded, the
    /// committed candidate universe verbatim.
    #[arg(long, num_args = 1.., value_delimiter = ',')]
    pub exclude_paths: Vec<String>,
    /// THE REGION-AS-INTERFACE FILTER (the owner's ruling: "detect
    /// any locus, any position, any allele"): a region `chr:start-end`
    /// (0-based half-open, the component's own reference frame, e.g.
    /// chr1:4700000-4800000) resolved against the component's axis
    /// rows to the covered windows - the scored loci, in place of
    /// --loci. Conflicts with --loci.
    #[arg(long, conflicts_with = "loci")]
    pub region: Option<String>,
    /// THE FLANK-OVERLAP CANDIDATE FILTER (the owner's design, the
    /// presence semantics ruled at the measured fork): N in
    /// bp (`--flank-overlap 10000`) or a percent of the window length
    /// (`--flank-overlap 10%` -> ceil(pct x length)). A member row is
    /// KEPT in a window's candidate universe only if it carries
    /// aligned material over BOTH FLANK intervals of the window's
    /// axis extent - the first N bp and the last N bp, in the focal
    /// reference's own coordinate frame; a window shorter than N
    /// degenerates to simple presence. The axis row is never
    /// excluded (the --exclude-paths precedent).
    /// Default (unset): no row is excluded - the committed candidate
    /// universe verbatim, the committed receipts byte-identical.
    #[arg(long)]
    pub flank_overlap: Option<FlankOverlap>,
    /// THE GLOBAL-PANEL FRAME (the campaign's lever 1, the
    /// global-panel call mode: docs/syng-gmem-bwt/global-panel-call-mode.md):
    /// --panel names the ONE global panel and --partition-graphs the
    /// territory-local substrate (member rows "path:start-end",
    /// slice-relative coordinates — the committed locality maps).
    /// When set, the maps' member rows and the component lane are
    /// translated to the global frame at load (the record header IS
    /// the frame statement — the re-anchor principle), the machinery
    /// runs against the global grid, and the product emission is
    /// translated back to the partition's own frame (the committed
    /// reference-anchored output shape). Default (unset): the
    /// per-territory frame verbatim — the committed path byte-for-byte.
    #[arg(long, default_value_t = false)]
    pub global_frame: bool,
    /// THE RECORD-ORIGIN CONVENTION (the global-panel call mode's
    /// content-origin parameter): the committed substrate's locus.fa
    /// content need not sit at its header's own interval — the L1-L4
    /// substrates' content sits ONE BASE LATER (d=1, the uniform
    /// delta=1 AGC-identity receipts), the pilot's at the header
    /// interval exactly (d=0). The global frame, which deletes
    /// locus.fa by design, must be TOLD the convention: the maps'
    /// member rows and the component lane translate at the content
    /// origin `start + d` (the frame table likewise), and the
    /// emission's back-translation cancels d exactly (the product
    /// output is origin-independent). Default 0 = the header's own
    /// interval. Fail-closed: a nonzero origin WITHOUT --global-frame
    /// is refused (the committed local frame is closed — it never
    /// needed the absolute truth).
    #[arg(long, default_value_t = 0, requires = "global_frame")]
    pub record_origin: u64,
    /// THE FAR-BIN PRUNE (the owner's own routing principle, the
    /// §10.6.17-named census lever, enforced at the call's census view:
    /// "we used to decide where to put each read, so it's basically
    /// being routed by its best MEMs. If they go into different bins
    /// unless they're basically adjacent, we don't want to be thinking
    /// about them."). A census record's own routing keeps it at the
    /// windows its own walk placements touch — one window, or two
    /// ADJACENT windows when the anchors straddle a boundary (both by
    /// pure interval overlap; no tolerance, nothing tuned). A pattern
    /// whose placements scatter into NON-ADJACENT windows — the
    /// measured far-bin class, e.g. L1's one 43,686-multiplicity
    /// single-anchor pattern at windows [12,16] carrying 31.1% of the
    /// restricted census mass — is the "different bins, not basically
    /// adjacent" class: there is no choice-free home bin to keep it at,
    /// so the rule refuses to think about it ANYWHERE and the call
    /// prunes it from its census view entirely (bind, variants, and
    /// every locus's scoring). The test is pure window-set geometry
    /// (the touched-window diameter against the census's own
    /// pair-bins adjacency constant 1) — never a tuned threshold.
    /// DEFAULT: ON at the global frame (the owner's ruling on the
    /// §10.6.18 gate receipts: zero truth losses at L1-L4, the
    /// pilot's four truth WINS); the opt-out form is
    /// `--far-bin-prune=false` (the committed unpruned census view,
    /// byte-identical to the pre-ruling receipts). The LOCAL frame
    /// is byte-frozen and refuses the prune loudly (fail-closed in
    /// run_injected, naming the opt-out).
    #[arg(
        long,
        value_parser = clap::value_parser!(bool),
        default_value_t = true,
        default_missing_value = "true",
        num_args = 0..=1
    )]
    pub far_bin_prune: bool,
    /// THE TEST-MODE FLAG (the owner's ruling, locked in: the
    /// called-vs-truth sequence QV is a VERY SPECIAL TESTING THING
    /// and lives ONLY here): the truth haplotypes as panel path
    /// names (e.g. S288C#0#chrMT SK1#0#chrMT). Requires --product-calls;
    /// emits the per-locus truth-referenced sequence QV (the
    /// committed convention: gap-affine End2End 0/4/6/2, error =
    /// (mismatches + gap columns) / total columns, QV = -10 log10,
    /// the QV-60 perfect-call cap) to a SEPARATE test artifact; never
    /// mixed into the product output.
    #[arg(long, num_args = 1.., requires = "product_calls")]
    pub truth_qv: Vec<String>,
    /// The pinned biWFA alignment helper for the --truth-qv test mode
    /// (the committed receipt-side binary; gap-affine End2End
    /// 0/4/6/2, Medium memory, no heuristic — the same machinery the
    /// validated QV receipts used).
    #[arg(long, default_value = "genome/instrumented/qv-biwfa/target/release/qv-biwfa")]
    pub qv_biwfa: PathBuf,
}

// ---------------------------------------------------------------------------
// The far-bin prune's frame gate (see run_injected): the LOCAL frame is
// byte-frozen; the prune (now the product default) never applies there —
// refuse loudly, naming the opt-out form.
// ---------------------------------------------------------------------------
fn far_bin_prune_frame_gate(global_frame: bool, far_bin_prune: bool) -> io::Result<()> {
    ensure(
        !far_bin_prune || global_frame,
        "--far-bin-prune requires --global-frame (the committed local frame \
         is byte-frozen; state --far-bin-prune=false for the committed local \
         form)",
    )
}

// ---------------------------------------------------------------------------
// THE TRUTH-QV TEST MODE (the owner's ruling, locked in: the
// called-vs-truth sequence QV is a VERY SPECIAL TESTING THING — it
// lives ONLY behind the explicit --truth-qv flag, which REQUIRES the
// truth haplotypes as input; it is never part of the default product
// output and never mixed into it). The convention is the committed
// QV stage's exactly (realign-sequence-qv.py of record): the two
// called homologs are matched to the two truth homologs by best
// assignment WITHOUT replacement (both orders scored by the
// assignment-wide per-base error rate — the sum of the two pairs'
// edits over the sum of the two pairs' alignment columns; lower
// rate wins, ties to the lower total edit count, then the direct
// order); each pair is aligned gap-affine End2End (match 0 /
// mismatch 4 / gap-open 6 / gap-extend 2, Medium memory, no
// heuristic) through the SAME pinned biWFA helper binary the
// validated receipts used (protocol: one "seqA\tseqB" request line
// per alignment, one "matches\tmismatches\tins\tdels" response);
// per-base error = (mismatches + ins + dels) / total columns; QV =
// -10 log10(error) with the QV-60 cap at error <= 1e-6 (a perfect
// call, zero edits). The material per fold is the fold's own spelled
// sequence — the panel fetch of the member row [start, end) window,
// the same convention-(a) material of the committed stage.
//
// BESIDE it, env-gated (IMPG_TRUTH_QV_ALLWAVE=1), THE OFFICIAL
// SEQUENCE-QV CONVENTION (the owner's ruling: the sequence QV is
// quantity 1, called-vs-truth sequence accuracy, in the likegt
// max-qv form EXACTLY — allwave aligns the truth haplotypes against
// the called haplotypes with likegt's own invocation `-i <fasta> -t
// 4 -p tree:5:0:0`, identity = matches / query_span per PAF line,
// QV = -10 log10(1 - identity) capped at 50, per-truth-hap BEST
// alignment among the called homologs, the diplotype QV the mean of
// the two per-hap QVs) — emitted as additive official_allwave fields
// in the same test artifact, the committed biWFA fields untouched
// and comparable. The env-unset default path is byte-identical to
// the committed receipts.
// ---------------------------------------------------------------------------

/// One persistent biWFA helper child (stdin/stdout pipes; one
/// request-response round trip per alignment). The child is closed
/// EXPLICITLY by `finish` (close our stdin end so the helper sees EOF
/// and exits, then reap it); the field-drop order alone would leave a
/// wait() deadlock (the helper blocks on stdin while we block on its
/// exit) — so no Drop impl waits, and run() always calls finish().
pub struct QvAligner {
    child: std::process::Child,
    stdin: Option<std::process::ChildStdin>,
    stdout: std::io::BufReader<std::process::ChildStdout>,
}

impl QvAligner {
    pub fn spawn(helper: &std::path::Path) -> io::Result<Self> {
        let mut child = std::process::Command::new(helper)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|error| {
                invalid(&format!(
                    "spawning the qv-biwfa helper {helper:?}: {error}"
                ))
            })?;
        let stdin = child
            .stdin
            .take()
            .ok_or_else(|| invalid("the qv-biwfa helper opened no stdin"))?;
        let stdout = std::io::BufReader::new(
            child
                .stdout
                .take()
                .ok_or_else(|| invalid("the qv-biwfa helper opened no stdout"))?,
        );
        Ok(QvAligner {
            child,
            stdin: Some(stdin),
            stdout,
        })
    }

    /// One alignment: the two sequences' (matches, mismatches, ins,
    /// dels) per the helper's walked-CIGAR protocol.
    pub fn align(&mut self, first: &[u8], second: &[u8]) -> io::Result<[u64; 4]> {
        ensure(
            !first.is_empty() && !second.is_empty(),
            "empty sequence in the truth-QV test mode",
        )?;
        let stdin = self
            .stdin
            .as_mut()
            .ok_or_else(|| invalid("the qv-biwfa helper is already closed"))?;
        stdin.write_all(first)?;
        stdin.write_all(b"\t")?;
        stdin.write_all(second)?;
        stdin.write_all(b"\n")?;
        stdin.flush()?;
        let mut line = String::new();
        self.stdout.read_line(&mut line)?;
        let line = line.trim_end();
        if let Some(reason) = line.strip_prefix("FAIL:") {
            return Err(invalid(&format!(
                "the qv-biwfa helper failed an alignment: {reason}"
            )));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        ensure(fields.len() == 4, "malformed qv-biwfa helper response")?;
        let mut counts = [0u64; 4];
        for (slot, field) in counts.iter_mut().zip(fields) {
            *slot = field
                .parse()
                .map_err(|_| invalid("non-numeric qv-biwfa helper response"))?;
        }
        Ok(counts)
    }

    /// Close the request stream (the helper exits at EOF) and reap
    /// the child.
    pub fn finish(&mut self) -> io::Result<()> {
        if let Some(mut stdin) = self.stdin.take() {
            stdin.flush()?;
        }
        self.child
            .wait()
            .map(|_| ())
            .map_err(|error| invalid(&format!("the qv-biwfa helper wait failed: {error}")))
    }
}

// ---------------------------------------------------------------------------
// THE OFFICIAL SEQUENCE-QV CONVENTION (the owner's ruling: the
// sequence QV is quantity 1, called-vs-truth sequence accuracy,
// computed by the likegt max-qv form EXACTLY — surveyed from
// likegt's own implementation, /home/erikg/likegt's `max-qv` command
// (src/commands/max_qv.rs, the README's allwave aligner), not
// invented here). THE CONVENTION: allwave aligns the TRUTH
// haplotypes (the queries) against the CALLED haplotypes (the
// targets) with likegt's own invocation — `allwave -i <fasta> -t
// <threads> -p tree:5:0:0` (the max-qv command's default
// sparsification; allwave's OWN default dual-piece gap-affine scores
// 0,5,8,2,24,1, no override) — then, per the max-qv PAF derivation,
// identity = matches / query_span (query_span = query_end -
// query_start: the per-base error is measured in TRUTH sequence
// space over the ALIGNED SPAN, no forced End2End overhang), QV =
// -10·log10(1 - identity) capped at 50 (a perfect alignment carries
// QV 50.0 exactly); per TRUTH haplotype the BEST alignment among
// the called homologs is taken (likegt's most-similar-haplotype
// selection — the max-qv "best" rule; alignment lines from one
// query to one target pool by (matches, query_span) sums, likegt's
// per-target aggregation); and the diplotype sequence QV is the
// MEAN of the two per-truth-hap QVs (likegt's avg_qv). This runs
// BESIDE the committed biWFA End2End injective-assignment
// computation, env-gated (IMPG_TRUTH_QV_ALLWAVE=1;
// IMPG_TRUTH_QV_ALLWAVE_BIN overrides the allwave binary, default
// "allwave" from PATH) — the committed fields stay byte-identical
// and comparable.
// ---------------------------------------------------------------------------

/// One allwave state for the env-gated official-convention path: the
/// binary and likegt's max-qv thread default (4).
pub struct AllwaveQv {
    binary: std::ffi::OsString,
    threads: usize,
}

impl AllwaveQv {
    /// The likegt max-qv QV derivation from an alignment identity:
    /// QV = -10·log10(1 - identity), clamped to [0, 50]; a perfect
    /// alignment (identity >= 1) is QV 50.0 exactly (likegt's own
    /// cap-and-floor, max_qv.rs parse_paf_line).
    fn qv_of_identity(identity: f64) -> f64 {
        if identity >= 1.0 {
            50.0
        } else {
            let error = 1.0 - identity;
            if error <= 0.0 {
                50.0
            } else {
                (-10.0 * error.log10()).clamp(0.0, 50.0)
            }
        }
    }

    /// Parse allwave's PAF stdout into (query, target, matches,
    /// query_span) tuples — likegt's own fields: matches = field 10,
    /// query_span = query_end - query_start (fields 3/4), the aligned
    /// query span in TRUTH sequence space.
    fn parse_paf(stdout: &str) -> io::Result<Vec<(String, String, u64, u64)>> {
        let mut lines = Vec::new();
        for line in stdout.lines() {
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                // allwave's progress chatter never reaches stdout, but
                // a malformed line is an error, not a silent skip.
                return Err(invalid(&format!(
                    "malformed allwave PAF line ({} fields): {line}",
                    fields.len()
                )));
            }
            let query_start: u64 = fields[2]
                .parse()
                .map_err(|_| invalid("non-numeric allwave PAF query start"))?;
            let query_end: u64 = fields[3]
                .parse()
                .map_err(|_| invalid("non-numeric allwave PAF query end"))?;
            let matches: u64 = fields[9]
                .parse()
                .map_err(|_| invalid("non-numeric allwave PAF matches"))?;
            lines.push((
                fields[0].to_string(),
                fields[5].to_string(),
                matches,
                query_end.saturating_sub(query_start),
            ));
        }
        Ok(lines)
    }

    /// One truth haplotype's best alignment among the called
    /// homologs, in likegt's own form: the (query, target) pairs pool
    /// their (matches, query_span) sums (max_qv.rs's per-target
    /// aggregation), the best target is the highest pooled identity
    /// (ties to the higher matches, then the higher query_span, then
    /// the first called homolog). Returns (identity, matches,
    /// query_span, best target name) — identity 0.0 when allwave
    /// emitted no alignment for this truth haplotype.
    fn best_for_hap(
        paf: &[(String, String, u64, u64)],
        truth: &str,
        called: [&str; 2],
    ) -> (f64, u64, u64, Option<String>) {
        let mut best: Option<(u64, u64, usize)> = None; // (matches, span, called index)
        for (index, target) in called.iter().enumerate() {
            let (mut matches, mut span) = (0u64, 0u64);
            for (query, t, m, s) in paf {
                if query == truth && t == target {
                    matches += m;
                    span += s;
                }
            }
            if span == 0 {
                continue; // no alignment line for this pair
            }
            let identity = matches as f64 / span as f64;
            let better = match best {
                None => true,
                Some((b_m, b_s, _)) => {
                    let b_identity = b_m as f64 / b_s as f64;
                    (identity, matches, span) > (b_identity, b_m, b_s)
                }
            };
            if better {
                best = Some((matches, span, index));
            }
        }
        match best {
            None => (0.0, 0, 0, None),
            Some((matches, span, index)) => (
                matches as f64 / span as f64,
                matches,
                span,
                Some(called[index].to_string()),
            ),
        }
    }

    /// One locus under the official convention: the two truth
    /// haplotypes' best alignments among the two called homologs, the
    /// per-hap QVs, and the diplotype QV as their mean (likegt's
    /// avg_qv form).
    fn locus_qv(
        &self,
        locus: u32,
        called0: &[u8],
        called1: &[u8],
        truth0: &[u8],
        truth1: &[u8],
    ) -> io::Result<serde_json::Value> {
        ensure(
            !called0.is_empty()
                && !called1.is_empty()
                && !truth0.is_empty()
                && !truth1.is_empty(),
            "empty sequence in the official-convention truth-QV path",
        )?;
        let fasta = std::env::temp_dir().join(format!(
            "impg-truth-qv-allwave-{}-{locus}.fa",
            std::process::id()
        ));
        {
            use std::io::Write;
            let mut out = std::fs::File::create(&fasta)?;
            for (name, seq) in [
                ("called0", called0),
                ("called1", called1),
                ("truth0", truth0),
                ("truth1", truth1),
            ] {
                write!(out, ">{name}\n")?;
                out.write_all(seq)?;
                out.write_all(b"\n")?;
            }
        }
        let output = std::process::Command::new(&self.binary)
            .arg("-i")
            .arg(&fasta)
            .arg("-t")
            .arg(self.threads.to_string())
            .arg("-p")
            .arg("tree:5:0:0")
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .output()
            .map_err(|error| {
                invalid(&format!("spawning allwave {binary:?}: {error}", binary = self.binary))
            })?;
        let _ = std::fs::remove_file(&fasta);
        ensure(
            output.status.success(),
            "allwave failed in the official-convention truth-QV path",
        )?;
        let stdout = String::from_utf8_lossy(&output.stdout);
        let paf = Self::parse_paf(&stdout)?;
        let called = ["called0", "called1"];
        let (id0, m0, s0, t0) = Self::best_for_hap(&paf, "truth0", called);
        let (id1, m1, s1, t1) = Self::best_for_hap(&paf, "truth1", called);
        let (qv0, qv1) = (Self::qv_of_identity(id0), Self::qv_of_identity(id1));
        let pairs = |truth: &str| {
            let mut value = serde_json::Map::new();
            for target in called {
                let (mut matches, mut span) = (0u64, 0u64);
                for (query, t, m, s) in &paf {
                    if query == truth && t == target {
                        matches += m;
                        span += s;
                    }
                }
                value.insert(
                    format!("truth_{truth}_vs_{target}"),
                    if span > 0 {
                        json!({
                            "matches": matches,
                            "query_span": span,
                            "identity": matches as f64 / span as f64,
                        })
                    } else {
                        serde_json::Value::Null
                    },
                );
            }
            serde_json::Value::Object(value)
        };
        Ok(json!({
            "convention": "likegt-max-qv-allwave (the owner's official ruling)",
            "invocation": format!(
                "allwave -i <locus fasta> -t {} -p tree:5:0:0 (likegt's own max-qv default)",
                self.threads
            ),
            "qv": (qv0 + qv1) / 2.0,
            "qv_truth_hap0": qv0,
            "qv_truth_hap1": qv1,
            "identity_truth_hap0": id0,
            "identity_truth_hap1": id1,
            "best_called_hap_truth_hap0": t0,
            "best_called_hap_truth_hap1": t1,
            "perfect": id0 >= 1.0 && id1 >= 1.0,
            "no_alignment_truth_hap0": t0.is_none(),
            "no_alignment_truth_hap1": t1.is_none(),
            "pairs": {
                "truth0": pairs("truth0"),
                "truth1": pairs("truth1"),
            },
            "best_matches": [m0, m1],
            "best_query_spans": [s0, s1],
        }))
    }
}

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, s.to_string())
}

fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

struct RssGuard {
    budget_bytes: Option<u64>,
}

impl RssGuard {    fn new(budget_gib: f64) -> Self {
        RssGuard {
            budget_bytes: (budget_gib > 0.0).then(|| (budget_gib * (1u64 << 30) as f64) as u64),
        }
    }
    fn probe(&self, stage: &str) -> io::Result<u64> {
        let text = std::fs::read_to_string("/proc/self/status")?;
        let field_kb = |prefix: &str| -> u64 {
            text.lines()
                .find(|l| l.starts_with(prefix))
                .and_then(|l| l.split_whitespace().nth(1).and_then(|n| n.parse::<u64>().ok()))
                .unwrap_or(0)
        };
        let rss_kb = field_kb("VmRSS:");
        // VmHWM is the monotone high-water mark: VmRSS samples instants,
        // so a short-lived spike inside a phase (the tie certificate's
        // Value tree is one) can be missed entirely between two probes.
        let hwm_kb = field_kb("VmHWM:");
        if let Some(budget) = self.budget_bytes {
            ensure(
                rss_kb * 1024 <= budget,
                &format!("RSS guard exceeded at {stage}: {rss_kb} kB"),
            )?;
        }
        eprintln!("[score] rss {stage}: {rss_kb} kB peak {hwm_kb} kB");
        Ok(rss_kb)
    }
}

/// THE RECEIPT WRITER WITH THE TIE CERTIFICATE SPLICED IN.
///
/// Byte-identical to `serde_json::to_writer` of the same object with that
/// member set: serde_json's default map is ORDERED (alphabetical - this
/// crate does not enable `preserve_order`), so a member's slot is purely
/// a function of its key, and every other member goes through the very
/// same writer the whole-object call would use. The difference is
/// memory: the certificate streams from the typed lists instead of being
/// materialized as a value tree first.
fn write_report_with_cut_indices(
    writer: &mut impl Write,
    report: &serde_json::Value,
    cut_indices_key: &str,
    lists: &[(usize, CutSet)],
) -> io::Result<()> {
    let serde_json::Value::Object(map) = report else {
        // (Not an object: nothing to splice. The plain writer's own
        // bytes, which is exactly what the caller would have got.)
        return serde_json::to_writer(writer, report).map_err(io::Error::from);
    };
    writer.write_all(b"{")?;
    let mut position = 0usize;
    for (key, value) in map {
        if position > 0 {
            writer.write_all(b",")?;
        }
        position += 1;
        serde_json::to_writer(&mut *writer, key)?;
        writer.write_all(b":")?;
        if key == cut_indices_key {
            // The `json!([index, list])`-equivalent form
            // `[[index,[spectrum...]],...]`, written straight from the
            // typed list (usize's decimal rendering IS serde_json's own
            // integer form).
            writer.write_all(b"[")?;
            for (list_position, (index, list)) in lists.iter().enumerate() {
                if list_position > 0 {
                    writer.write_all(b",")?;
                }
                write!(writer, "[{},[", index)?;
                for (entry_position, spectrum_index) in list.iter().enumerate() {
                    if entry_position > 0 {
                        writer.write_all(b",")?;
                    }
                    write!(writer, "{}", spectrum_index)?;
                }
                writer.write_all(b"]]")?;
            }
            writer.write_all(b"]")?;
        } else {
            serde_json::to_writer(&mut *writer, value)?;
        }
    }
    writer.write_all(b"}")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Receipt inputs (slice A structures, verbatim).
// ---------------------------------------------------------------------------
/// One verified occurrence in the committed multi-matching census.
#[derive(Deserialize)]
struct CensusOccurrence {
    path: usize,
    start: u64,
    orientation: u8,
    /// The touched component windows (the routing's own territory-touch
    /// rule; these are WINDOW ids, not partition ids).
    #[serde(default)]
    partitions: Vec<u32>,
    /// Per read-matched interval, the contained covered node abs ids in
    /// path order (the committed containment rule).
    #[serde(default)]
    intervals: Vec<Vec<u32>>,
}

#[derive(Deserialize)]
struct CensusRecordLine {
    record: u64,
    multiplicity: u64,
    #[serde(default)]
    anchors: u64,
    #[serde(default)]
    occurrences: Vec<CensusOccurrence>,
}

/// THE OWNER'S FAR-BIN RULE over one census record (the routing
/// principle the §10.6.17 autopsy named, now the prune's own
/// geometry): "we used to decide where to put each read, so it's
/// basically being routed by its best MEMs. If they go into
/// different bins unless they're basically adjacent, we don't want
/// to be thinking about them." A record's own routing keeps it at
/// the windows its own walk placements touch — ONE window, or TWO
/// ADJACENT windows when the anchors straddle a boundary (both by
/// pure interval overlap, the routing's `partitions_at` rule: no
/// tolerance, nothing to tune). A pattern whose placements scatter
/// into NON-ADJACENT windows — the measured far-bin class, e.g.
/// L1's one 43,686-multiplicity single-anchor pattern routed to
/// windows [12,16] — is exactly the owner's "different bins, not
/// basically adjacent" class: there is no choice-free home bin to
/// keep it at, so the rule refuses to think about it anywhere. The
/// test is pure window-set geometry: the record's touched-window
/// DIAMETER against the adjacency constant 1 (adjacent-or-equal) —
/// the census's own pair-bins taxonomy constant, never a tuned
/// threshold. A record touching no window is not scattered (it
/// votes nowhere by its own routing; the pilot rule already
/// excludes it).
fn census_record_far_bin_scattered(line: &CensusRecordLine) -> bool {
    let mut lo = u32::MAX;
    let mut hi = u32::MIN;
    for occ in &line.occurrences {
        for &window in &occ.partitions {
            lo = lo.min(window);
            hi = hi.max(window);
        }
    }
    lo != u32::MAX && hi - lo > 1
}

/// One member row of a partition graph map (the alignment-induced
/// membership; the partition BEDs are source-forward).
#[derive(Clone, Debug, Deserialize, PartialEq)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

#[derive(Clone, Deserialize)]
struct PartitionMap {
    partition: u32,
    members: Vec<MemberRow>,
}

// ---------------------------------------------------------------------------
// The projection core (slice A's committed machinery, verbatim).
// ---------------------------------------------------------------------------

/// The canonical decoded walk of a record's token pattern: anchors
/// (signed node, canonical position), first anchor at position 0.
fn decode_record_tokens(tokens: &[u64]) -> io::Result<Vec<(i32, u64)>> {
    ensure(!tokens.is_empty() && tokens.len() % 2 == 1, "invalid record tokens")?;
    let mut anchors = Vec::with_capacity(tokens.len() / 2 + 1);
    let mut position = 0u64;
    for (i, &token) in tokens.iter().enumerate() {
        if i % 2 == 0 {
            let zigzag = token
                .checked_sub(2)
                .ok_or_else(|| invalid("record node token"))?
                / 2;
            let node = ((zigzag >> 1) as i64 ^ -(zigzag as i64 & 1)) as i32;
            anchors.push((node, position));
        } else {
            let gap = token
                .checked_sub(1)
                .ok_or_else(|| invalid("record gap token"))?
                / 2;
            position += gap;
        }
    }
    Ok(anchors)
}

/// The record's span: the last anchor's window end in canonical
/// coordinates (the pattern's full hull length).
fn walk_span(walk: &[(i32, u64)], k: u64) -> u64 {
    walk.last().map(|&(_, p)| p + k).unwrap_or(k)
}

/// The read's OWN anchor positions for one mirror state, per canonical
/// anchor index: an unmirrored read spells the canonical pattern (own
/// anchor j at w_lo + c_j); a mirrored read spells its reverse
/// complement (canonical anchor j at w_lo + (c_last - c_j)).
fn own_positions_of(canonical: &[(i32, u64)], mirror: bool, w_lo: u64) -> Vec<u64> {
    let last = canonical.last().map(|&(_, p)| p).unwrap_or(0);
    canonical
        .iter()
        .map(|&(_, pos)| {
            if mirror {
                w_lo + (last - pos)
            } else {
                w_lo + pos
            }
        })
        .collect()
}

/// The read's physical strand at one occurrence: the census orientation
/// refers to the canonical token walk, which is the read's own walk
/// exactly when the read is not mirrored.
fn physical_forward(orientation: u8, mirror: bool) -> bool {
    (orientation == 0) != mirror
}

/// The anchor window start under one placement (slice A, verbatim).
fn anchor_window_start(canonical_pos: u64, start: u64, span: u64, k: u64, orientation: u8) -> u64 {
    if orientation == 0 {
        start + canonical_pos
    } else {
        start + span - k - canonical_pos
    }
}

/// One candidate key's canonical anchor chain with the projection
/// geometry the binding verifies against (slice A).
struct KeyShape {
    canonical: Vec<(i32, u64)>,
    span: u64,
    read: usize,
    mirrored: bool,
    own_positions: Vec<u64>,
}

/// The per-shift crop of a hoisted (union-fetched) occurrence range
/// (phase 1, lever 1): exactly the bytes the pre-rebuild per-shift
/// fetch returned. All three origin shifts share the high bound
/// `start + span + 1` (clamped to the path length by `fetch_seq`
/// before the call) and differ only in the low bound
/// `max(0, start + shift - 1)`, so from the union buffer — fetched
/// once per (occurrence, candidate) starting at panel coordinate
/// `seq_lo` — the per-shift crop is the buffer's tail from the
/// shift's low bound (empty when that low bound reaches past the
/// clamped high end, the empty-crop case of the old fetch).
fn verify_shift_crop<'a>(seq: &'a [u8], seq_lo: u64, occ_start: u64, shift: i64) -> &'a [u8] {
    let lo = (occ_start as i64 + shift - 1).max(0) as u64;
    let from = lo.saturating_sub(seq_lo) as usize;
    seq.get(from..).unwrap_or(&[])
}

/// Sequence verification of one key at one occurrence (slice A,
/// verbatim semantics; phase 1 hoisted the fetch): every anchor
/// window must sit at its projected panel position and equal the
/// representative read's own k-mer. `seq` is the once-per-
/// (occurrence, candidate) union range fetch covering all three
/// origin shifts, starting at panel coordinate `seq_lo`; each shift
/// verifies against the exact slice its own per-shift fetch
/// returned before the rebuild (unit-proven; the identity gate
/// proves the receipts unchanged end-to-end).
fn verify_key_at(
    shape: &KeyShape,
    occ: &CensusOccurrence,
    reads: &[Vec<u8>],
    seq: &[u8],
    seq_lo: u64,
    k: u64,
    shift: i64,
) -> bool {
    let lo = (occ.start as i64 + shift - 1).max(0) as u64;
    let seq = verify_shift_crop(seq, seq_lo, occ.start, shift);
    let origin = occ.start as i64 + shift;
    let read = &reads[shape.read];
    let forward = physical_forward(occ.orientation, shape.mirrored);
    for (j, &(node, canonical_pos)) in shape.canonical.iter().enumerate() {
        let _ = node;
        let window_lo =
            anchor_window_start(canonical_pos, origin as u64, shape.span, k, occ.orientation);
        let rel_lo = window_lo as i64 - lo as i64;
        if rel_lo < 0 || rel_lo + k as i64 > seq.len() as i64 {
            return false;
        }
        let window = &seq[rel_lo as usize..(rel_lo + k as i64) as usize];
        let own_pos = shape.own_positions[j] as usize;
        let kmer = &read[own_pos..own_pos + k as usize];
        let ok = if forward {
            window == kmer
        } else {
            revcomp(window) == kmer
        };
        if !ok {
            return false;
        }
    }
    true
}

/// Merge windows (lo, hi): overlapping and exactly-abutting windows
/// merge; a positive gap splits. Returns the merged intervals.
fn merge_windows(mut windows: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    windows.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(windows.len());
    for (lo, hi) in windows {
        match merged.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => merged.push((lo, hi)),
        }
    }
    merged
}

/// The complement of the merged own-coordinate anchor windows within
/// the read: the SKIPPED base positions (the bases the MEMs do not
/// cover).
fn skipped_positions(own_positions: &[u64], k: u64, read_len: u64) -> Vec<(u64, u64)> {
    let windows: Vec<(u64, u64)> = own_positions.iter().map(|&p| (p, p + k)).collect();
    let merged = merge_windows(windows);
    let mut skipped: Vec<(u64, u64)> = Vec::new();
    let mut cursor = 0u64;
    for (lo, hi) in merged {
        if lo > cursor {
            skipped.push((cursor, lo));
        }
        cursor = cursor.max(hi);
    }
    if cursor < read_len {
        skipped.push((cursor, read_len));
    }
    skipped
}

/// The required panel-path sign of canonical anchor j at one
/// occurrence: orientation 0 spells the pattern (sign as decoded);
/// orientation 1 spells its reverse complement (flipped).
fn required_sign(node: i32, orientation: u8) -> i32 {
    if orientation == 0 {
        node
    } else {
        -node
    }
}

/// The anchor correspondence of one placement against one candidate
/// row's node positions (node abs id -> (window-start bp, sign)):
/// per canonical anchor, the row-side window-start position with the
/// required sign. A node occurring several times admits several
/// monotone assignments (strictly increasing for orientation 0,
/// decreasing for orientation 1, in canonical anchor order); a
/// position is PINNED iff it is identical across ALL monotone
/// assignments, otherwise it is AMBIGUOUS (None, named). Anchors
/// absent from the row are None but not ambiguous.
fn anchor_correspondence(
    canonical: &[(i32, u64)],
    orientation: u8,
    row_positions: &dyn Fn(u32) -> Vec<(u64, i32)>,
) -> (Vec<Option<u64>>, Vec<usize>) {
    let m = canonical.len();
    let mut candidates: Vec<Vec<u64>> = Vec::with_capacity(m);
    for &(node, _) in canonical {
        let required = required_sign(node, orientation);
        let positions = row_positions(required.unsigned_abs());
        candidates.push(
            positions
                .iter()
                .filter(|&&(_, sign)| sign == required.signum())
                .map(|&(bp, _)| bp)
                .collect(),
        );
    }
    let mut pinned: Vec<Option<u64>> = vec![None; m];
    // (The exact feasibility DP — the chrIV repeat-domain repair.
    // The committed enumeration descended every monotone assignment;
    // at chrIV's repeat-dense rows the leaf count is exponential in the
    // anchors-per-repeat (measured: locus 39's axis partition carries
    // rows with a syncmer node at 34 positions, and the serial run
    // spent 760.7s in that ONE locus's scoring phase on exactly this
    // path; sibling partitions carry scaffold rows with nodes at 1356
    // positions). The DP computes the SAME pinned vector exactly:
    // every complete monotone chain contains every anchor that has
    // candidates (the recursion dies otherwise), so a position is
    // pinned iff it is the anchor's UNIQUE feasible candidate —
    // feasible forward (a monotone prefix of the earlier candidate-
    // bearing anchors ends below it) and backward (a suffix of the
    // later ones starts above it); two feasible candidates disagree
    // across leaves (ambiguous, None); an anchor with candidates but
    // NO feasible position means ZERO complete chains, and the
    // enumeration's zero-leaf outcome leaves every pin None — which
    // the DP reproduces because no (anchor, position) is then both
    // forward- and backward-feasible. Anchors with no candidates
    // never enter a chain and stay None, exactly as before. The
    // equivalence is unit-proven against a copy of the committed
    // enumeration over randomized synthetic walks, and the chrMT/chrI
    // receipts are byte-identity-gated against the committed runs.)
    let idx: Vec<usize> = (0..m).filter(|&j| !candidates[j].is_empty()).collect();
    if !idx.is_empty() {
        let increasing = orientation == 0;
        let n = idx.len();
        // forward[k][t]: some monotone prefix over the idx chain
        // 0..=k ends at candidate t of idx[k] (candidates are sorted
        // ascending — the fold walks are).
        let mut forward: Vec<Vec<bool>> = Vec::with_capacity(n);
        for k in 0..n {
            forward.push(vec![false; candidates[idx[k]].len()]);
        }
        for t in 0..candidates[idx[0]].len() {
            forward[0][t] = true;
        }
        for k in 1..n {
            let prev = &candidates[idx[k - 1]];
            // the prefix/suffix-OR of the previous anchor's forward
            let mut below = vec![false; prev.len() + 1];
            if increasing {
                for t in 0..prev.len() {
                    below[t + 1] = below[t] || forward[k - 1][t];
                }
            } else {
                let mut above = vec![false; prev.len() + 1];
                for t in (0..prev.len()).rev() {
                    above[t] = above[t + 1] || forward[k - 1][t];
                }
                below = above;
            }
            let here = &candidates[idx[k]];
            for (t, &bp) in here.iter().enumerate() {
                if increasing {
                    let cut = prev.partition_point(|&p| p < bp);
                    forward[k][t] = below[cut];
                } else {
                    let start = prev.partition_point(|&p| p <= bp);
                    forward[k][t] = below[start];
                }
            }
        }
        // backward[k][t]: some monotone suffix over idx[k+1..] starts
        // on the far side of candidate t of idx[k].
        let mut backward: Vec<Vec<bool>> = Vec::with_capacity(n);
        for k in 0..n {
            backward.push(vec![false; candidates[idx[k]].len()]);
        }
        for t in 0..candidates[idx[n - 1]].len() {
            backward[n - 1][t] = true;
        }
        for k in (0..n.saturating_sub(1)).rev() {
            let next = &candidates[idx[k + 1]];
            let mut beyond = vec![false; next.len() + 1];
            if increasing {
                for t in (0..next.len()).rev() {
                    beyond[t] = beyond[t + 1] || backward[k + 1][t];
                }
            } else {
                let mut under = vec![false; next.len() + 1];
                for t in 0..next.len() {
                    under[t + 1] = under[t] || backward[k + 1][t];
                }
                beyond = under;
            }
            let here = &candidates[idx[k]];
            for (t, &bp) in here.iter().enumerate() {
                if increasing {
                    let start = next.partition_point(|&p| p <= bp);
                    backward[k][t] = beyond[start];
                } else {
                    let cut = next.partition_point(|&p| p < bp);
                    backward[k][t] = beyond[cut];
                }
            }
        }
        // THE PIN: the unique forward-AND-backward-feasible candidate
        // position (DISTINCT bp — the enumeration's leaves agree on a
        // position, and duplicate candidate entries at one bp are the
        // same agreement) per anchor — the position every complete
        // chain uses; zero or >=2 distinct feasible positions leave
        // the pin None, and when no complete chain exists no candidate
        // is both-way feasible so every anchor stays None — the
        // zero-leaf outcome. (Candidates are sorted, so equal bps are
        // adjacent.)
        for k in 0..n {
            let j = idx[k];
            let mut feasible: Option<u64> = None;
            let mut feasible_count = 0u32;
            for (t, &bp) in candidates[j].iter().enumerate() {
                if forward[k][t] && backward[k][t] && feasible != Some(bp) {
                    feasible = Some(bp);
                    feasible_count += 1;
                }
            }
            if feasible_count == 1 {
                pinned[j] = feasible;
            }
        }
    }
    let ambiguous: Vec<usize> = (0..m)
        .filter(|&j| !candidates[j].is_empty() && pinned[j].is_none())
        .collect();
    (pinned, ambiguous)
}

// ---------------------------------------------------------------------------
// FASTQ streaming and the record-derivation cache (slice A, verbatim;
// the cache is a pure function of the FASTQ and the panel).
// ---------------------------------------------------------------------------

fn stream_fastq(
    path: &PathBuf,
    mut visit: impl FnMut(&[u8], &[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let (reader, _) = niffler::get_reader(Box::new(File::open(path)?)).map_err(io::Error::other)?;
    let mut lines = BufReader::new(reader).lines();
    let mut current = lines
        .next()
        .transpose()?
        .ok_or_else(|| invalid("empty reads file"))?;
    loop {
        ensure(current.starts_with('@'), "reads file is not FASTQ")?;
        let seq = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ sequence"))?;
        let plus = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ separator"))?;
        let qual = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ quality"))?;
        ensure(
            plus.starts_with('+') && seq.len() == qual.len(),
            "invalid FASTQ record",
        )?;
        visit(seq.as_bytes(), qual.as_bytes())?;
        match lines.next().transpose()? {
            Some(line) => current = line,
            None => break,
        }
    }
    Ok(())
}

#[allow(clippy::type_complexity)]
fn read_derive_cache(
    path: &std::path::Path,
) -> io::Result<(Vec<Vec<u64>>, Vec<Vec<usize>>, Vec<Vec<u8>>, Vec<u32>, Vec<Vec<(u32, Vec<(i32, u32)>)>>)> {
    let mut file = BufReader::new(File::open(path)?);
    let mut u64buf = [0u8; 8];
    let mut u32buf = [0u8; 4];
    let mut take_u64 = |file: &mut BufReader<File>| -> io::Result<u64> {
        file.read_exact(&mut u64buf)?;
        Ok(u64::from_le_bytes(u64buf))
    };
    let mut take_u32 = |file: &mut BufReader<File>| -> io::Result<u32> {
        file.read_exact(&mut u32buf)?;
        Ok(u32::from_le_bytes(u32buf))
    };
    let n_keys = take_u64(&mut file)?;
    let mut key_tokens = Vec::with_capacity(n_keys as usize);
    for _ in 0..n_keys {
        let len = take_u32(&mut file)? as usize;
        let mut tokens = Vec::with_capacity(len);
        for _ in 0..len {
            tokens.push(take_u64(&mut file)?);
        }
        key_tokens.push(tokens);
    }
    let n_lists = take_u64(&mut file)?;
    let mut key_reads = Vec::with_capacity(n_lists as usize);
    for _ in 0..n_lists {
        let len = take_u32(&mut file)? as usize;
        let mut readers = Vec::with_capacity(len);
        for _ in 0..len {
            readers.push(take_u32(&mut file)? as usize);
        }
        key_reads.push(readers);
    }
    let n_reads = take_u64(&mut file)?;
    let mut reads = Vec::with_capacity(n_reads as usize);
    for _ in 0..n_reads {
        let len = take_u32(&mut file)? as usize;
        let mut seq = vec![0u8; len];
        file.read_exact(&mut seq)?;
        reads.push(seq);
    }
    let mut read_multiplicity = Vec::with_capacity(n_reads as usize);
    for _ in 0..n_reads {
        read_multiplicity.push(take_u32(&mut file)?);
    }
    let n_record_lists = take_u64(&mut file)?;
    let mut read_records = Vec::with_capacity(n_record_lists as usize);
    for _ in 0..n_record_lists {
        let len = take_u32(&mut file)? as usize;
        let mut records = Vec::with_capacity(len);
        for _ in 0..len {
            let key = take_u32(&mut file)?;
            let walk_len = take_u32(&mut file)? as usize;
            let mut walk = Vec::with_capacity(walk_len);
            for _ in 0..walk_len {
                let mut local = [0u8; 4];
                file.read_exact(&mut local)?;
                let node = i32::from_le_bytes(local);
                let pos = take_u32(&mut file)?;
                walk.push((node, pos));
            }
            records.push((key, walk));
        }
        read_records.push(records);
    }
    Ok((key_tokens, key_reads, reads, read_multiplicity, read_records))
}

// ---------------------------------------------------------------------------
// The canonical-scheme pin skeleton (slice D's frame repair).
//
// Slice A's binding lesson, applied to the slice-B pin machinery: the
// syng's stored path walk (`walk_path_range`, what the partition GFAs'
// P lines spell) carries only ONE frame's syncmer selection, so anchor
// k-mers that are rc-frame-qualified on a path are ABSENT from that
// walk — a fold pin skeleton built from it is frame-blind and cannot
// pin records whose anchors qualify in the rc frame (the slice-C
// anatomy's named cause: 2,059 in-window full-match units at chrI L7
// with no valid truth placement). The repair reproduces the routing's
// own canonical-scheme territory step extraction (ported from
// build_territory_index_rows): per position, the frame that spells the
// k-mer's canonical form min(K, rc(K)) forward decides; the chosen
// frame's step node is kept (node identity is frame-independent — the
// same interned syncmer id, the sign carries the frame's orientation)
// and a position whose canonical frame did not qualify carries NO
// step. The census records are derived under this same scheme, so
// every anchor of every record is present at every true occurrence of
// either orientation.
// ---------------------------------------------------------------------------

/// The per-position canonical-scheme selection over the two frames'
/// step maps (positions absolute bp, step nodes signed). Frame tag:
/// 0 = the stored forward frame's selection, 1 = the rc frame's
/// qualification.
fn canonical_scheme_steps(
    forward_map: &BTreeMap<u64, i32>,
    reverse_map: &BTreeMap<u64, i32>,
    seq: &[u8],
    seq_lo: u64,
    k: u64,
) -> Vec<(u64, i32, u8)> {
    let mut positions: BTreeSet<u64> = forward_map.keys().copied().collect();
    positions.extend(reverse_map.keys().copied());
    let mut steps = Vec::with_capacity(positions.len());
    for bp in positions {
        let rel = (bp - seq_lo) as usize;
        let canonical_forward = seq
            .get(rel..rel + k as usize)
            .map(|window| {
                crate::genome_inference::mem_records::window_is_canonical_forward(window)
            })
            .unwrap_or(false);
        let chosen = if canonical_forward {
            forward_map.get(&bp)
        } else {
            reverse_map.get(&bp)
        };
        if let Some(&node) = chosen {
            let frame = if canonical_forward { 0u8 } else { 1u8 };
            steps.push((bp, node, frame));
        }
    }
    steps
}

/// The canonical-scheme steps overlapping one row range [start, end)
/// — the pin skeleton's step source. Every emitted step is verified BY
/// SEQUENCE against the AGC-fetched panel sequence (slice A's
/// exactness discipline): the window at the claimed position must
/// equal the claimed node's interned syncmer sequence, orientation-
/// aware per the sign (`syncmer_seq(negative)` is the reverse
/// complement), and the frame decision must agree with the window's
/// canonical form. Returns the steps plus the census of the stored
/// walk's diff (the diagnosis): `added` = positions the stored walk
/// lacks (the rc-frame-only anchors — the blinded class),
/// `replaced` = positions where the canonical scheme keeps a different
/// node than the stored walk, `dropped` = stored positions the
/// canonical scheme excludes (their canonical frame did not qualify —
/// no census record can anchor there).
#[allow(clippy::too_many_arguments)]
fn canonical_steps_overlapping(
    panel: &SyngIndex,
    fetch: &(dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
    path_name: &str,
    start: u64,
    end: u64,
    k: u64,
    stored: &[(u64, i32)],
) -> io::Result<(Vec<(u64, i32, u8)>, u64, u64, u64)> {
    // The rc-symmetric augmentation (the raw single-frame extraction on
    // the fetched range's reverse complement, mapped back to the
    // path's own coordinates with the frame's orientation on the
    // sign). Fetch with one syncmer length of context on each side so
    // every overlapping step's full k-mer is inside the fetched
    // sequence for the canonical-form check.
    let seq_lo = start.saturating_sub(k);
    let seq_hi = end + k;
    let seq = fetch(path_name, seq_lo, seq_hi)?;
    let rc_seq = revcomp(&seq);
    let reverse: Vec<(u64, i32)> =
        crate::genome_inference::mem_records::raw_matched_syncmers(panel, &rc_seq)?
            .into_iter()
            .map(|(signed, q)| (seq_lo + seq.len() as u64 - k - q, -signed))
            .collect();
    // Both frames' selections restricted to the stored walk's own
    // domain: steps whose windows overlap [start, end).
    let keep = |bp: u64| bp + k > start && bp < end;
    let forward_map: BTreeMap<u64, i32> = stored
        .iter()
        .copied()
        .filter(|&(bp, _)| keep(bp))
        .collect();
    let reverse_map: BTreeMap<u64, i32> =
        reverse.into_iter().filter(|&(bp, _)| keep(bp)).collect();
    let steps = canonical_scheme_steps(&forward_map, &reverse_map, &seq, seq_lo, k);
    let mut verified = 0u64;
    let mut added = 0u64;
    let mut replaced = 0u64;
    for &(bp, node, frame) in &steps {
        let rel = (bp - seq_lo) as usize;
        let window = &seq[rel..rel + k as usize];
        let mut interned = panel.syncmer_seq(node);
        // The syng's kmerHashSeq may return lowercase bases; the GFA
        // writer's own convention uppercases them (write_segment's
        // make_ascii_uppercase) — match it.
        interned.make_ascii_uppercase();
        ensure(
            window == interned.as_slice(),
            &format!(
                "canonical-scheme step fails sequence verification against the fetched \
                 sequence: {path_name} bp {bp} node {node}"
            ),
        )?;
        let canonical_forward =
            crate::genome_inference::mem_records::window_is_canonical_forward(window);
        ensure(
            canonical_forward == (frame == 0),
            "canonical-scheme frame decision disagrees with the window's canonical form",
        )?;
        verified += 1;
        match forward_map.get(&bp) {
            Some(&stored_node) if stored_node == node => {}
            Some(_) => replaced += 1,
            None => added += 1,
        }
    }
    // Completeness of the forward part: every stored step at a
    // canonical-forward position must survive with the stored node.
    for (&bp, &node) in forward_map.iter() {
        let rel = (bp - seq_lo) as usize;
        let canonical_forward = seq
            .get(rel..rel + k as usize)
            .map(|window| {
                crate::genome_inference::mem_records::window_is_canonical_forward(window)
            })
            .unwrap_or(false);
        if canonical_forward {
            ensure(
                steps.iter().any(|&(s_bp, s_node, _)| s_bp == bp && s_node == node),
                "canonical-scheme skeleton dropped a canonical-forward stored step",
            )?;
        }
    }
    let chosen_positions: BTreeSet<u64> = steps.iter().map(|&(bp, _, _)| bp).collect();
    let dropped = forward_map
        .keys()
        .filter(|&&bp| !chosen_positions.contains(&bp))
        .count() as u64;
    Ok((steps, verified + added, added + replaced, dropped))
}

// ---------------------------------------------------------------------------
// Candidate rows (slice A's RowStore, verbatim semantics; slice D: the
// pin skeleton's step source is the canonical scheme by default).
// ---------------------------------------------------------------------------

/// The per-row skeleton diff over the CONTAINED steps (windows fully
/// inside the row extent — the pin-skeleton domain): the canonical-
/// scheme steps with frame tags beside the stored walk's steps, plus
/// the position census (added = the rc-frame-only anchors the stored
/// walk lacks; replaced = positions where the canonical scheme keeps
/// a different node; dropped = stored positions the canonical scheme
/// excludes).
struct SkeletonDiff {
    canonical: Vec<(u64, i32)>,
    frames: Vec<(u64, u8)>,
    stored: Vec<(u64, i32)>,
    added: u64,
    replaced: u64,
    dropped: u64,
    kept_forward: u64,
    kept_reverse: u64,
}

fn contained_skeleton_diff(row: &RowFold, k: u64) -> SkeletonDiff {
    let contained_of = |walk: &[(u64, i32)]| {
        walk.iter()
            .filter(|&&(bp, _)| bp >= row.start && bp + k <= row.end)
            .map(|&(bp, node)| (bp - row.start, node))
            .collect::<Vec<(u64, i32)>>()
    };
    let stored = contained_of(&row.stored_walk);
    let canonical = contained_of(&row.walk);
    let frames = row
        .walk
        .iter()
        .zip(row.frames.iter())
        .filter(|(&(bp, _), _)| bp >= row.start && bp + k <= row.end)
        .map(|(&(bp, _), &frame)| (bp - row.start, frame))
        .collect::<Vec<(u64, u8)>>();
    let stored_positions: BTreeSet<u64> = stored.iter().map(|&(bp, _)| bp).collect();
    let canonical_positions: BTreeSet<u64> = canonical.iter().map(|&(bp, _)| bp).collect();
    let stored_nodes: BTreeMap<u64, i32> = stored.iter().copied().collect();
    let mut added = 0u64;
    let mut replaced = 0u64;
    let mut dropped = 0u64;
    let mut kept_forward = 0u64;
    let mut kept_reverse = 0u64;
    for &(bp, node) in &canonical {
        match stored_nodes.get(&bp) {
            None => added += 1,
            Some(&stored_node) if stored_node == node => {}
            Some(_) => replaced += 1,
        }
    }
    for &bp in &stored_positions {
        if !canonical_positions.contains(&bp) {
            dropped += 1;
        }
    }
    for &(_, frame) in &frames {
        if frame == 0 {
            kept_forward += 1;
        } else {
            kept_reverse += 1;
        }
    }
    SkeletonDiff {
        canonical,
        frames,
        stored,
        added,
        replaced,
        dropped,
        kept_forward,
        kept_reverse,
    }
}

#[derive(Clone)]
struct RowFold {
    #[allow(dead_code)]
    path: usize,
    start: u64,
    end: u64,
    seq: Vec<u8>,
    /// The pin skeleton's steps overlapping the row range (absolute
    /// bp, signed node), sorted by bp — the CANONICAL-SCHEME selection
    /// by default (slice D's frame repair), or the stored path walk
    /// under the identity gate (the slice-B/C before-record).
    walk: Vec<(u64, i32)>,
    /// The frame tag per walk step (0 = the stored forward frame's
    /// selection, 1 = the rc frame's qualification) — the audit view;
    /// empty under the identity gate.
    frames: Vec<u8>,
    /// The stored path walk (the committed partition GFAs' own P-line
    /// steps) — the audit's comparison baseline.
    stored_walk: Vec<(u64, i32)>,
    /// (partition, member index) pairs of every partition holding this
    /// row.
    members: Vec<(u32, usize)>,
}

/// The axis-partition row store: rows loaded lazily per partition,
/// deduplicated by (path, start, end) across partitions.
struct RowStore<'a> {
    panel: &'a SyngIndex,
    path_of_name: &'a HashMap<String, usize>,
    fetch: &'a (dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
    k: u64,
    /// The identity gate: keep the stored path walk as the pin
    /// skeleton (env IMPG_REALIGN_STORED_WALK_SKELETON; assessment-side
    /// diagnostic for the before/after pairing — the repair is the
    /// default).
    stored_frame: bool,
    folds: Vec<RowFold>,
    dedup: HashMap<(usize, u64, u64), usize>,
    by_partition: HashMap<u32, Vec<usize>>,
    by_partition_path: HashMap<(u32, usize), Vec<usize>>,
}

impl<'a> RowStore<'a> {
    fn new(
        panel: &'a SyngIndex,
        path_of_name: &'a HashMap<String, usize>,
        fetch: &'a (dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
        k: u64,
        stored_frame: bool,
    ) -> Self {
        RowStore {
            panel,
            path_of_name,
            fetch,
            k,
            stored_frame,
            folds: Vec::new(),
            dedup: HashMap::new(),
            by_partition: HashMap::new(),
            by_partition_path: HashMap::new(),
        }
    }

    fn ensure_partition(&mut self, members: &[MemberRow], partition: u32) -> io::Result<()> {
        if self.by_partition.contains_key(&partition) {
            return Ok(());
        }
        let mut indices = Vec::with_capacity(members.len());
        for (member_index, member) in members.iter().enumerate() {
            let path_idx = *self
                .path_of_name
                .get(&member.path_name)
                .ok_or_else(|| invalid("member path absent from the syng"))?;
            let end = member.end.min(self.panel.name_map.path_to_length[path_idx]);
            ensure(member.start < end, "degenerate member row")?;
            let dedup_key = (path_idx, member.start, end);
            let fold_index = match self.dedup.get(&dedup_key) {
                Some(&index) => index,
                None => {
                    let seq = (self.fetch)(&member.path_name, member.start, end)?;
                    let mut stored_walk: Vec<(u64, i32)> = self
                        .panel
                        .walk_path_range(path_idx, member.start, end)?
                        .into_iter()
                        .map(|(node, bp)| (bp, node))
                        .collect();
                    stored_walk.sort_unstable_by_key(|&(bp, _)| bp);
                    let (walk, frames) = if self.stored_frame {
                        (stored_walk.clone(), Vec::new())
                    } else {
                        let (steps, _verified, _changed, _dropped) = canonical_steps_overlapping(
                            self.panel,
                            self.fetch,
                            &member.path_name,
                            member.start,
                            end,
                            self.k,
                            &stored_walk,
                        )?;
                        (
                            steps.iter().map(|&(bp, node, _)| (bp, node)).collect(),
                            steps.iter().map(|&(_, _, frame)| frame).collect(),
                        )
                    };
                    let index = self.folds.len();
                    self.folds.push(RowFold {
                        path: path_idx,
                        start: member.start,
                        end,
                        seq,
                        walk,
                        frames,
                        stored_walk,
                        members: Vec::new(),
                    });
                    self.dedup.insert(dedup_key, index);
                    index
                }
            };
            self.folds[fold_index].members.push((partition, member_index));
            self.by_partition_path
                .entry((partition, path_idx))
                .or_default()
                .push(fold_index);
            indices.push(fold_index);
        }
        self.by_partition.insert(partition, indices);
        Ok(())
    }

    fn rows_of_partition_path(&self, partition: u32, path: usize) -> &[usize] {
        self.by_partition_path
            .get(&(partition, path))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
}

// ---------------------------------------------------------------------------
// The cluster machinery (mirrored pure functions from
// panel_route_spine/cosine_probe.rs — the existing instrument of
// record's own functions; verbatim semantics, no threshold altered).
// ---------------------------------------------------------------------------

type GraphNode = u32;
type PackedEdge = u64;

// ---------------------------------------------------------------------------
// (THE PANEL-SIDE PER-LOCUS PRECOMPUTE, SS10.6.23 slice 1 — the
// ONE-SPELLING helpers. The fold's derived content, extracted
// VERBATIM from the live fold build so the live path, the precompute
// builder, and the serve all derive a fold's walk / node positions /
// multisets through THIS code — the served structure and the live-built
// one are the same derivation by construction, never two spellings
// that could drift.)
// ---------------------------------------------------------------------------

/// The row's CONTAINED walk steps (windows fully inside the row
/// extent), relative to the row start — the live fold build's
/// `relative` (the grouping key's walk half, and the fold's own
/// pinnable skeleton).
fn contained_relative_walk(row: &RowFold, k: u64) -> Vec<(u64, i32)> {
    row.walk
        .iter()
        .filter(|&&(bp, _)| bp >= row.start && bp + k <= row.end)
        .map(|&(bp, node)| (bp - row.start, node))
        .collect()
}

/// The fold's node-position table (node abs id -> the relative
/// (bp, sign) occurrences of the contained steps, in walk order).
fn fold_node_positions(row: &RowFold, k: u64) -> HashMap<u32, Vec<(u64, i32)>> {
    let mut node_positions: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
    for &(bp, node) in &row.walk {
        if bp >= row.start && bp + k <= row.end {
            node_positions
                .entry(node.unsigned_abs())
                .or_default()
                .push((bp - row.start, node.signum()));
        }
    }
    node_positions
}

/// The fold's node-usage and edge-usage multisets over the contained
/// relative walk (the QUAL spectrum's material signature), sorted by
/// key with the counts coalesced.
fn fold_multisets(contained: &[(u64, i32)]) -> (Vec<(u64, u32)>, Vec<(u64, u32)>) {
    let mut nodes: Vec<(u64, u32)> = contained
        .iter()
        .map(|&(_, node)| (node.unsigned_abs() as u64, 1))
        .collect();
    nodes.sort_unstable();
    nodes.dedup_by(|a, b| {
        if a.0 == b.0 {
            b.1 += a.1;
            true
        } else {
            false
        }
    });
    let mut edges: Vec<(u64, u32)> = contained
        .windows(2)
        .map(|w| (pack_edge(w[0].1.unsigned_abs(), w[1].1.unsigned_abs()), 1))
        .collect();
    edges.sort_unstable();
    edges.dedup_by(|a, b| {
        if a.0 == b.0 {
            b.1 += a.1;
            true
        } else {
            false
        }
    });
    (nodes, edges)
}

/// THE LIVE FOLD LAYOUT BUILD (the per-window fold grouping, the
/// slice-1 pure function extracted from `compute_locus` VERBATIM — the
/// precompute builder runs THIS function unmasked to record the
/// grouping; the live call runs it with its masking filters; the serve
/// reproduces its output from the recorded grouping; the unit tests
/// prove serve == live on randomized fixtures). Returns the folds in
/// first-appearance order, the representative ROW id per fold, the
/// fold id per member (the builder's recorded grouping), and the
/// filter counts.
#[allow(clippy::too_many_arguments)]
fn build_window_folds_live(
    rows: &[RowFold],
    dedup: &HashMap<(usize, u64, u64), usize>,
    members: &[MemberRow],
    panel_component: &str,
    axis_start: u64,
    axis_end: u64,
    flank_overlap: Option<FlankOverlap>,
    exclude_paths: &[String],
    path_of_name: &HashMap<String, usize>,
    path_len: &dyn Fn(usize) -> u64,
    k: u64,
) -> io::Result<(Vec<Fold>, Vec<usize>, Vec<usize>, u64, u64)> {
    let mut fold_map: BTreeMap<(Vec<u8>, Vec<(u64, i32)>), usize> = BTreeMap::new();
    let mut folds: Vec<Fold> = Vec::new();
    let mut fold_rep_rows: Vec<usize> = Vec::new();
    let mut member_fold: Vec<usize> = Vec::with_capacity(members.len());
    let mut excluded_rows: u64 = 0;
    let mut flank_excluded_rows: u64 = 0;
    for member in members {
        // THE FLANK-OVERLAP CANDIDATE FILTER (the owner's design):
        // a member row that does not cover BOTH FLANKS of the
        // window's axis extent - the first N bp and the last N bp,
        // in this frame (the partition maps' own coordinates) -
        // leaves the candidate universe here, at fold construction,
        // exactly as absent rows would. THE AXIS ROW IS NEVER
        // EXCLUDED (it is the coordinate frame; the --exclude-paths
        // precedent). Loudly counted per locus.
        if let Some(flank) = flank_overlap {
            if member.path_name != panel_component
                && !row_covers_flanks(member.start, member.end, axis_start, axis_end, flank)
            {
                flank_excluded_rows += 1;
                member_fold.push(usize::MAX);
                continue;
            }
        }
        // THE MASKING OPERATOR (the call-locus --exclude-paths
        // surface): the stranger filter folded onto the call path
        // - matching rows leave the candidate universe here, at
        // fold construction, exactly as absent rows would (the
        // axis row is never excluded: it is the frame). Loudly
        // counted per locus.
        if row_is_excluded(&member.path_name, panel_component, exclude_paths) {
            excluded_rows += 1;
            member_fold.push(usize::MAX);
            continue;
        }
        let path_idx = *path_of_name
            .get(&member.path_name)
            .ok_or_else(|| invalid("member path absent from the syng"))?;
        let row_index = *dedup
            .get(&(path_idx, member.start, member.end.min(path_len(path_idx))))
            .ok_or_else(|| invalid("member row missing from the store"))?;
        let row = &rows[row_index];
        // The CONTAINED steps only (windows fully inside the row
        // extent): an anchor whose window overlaps the row edge
        // cannot be a pinned backbone anchor (the read must lie
        // fully inside), and keeping it would place its window
        // start at a NEGATIVE relative position - the u64
        // subtraction would wrap and poison the monotone
        // enumeration. Edge-overlapping anchors abstain (the named
        // remainder class).
        let relative = contained_relative_walk(row, k);
        let key = (row.seq.clone(), relative.clone());
        let index = *fold_map.entry(key).or_insert_with(|| {
            let node_positions = fold_node_positions(row, k);
            let (nodes, edges) = fold_multisets(&relative);
            folds.push(Fold {
                seq: row.seq.clone(),
                len: (row.end - row.start) as u64,
                walk: relative.clone(),
                node_positions,
                members: Vec::new(),
                nodes,
                edges,
            });
            fold_rep_rows.push(row_index);
            folds.len() - 1
        });
        folds[index].members.push(member.clone());
        member_fold.push(index);
    }
    Ok((
        folds,
        fold_rep_rows,
        member_fold,
        excluded_rows,
        flank_excluded_rows,
    ))
}

// ---------------------------------------------------------------------------
// (THE BIND SCAFFOLDING, SS10.6.25 slice 2 — the canonical-scheme step
// derivation over one panel range, ONE SPELLING: the bind's former
// `canonical_steps_near` closure extracted VERBATIM with the window
// bounds as parameters, so the live center call, the precompute
// builder, and the serve's uncovered-position fallback all run THE
// SAME code — no second spelling of the derivation exists anywhere.
// The step value at a position P is a function of the panel's stored
// walk entry at P, the path bytes in [P-(w-1), P+(w-1)+syncmer_len)
// clamped to the path (the reverse-strand syncmer window's own k-mer
// span; the canonical orientation check is local), and the panel
// dictionary — the LOCALITY LEMMA the serve's served band is derived
// from, every constant the syncmer window's own parameter, never
// tuned.)
// ---------------------------------------------------------------------------

/// The canonical-scheme territory steps over one panel range
/// `[lo, hi)` — `canonical_steps_near(path, around)` verbatim for
/// `lo = around-150`, `hi = around+150` (the live bind's own window;
/// the extraction is code motion and the answer gates are its
/// proof). The forward map is the panel's own stored walk slice; the
/// reverse map is the reverse-complement view's matched syncmers,
/// mapped back to forward positions; each position chooses its frame
/// by the canonical orientation rule.
fn canonical_steps_over_range(
    panel: &SyngIndex,
    fetch_seq: &dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
    k: u64,
    path: usize,
    lo: u64,
    hi: u64,
) -> io::Result<Vec<(u64, i32)>> {
    let forward: Vec<(u64, i32)> = panel
        .walk_path_range(path, lo, hi)?
        .into_iter()
        .map(|(node, bp)| (bp, node))
        .collect();
    let seq_lo = lo.saturating_sub(k);
    let seq_hi = hi + k;
    let seq = fetch_seq(&panel.name_map.path_to_name[path], seq_lo, seq_hi)?;
    let rc_seq = revcomp(&seq);
    let reverse: Vec<(u64, i32)> =
        crate::genome_inference::mem_records::raw_matched_syncmers(panel, &rc_seq)?
            .into_iter()
            .map(|(signed, q)| {
                (
                    seq_lo + seq.len() as u64 - k - q,
                    -signed,
                )
            })
            .collect();
    let forward_map: BTreeMap<u64, i32> = forward.into_iter().collect();
    let reverse_map: BTreeMap<u64, i32> = reverse.into_iter().collect();
    let mut positions: BTreeSet<u64> = forward_map.keys().copied().collect();
    positions.extend(reverse_map.keys().copied());
    let mut steps = Vec::with_capacity(positions.len());
    for bp in positions {
        let rel = (bp - seq_lo) as usize;
        let canonical_forward = seq
            .get(rel..rel + k as usize)
            .map(|window| {
                crate::genome_inference::mem_records::window_is_canonical_forward(window)
            })
            .unwrap_or(false);
        let chosen = if canonical_forward {
            forward_map.get(&bp)
        } else {
            reverse_map.get(&bp)
        };
        if let Some(&node) = chosen {
            steps.push((bp, node));
        }
    }
    Ok(steps)
}

/// THE SERVED BAND of one derivation range `[lo, hi)` (the locality
/// lemma's own predicate, every constant derived from the syncmer
/// window's parameters: `reach` = w-1, `syncmer_len` = w+smer): the
/// position's decision context — the path bytes
/// `[P-reach, P+reach+syncmer_len)` clamped to the path — must lie
/// inside the range's fetch `[lo-syncmer_len, hi+syncmer_len)` (the
/// walk-slice clause `P in [lo,hi)` is checked by the caller), and at
/// the path's end the clamped context is served when the range's
/// fetch reaches the path end (the clamped-context clause — the last
/// syncmer window of a path truncates identically in both
/// derivations).
fn bind_range_covers(
    lo: u64,
    hi: u64,
    position: u64,
    path_len: u64,
    reach: u64,
    syncmer_len: u64,
) -> bool {
    position >= lo
        && position < hi
        && (position + reach <= hi
            || (position + reach + syncmer_len >= path_len
                && path_len <= hi + syncmer_len))
}

/// THE BIND SCAFFOLDING BUILDER (SS10.6.25 slice 2 — the once-per-panel
/// price): per panel path, the member-row extents (the maps' own rows —
/// the panel-knowable coordinate set; the census's positions are
/// sample-side, A4's own statement) merged when they sit within the
/// syncmer window's own reach, each merged extent `[s,e)` derived over
/// `[s, e+reach)` so the served band covers the extent exactly —
/// through the SAME `canonical_steps_over_range` the live bind runs
/// (the builder IS the live derivation, the one-spelling guarantee).
#[allow(clippy::too_many_arguments)]
fn build_bind_scaffold(
    panel: &SyngIndex,
    fetch_seq: &dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
    k: u64,
    member_rows: impl Iterator<Item = (usize, u64, u64)>,
) -> io::Result<Vec<BindScaffoldPathAsset>> {
    let reach = panel.params.w as u64 - 1;
    let mut by_path: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    for (path, start, end) in member_rows {
        let path_len = panel.name_map.path_to_length[path];
        by_path
            .entry(path)
            .or_default()
            .push((start, end.min(path_len)));
    }
    let mut paths: Vec<BindScaffoldPathAsset> = Vec::with_capacity(by_path.len());
    for (path, mut extents) in by_path {
        extents.sort_unstable();
        // Merge the extents the window's own reach apart (contiguous
        // coverage without gaps), then derive over `[s, e+reach)` so
        // the served band `[s, e+reach)` with `P+reach <= hi` covers
        // the extent exactly.
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (start, end) in extents {
            match merged.last_mut() {
                Some(last) if start <= last.1 + reach => last.1 = last.1.max(end),
                _ => merged.push((start, end)),
            }
        }
        let ranges: Vec<(u64, u64)> = merged
            .iter()
            .map(|&(start, end)| (start, end + reach))
            .collect();
        let mut steps: Vec<(u64, i32)> = Vec::new();
        for &(lo, hi) in &ranges {
            steps.extend(canonical_steps_over_range(
                panel, fetch_seq, k, path, lo, hi,
            )?);
        }
        steps.sort_unstable_by_key(|&(bp, _)| bp);
        paths.push(BindScaffoldPathAsset {
            path_name: panel.name_map.path_to_name[path].clone(),
            ranges,
            steps,
        });
    }
    Ok(paths)
}
// ---------------------------------------------------------------------------
// (THE PANEL-SIDE PER-LOCUS PRECOMPUTE ASSET — SS10.6.23 slice 1's
// panel-side artifact: the row context (the RowStore rows, the
// sequences replaced by their FNV-64 CONTENT KEY) and the per-window
// fold groupings (representative row ids + per-member fold ids),
// serialized exact-width little-endian with NO FLOATS anywhere, so
// the serve reconstructs bit-identical structures by construction.
// The asset's manifest keys on CONTENT, never filenames: the panel
// identity object, the FNV digest over the sorted *.gfa.map.json
// bytes, the routes graph.json digest, the component, the frame +
// record origin, the stored-frame convention, k.)
// ---------------------------------------------------------------------------

/// The FNV-1a over the partition-graph maps' CONTENT, in sorted file
/// name order — the manifest's own drift detector (never a filename
/// alone; the vendor-pointer-drift class is the cautionary receipt).
fn partition_maps_digest(dir: &Path) -> io::Result<u64> {
    let mut names: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        if name.ends_with(".gfa.map.json") {
            names.push(name);
        }
    }
    names.sort();
    let mut digest = fnv1a64(b"");
    for name in &names {
        let bytes = std::fs::read(dir.join(name))?;
        hash_bytes(&mut digest, name.as_bytes());
        hash_bytes(&mut digest, &[0xFF]);
        hash_bytes(&mut digest, &bytes);
    }
    Ok(digest)
}

fn routes_graph_digest(path: &Path) -> io::Result<u64> {
    let bytes = std::fs::read(path)?;
    Ok(fnv1a64(&bytes))
}

fn hash_bytes(hash: &mut u64, bytes: &[u8]) {
    for &byte in bytes {
        *hash ^= byte as u64;
        *hash = hash.wrapping_mul(0x100_0000_01b3);
    }
}

#[derive(Serialize, Deserialize)]
pub struct PanelPrecomputeManifest {
    format: String,
    panel: PanelIdentity,
    k: u64,
    component: String,
    frame: String,
    record_origin: u64,
    stored_frame: bool,
    maps_digest: u64,
    routes_digest: u64,
    windows: u64,
    rows: u64,
    folds: u64,
    // (SS10.6.25 slice 2 — the bind-scaffolding section's counters,
    // absent in the slice-1 /1 manifests and defaulted there: the
    // section is validated against them at load.)
    #[serde(default)]
    bind_paths: u64,
    #[serde(default)]
    bind_ranges: u64,
    #[serde(default)]
    bind_steps: u64,
}

impl PanelPrecomputeManifest {
    /// The asset's format string — the reuse check's own key (the
    /// panel-prep step upgrades a slice-1 /1 asset loudly rather than
    /// reusing an asset without the bind-scaffolding section).
    pub fn format_text(&self) -> String {
        self.format.clone()
    }
}

/// One path's bind scaffolding (SS10.6.25 slice 2): the derivation
/// ranges the builder ran (the served band is derived from them), and
/// the canonical-scheme steps the derivation returned, sorted by bp —
/// exact-width integers only, NO FLOATS anywhere.
#[derive(Clone)]
struct BindScaffoldPathAsset {
    path_name: String,
    ranges: Vec<(u64, u64)>,
    steps: Vec<(u64, i32)>,
}

struct PrecomputeRowAsset {
    path_name: String,
    start: u64,
    end: u64,
    seq_fnv: u64,
    stored_walk: Vec<(u64, i32)>,
    walk: Vec<(u64, i32)>,
    frames: Vec<u8>,
    members: Vec<(u32, usize)>,
}

struct PrecomputeWindowAsset {
    partition: u32,
    axis_start: u64,
    axis_end: u64,
    fold_rep: Vec<usize>,
    member_fold: Vec<usize>,
}

/// The loaded, validated asset beside the call's own partition maps.
struct PanelPrecomputeAsset {
    manifest: PanelPrecomputeManifest,
    rows: Vec<PrecomputeRowAsset>,
    windows: Vec<PrecomputeWindowAsset>,
    // (SS10.6.25 slice 2 — None for the slice-1 /1 format: the asset
    // predates the bind-scaffolding section and the bind derives
    // live, loudly stated.)
    bind: Option<Vec<BindScaffoldPathAsset>>,
}

// (The exact-width little-endian encoding: length-prefixed vectors,
// u64 bp / i32 signed node / u8 frame; NO FLOATS anywhere in the
// asset, so the round trip is lossless by construction.)
struct AssetWriter {
    buf: Vec<u8>,
}
impl AssetWriter {
    fn u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    fn i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    fn u8(&mut self, v: u8) {
        self.buf.push(v);
    }
    fn str(&mut self, v: &str) {
        self.u64(v.len() as u64);
        self.buf.extend_from_slice(v.as_bytes());
    }
    fn walk(&mut self, v: &[(u64, i32)]) {
        self.u64(v.len() as u64);
        for &(bp, node) in v {
            self.u64(bp);
            self.i32(node);
        }
    }
    fn bytes(&mut self, v: &[u8]) {
        self.u64(v.len() as u64);
        self.buf.extend_from_slice(v);
    }
    fn row_members(&mut self, v: &[(u32, usize)]) {
        self.u64(v.len() as u64);
        for &(partition, member) in v {
            self.u64(partition as u64);
            self.u64(member as u64);
        }
    }
    fn usize_vec(&mut self, v: &[usize]) {
        self.u64(v.len() as u64);
        for &value in v {
            self.u64(value as u64);
        }
    }
}

struct AssetReader<'a> {
    buf: &'a [u8],
    pos: usize,
}
impl<'a> AssetReader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        AssetReader { buf, pos: 0 }
    }
    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        ensure(
            self.pos + n <= self.buf.len(),
            "the panel-precompute asset is truncated",
        )?;
        let slice = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
    }
    fn u64(&mut self) -> io::Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u8(&mut self) -> io::Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn str(&mut self) -> io::Result<String> {
        let len = self.u64()? as usize;
        Ok(String::from_utf8(self.take(len)?.to_vec())
            .map_err(|_| invalid("the panel-precompute asset carries an invalid path name"))?)
    }
    fn walk(&mut self) -> io::Result<Vec<(u64, i32)>> {
        let len = self.u64()? as usize;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push((self.u64()?, self.i32()?));
        }
        Ok(out)
    }
    fn bytes(&mut self) -> io::Result<Vec<u8>> {
        let len = self.u64()? as usize;
        Ok(self.take(len)?.to_vec())
    }
    fn row_members(&mut self) -> io::Result<Vec<(u32, usize)>> {
        let len = self.u64()? as usize;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            let partition = self.u64()? as u32;
            let member = self.u64()? as usize;
            out.push((partition, member));
        }
        Ok(out)
    }
    fn usize_vec(&mut self) -> io::Result<Vec<usize>> {
        let len = self.u64()? as usize;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(self.u64()? as usize);
        }
        Ok(out)
    }
    fn done(&self) -> io::Result<()> {
        ensure(
            self.pos == self.buf.len(),
            "the panel-precompute asset carries trailing bytes",
        )
    }
}

fn encode_rows(rows: &[PrecomputeRowAsset]) -> Vec<u8> {
    let mut w = AssetWriter { buf: Vec::new() };
    w.u64(rows.len() as u64);
    for row in rows {
        w.str(&row.path_name);
        w.u64(row.start);
        w.u64(row.end);
        w.u64(row.seq_fnv);
        w.walk(&row.stored_walk);
        w.walk(&row.walk);
        w.bytes(&row.frames);
        w.row_members(&row.members);
    }
    w.buf
}

fn decode_rows(buf: &[u8]) -> io::Result<Vec<PrecomputeRowAsset>> {
    let mut r = AssetReader::new(buf);
    let len = r.u64()? as usize;
    let mut rows = Vec::with_capacity(len);
    for _ in 0..len {
        rows.push(PrecomputeRowAsset {
            path_name: r.str()?,
            start: r.u64()?,
            end: r.u64()?,
            seq_fnv: r.u64()?,
            stored_walk: r.walk()?,
            walk: r.walk()?,
            frames: r.bytes()?,
            members: r.row_members()?,
        });
    }
    r.done()?;
    Ok(rows)
}

fn encode_windows(windows: &[PrecomputeWindowAsset]) -> Vec<u8> {
    let mut w = AssetWriter { buf: Vec::new() };
    w.u64(windows.len() as u64);
    for window in windows {
        w.u64(window.partition as u64);
        w.u64(window.axis_start);
        w.u64(window.axis_end);
        w.usize_vec(&window.fold_rep);
        w.usize_vec(&window.member_fold);
    }
    w.buf
}

fn decode_windows(buf: &[u8]) -> io::Result<Vec<PrecomputeWindowAsset>> {
    let mut r = AssetReader::new(buf);
    let len = r.u64()? as usize;
    let mut windows = Vec::with_capacity(len);
    for _ in 0..len {
        windows.push(PrecomputeWindowAsset {
            partition: r.u64()? as u32,
            axis_start: r.u64()?,
            axis_end: r.u64()?,
            fold_rep: r.usize_vec()?,
            member_fold: r.usize_vec()?,
        });
    }
    r.done()?;
    Ok(windows)
}

// (SS10.6.25 slice 2 — the bind-scaffolding section's exact-width
// little-endian encoding: length-prefixed vectors, u64 bp / i32
// signed node / u64 range bounds, NO FLOATS anywhere in the asset,
// so the round trip is lossless by construction.)
fn encode_bind_scaffolding(paths: &[BindScaffoldPathAsset]) -> Vec<u8> {
    let mut w = AssetWriter { buf: Vec::new() };
    w.u64(paths.len() as u64);
    for path in paths {
        w.str(&path.path_name);
        w.u64(path.ranges.len() as u64);
        for &(lo, hi) in &path.ranges {
            w.u64(lo);
            w.u64(hi);
        }
        w.u64(path.steps.len() as u64);
        for &(bp, node) in &path.steps {
            w.u64(bp);
            w.i32(node);
        }
    }
    w.buf
}

fn decode_bind_scaffolding(buf: &[u8]) -> io::Result<Vec<BindScaffoldPathAsset>> {
    let mut r = AssetReader::new(buf);
    let len = r.u64()? as usize;
    let mut paths = Vec::with_capacity(len);
    for _ in 0..len {
        let path_name = r.str()?;
        let ranges_len = r.u64()? as usize;
        let mut ranges = Vec::with_capacity(ranges_len);
        for _ in 0..ranges_len {
            let lo = r.u64()?;
            let hi = r.u64()?;
            ranges.push((lo, hi));
        }
        let steps_len = r.u64()? as usize;
        let mut steps = Vec::with_capacity(steps_len);
        for _ in 0..steps_len {
            let bp = r.u64()?;
            let node = r.i32()?;
            steps.push((bp, node));
        }
        paths.push(BindScaffoldPathAsset {
            path_name,
            ranges,
            steps,
        });
    }
    r.done()?;
    Ok(paths)
}

/// The expected manifest keys, computed by the caller from its OWN
/// live inputs — the serve compares field by field and names the
/// mismatching key in its error (never a silent divergence, never a
/// stale serve).
struct PanelPrecomputeKeys {
    panel: PanelIdentity,
    k: u64,
    component: String,
    frame: String,
    record_origin: u64,
    stored_frame: bool,
    maps_digest: u64,
    routes_digest: u64,
}

fn load_panel_precompute(
    dir: &Path,
    keys: &PanelPrecomputeKeys,
) -> io::Result<PanelPrecomputeAsset> {
    let manifest: PanelPrecomputeManifest = read_json(&dir.join("manifest.json"))?;
    // (SS10.6.25 slice 2: the slice-1 /1 format still serves the row
    // context + fold layouts under the default toggle — its bind
    // derives live, loudly stated; the /2 format adds the
    // bind-scaffolding section. Any other format string is unknown.)
    ensure(
        manifest.format == "impg-panel-precompute-loci/1"
            || manifest.format == "impg-panel-precompute-loci/2",
        &format!(
            "the panel-precompute asset {} carries an unknown format {:?}",
            dir.display(),
            manifest.format
        ),
    )?;
    let named = |field: &str, got: String, want: String| {
        format!(
            "the panel-precompute asset {} mismatches the live inputs at {field}: \
             asset {got:?} vs call {want:?}",
            dir.display()
        )
    };
    ensure(
        manifest.panel == keys.panel,
        &named(
            "the panel identity",
            format!("{:?}", manifest.panel),
            format!("{:?}", keys.panel),
        ),
    )?;
    ensure(
        manifest.k == keys.k,
        &named("k", manifest.k.to_string(), keys.k.to_string()),
    )?;
    ensure(
        manifest.component == keys.component,
        &named("the component", manifest.component.clone(), keys.component.clone()),
    )?;
    ensure(
        manifest.frame == keys.frame,
        &named("the frame", manifest.frame.clone(), keys.frame.clone()),
    )?;
    ensure(
        manifest.record_origin == keys.record_origin,
        &named(
            "the record origin",
            manifest.record_origin.to_string(),
            keys.record_origin.to_string(),
        ),
    )?;
    ensure(
        manifest.stored_frame == keys.stored_frame,
        &named(
            "the stored-frame convention",
            manifest.stored_frame.to_string(),
            keys.stored_frame.to_string(),
        ),
    )?;
    ensure(
        manifest.maps_digest == keys.maps_digest,
        &named(
            "the partition maps digest",
            format!("{:016x}", manifest.maps_digest),
            format!("{:016x}", keys.maps_digest),
        ),
    )?;
    ensure(
        manifest.routes_digest == keys.routes_digest,
        &named(
            "the routes graph digest",
            format!("{:016x}", manifest.routes_digest),
            format!("{:016x}", keys.routes_digest),
        ),
    )?;
    let rows = decode_rows(&std::fs::read(dir.join("rows.bin"))?)?;
    let windows = decode_windows(&std::fs::read(dir.join("windows.bin"))?)?;
    ensure(
        rows.len() as u64 == manifest.rows,
        "the panel-precompute asset's row count disagrees with its manifest",
    )?;
    ensure(
        windows.len() as u64 == manifest.windows,
        "the panel-precompute asset's window count disagrees with its manifest",
    )?;
    let bind = if manifest.format == "impg-panel-precompute-loci/2" {
        let bind = decode_bind_scaffolding(&std::fs::read(dir.join("bind-scaffolding.bin"))?)?;
        let ranges: u64 = bind.iter().map(|path| path.ranges.len() as u64).sum();
        let steps: u64 = bind.iter().map(|path| path.steps.len() as u64).sum();
        ensure(
            bind.len() as u64 == manifest.bind_paths
                && ranges == manifest.bind_ranges
                && steps == manifest.bind_steps,
            "the panel-precompute asset's bind-scaffolding counters disagree \
             with the section",
        )?;
        Some(bind)
    } else {
        None
    };
    Ok(PanelPrecomputeAsset {
        manifest,
        rows,
        windows,
        bind,
    })
}

/// The serve toggle's own modes (the labeled-toggle doctrine's
/// shape — see the serve block in `run_injected`).
#[derive(Clone, Copy, PartialEq)]
enum PanelPrecomputeMode {
    Auto,
    Off,
    On,
}

/// The serve-side state: the validated window groupings, the rows'
/// sequence content keys, the per-row sequence cache (the fetched
/// representative-row bytes, FNV-verified against the asset's key on
/// first touch — the drift detector that never hashes the whole AGC),
/// and the bind scaffolding (SS10.6.25 slice 2 — None when the asset
/// predates the section; the bind then derives live, loudly).
struct PrecomputeServe {
    windows: Vec<PrecomputeWindowAsset>,
    seq_fnv: Vec<u64>,
    seqs: Mutex<HashMap<usize, Arc<Vec<u8>>>>,
    bind: Option<BindScaffoldServe>,
}

/// THE BIND SCAFFOLDING SERVE (SS10.6.25 slice 2): per panel path, the
/// builder's derivation ranges and the recorded canonical-scheme
/// steps (sorted by bp), plus the serve's own loud counters — how many
/// record lookups the asset answered and how many positions fell
/// back to the live derivation (the census's own positions outside
/// the panel-knowable coverage; exact by construction — the fallback
/// IS the committed path). Read-only shared state under the bind's
/// parallel seam: Sync by construction (immutable maps + atomics).
struct BindScaffoldServe {
    by_path: HashMap<usize, BindScaffoldPathServed>,
    reach: u64,
    syncmer_len: u64,
    served: AtomicU64,
    fallbacks: AtomicU64,
}

struct BindScaffoldPathServed {
    ranges: Vec<(u64, u64)>,
    steps: Vec<(u64, i32)>,
}

impl BindScaffoldServe {
    /// Build the serve state from the asset's section: every path name
    /// resolves against the call's own panel (a name that does not is a
    /// content-key mismatch, loud); the band constants come from the
    /// panel's own syncmer parameters (never tuned).
    fn from_asset(
        paths: Vec<BindScaffoldPathAsset>,
        path_of_name: &HashMap<String, usize>,
        panel: &SyngIndex,
    ) -> io::Result<Self> {
        let mut by_path = HashMap::with_capacity(paths.len());
        for path in paths {
            let path_idx = *path_of_name
                .get(&path.path_name)
                .ok_or_else(|| {
                    invalid(&format!(
                        "the panel-precompute asset's bind scaffolding names a path \
                         absent from the syng: {:?}",
                        path.path_name
                    ))
                })?;
            let previous = by_path.insert(
                path_idx,
                BindScaffoldPathServed {
                    ranges: path.ranges,
                    steps: path.steps,
                },
            );
            ensure(
                previous.is_none(),
                "the panel-precompute asset's bind scaffolding names a path twice",
            )?;
        }
        Ok(BindScaffoldServe {
            by_path,
            reach: panel.params.w as u64 - 1,
            syncmer_len: panel.syncmer_length_bp() as u64,
            served: AtomicU64::new(0),
            fallbacks: AtomicU64::new(0),
        })
    }

    /// The canonical step at (path, position) — `Some(Some(step))`
    /// when the position lies inside a derivation range's SERVED BAND
    /// (the recorded step at the position, or its absence — the bind's
    /// own no-step error fires exactly as the live derivation's
    /// would); `None` when the position is UNCOVERED (the caller
    /// derives live, verbatim). The locality lemma is the exactness
    /// argument: the builder ran the same derivation over a range
    /// whose fetch contains the position's own decision context, so
    /// the recorded value is the live center derivation's own bits.
    fn step_at(&self, path: usize, position: u64, path_len: u64) -> Option<Option<(u64, i32)>> {
        let served = self
            .by_path
            .get(&path)
            .filter(|served| {
                served.ranges.iter().any(|&(lo, hi)| {
                    bind_range_covers(
                        lo,
                        hi,
                        position,
                        path_len,
                        self.reach,
                        self.syncmer_len,
                    )
                })
            })
            .map(|served| {
                served
                    .steps
                    .binary_search_by_key(&position, |&(bp, _)| bp)
                    .ok()
                    .map(|index| served.steps[index])
            });
        match served {
            Some(step) => {
                self.served.fetch_add(1, Ordering::Relaxed);
                Some(step)
            }
            None => {
                self.fallbacks.fetch_add(1, Ordering::Relaxed);
                None
            }
        }
    }
}

impl PrecomputeServe {
    /// The representative row's sequence, fetched LIVE through the
    /// same validated fetch path (the sequences stay in the sources —
    /// the design's own fetch-by-coordinate answer) and verified
    /// against the asset's recorded FNV-64 content key.
    #[allow(clippy::too_many_arguments)]
    fn seq_of(
        &self,
        row_index: usize,
        rows: &[RowFold],
        fetch: &dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
        path_names: &[String],
    ) -> io::Result<Arc<Vec<u8>>> {
        if let Some(seq) = self.seqs.lock().unwrap().get(&row_index) {
            return Ok(seq.clone());
        }
        let row = &rows[row_index];
        let seq = fetch(&path_names[row.path], row.start, row.end)?;
        ensure(
            fnv1a64(&seq) == self.seq_fnv[row_index],
            &format!(
                "the panel-precompute asset's sequence content key for row {} mismatches \
                 the fetched bytes (the sources drifted from the panel the asset was \
                 built against - loud, never a silent divergence)",
                row.path
            ),
        )?;
        let seq = Arc::new(seq);
        self.seqs.lock().unwrap().insert(row_index, seq.clone());
        Ok(seq)
    }

    /// THE SERVE (the masking commute lemma in code): iterate the
    /// window's member list in order, apply the SAME per-member
    /// predicates the live build applies, and group by the asset's
    /// recorded fold ids — the served fold order is first-appearance-
    /// among-SURVIVORS exactly as the live build's is, the fold
    /// content is derived through the same one-spelling helpers from
    /// the representative row, and the representative row's sequence
    /// is fetched live and FNV-verified.
    #[allow(clippy::too_many_arguments)]
    fn serve_window_folds(
        &self,
        rows: &[RowFold],
        window: usize,
        members: &[MemberRow],
        panel_component: &str,
        axis_start: u64,
        axis_end: u64,
        flank_overlap: Option<FlankOverlap>,
        exclude_paths: &[String],
        k: u64,
        fetch: &dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
        path_names: &[String],
    ) -> io::Result<(Vec<Fold>, u64, u64)> {
        let asset_window = self
            .windows
            .get(window)
            .ok_or_else(|| invalid("the window is outside the panel-precompute asset"))?;
        ensure(
            asset_window.member_fold.len() == members.len(),
            &format!(
                "the panel-precompute asset's member-fold map for window {window} does \
                 not match the partition maps (content keys drifted)",
            ),
        )?;
        let mut folds: Vec<Fold> = Vec::new();
        let mut served_to_live: Vec<Option<usize>> = vec![None; asset_window.fold_rep.len()];
        let mut excluded_rows: u64 = 0;
        let mut flank_excluded_rows: u64 = 0;
        for (member_index, member) in members.iter().enumerate() {
            if let Some(flank) = flank_overlap {
                if member.path_name != panel_component
                    && !row_covers_flanks(member.start, member.end, axis_start, axis_end, flank)
                {
                    flank_excluded_rows += 1;
                    continue;
                }
            }
            if row_is_excluded(&member.path_name, panel_component, exclude_paths) {
                excluded_rows += 1;
                continue;
            }
            let served_fold = asset_window.member_fold[member_index];
            ensure(
                served_fold < asset_window.fold_rep.len(),
                "the panel-precompute asset carries an out-of-range fold id",
            )?;
            let live_index = match served_to_live[served_fold] {
                Some(index) => index,
                None => {
                    let rep_row = asset_window.fold_rep[served_fold];
                    let row = rows
                        .get(rep_row)
                        .ok_or_else(|| invalid("the panel-precompute asset names an absent row"))?;
                    let seq = self.seq_of(rep_row, rows, fetch, path_names)?;
                    let relative = contained_relative_walk(row, k);
                    let node_positions = fold_node_positions(row, k);
                    let (nodes, edges) = fold_multisets(&relative);
                    folds.push(Fold {
                        seq: (*seq).clone(),
                        len: (row.end - row.start) as u64,
                        walk: relative,
                        node_positions,
                        members: Vec::new(),
                        nodes,
                        edges,
                    });
                    let index = folds.len() - 1;
                    served_to_live[served_fold] = Some(index);
                    index
                }
            };
            folds[live_index].members.push(member.clone());
        }
        Ok((folds, excluded_rows, flank_excluded_rows))
    }
}

/// THE PRECOMPUTE STEP (SS10.6.23 slice 1 — the once-per-panel price,
/// paid at panel-prep or by the standalone step, NEVER per sample): run
/// the SAME row-context build and the SAME unmasked fold-layout
/// grouping the live call runs (the one-spelling guarantee — this
/// builder IS the live derivation), serialize them into the panel
/// asset's new section, and receipt the price honestly (walls and
/// counts on stderr, beside the asset it builds). The asset keys on
/// CONTENT, never filenames: the panel identity object, the maps
/// digest, the routes digest, the component, the frame + origin, the
/// stored-frame convention, k.
pub fn precompute_panel_loci(
    panel_prefix: &str,
    routes_dir: &Path,
    partition_graphs: &Path,
    component: &str,
    global_frame: bool,
    record_origin: u64,
    out_dir: &Path,
) -> io::Result<()> {
    let started = Instant::now();
    ensure(
        record_origin == 0 || global_frame,
        "--record-origin requires --global-frame (the committed local frame is closed)",
    )?;
    let identity = PanelIdentity::read(panel_prefix)?;
    let panel = SyngIndex::load(panel_prefix, SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let graph_path = routes_dir.join("graph.json");
    let graph: routes::Graph = read_json(&graph_path)?;
    if graph.panel != identity || !graph.generation_complete {
        return Err(invalid("incompatible route inventory"));
    }
    let routes_digest = routes_graph_digest(&graph_path)?;
    let lanes = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect::<Vec<_>>();
    let sources = routes::Sources::open(&graph.source_paths, lanes)?;
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
    let source_of_name: HashMap<&str, usize> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.as_str(), lane.id))
        .collect();
    let lane_cache = LaneCache {
        lane_lengths: sources.lanes.iter().map(|&(_, len)| len).collect(),
        sources: std::sync::Arc::new(sources),
        whole_lane: !global_frame,
        seqs: (0..graph.lanes.len())
            .map(|_| std::sync::Mutex::new(None))
            .collect(),
        loads: AtomicU64::new(0),
        bytes: AtomicU64::new(0),
        nanos: AtomicU64::new(0),
        lane_trace: Default::default(),
    };
    let fetch_seq = |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
        let path_idx = *path_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the syng name map"))?;
        let path_len = panel.name_map.path_to_length[path_idx];
        let hi = hi.min(path_len);
        let source = *source_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the route sources"))?;
        if lo >= hi {
            return Ok(Vec::new());
        }
        lane_cache.fetch(source, lo, hi)
    };
    // The maps load + the frame translation — the call's own inputs
    // resolution, mirrored (the DERIVATIONS below are the one-spelling
    // functions; this resolution is I/O plumbing and states itself).
    let maps_digest = partition_maps_digest(partition_graphs)?;
    let mut maps: BTreeMap<u32, PartitionMap> = BTreeMap::new();
    let mut axis_rows: Vec<(u32, u64, u64)> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(partition_graphs)? {
        let path = entry?.path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        if name.ends_with(".gfa.map.json") {
            names.push(name);
        }
    }
    names.sort();
    for name in &names {
        let text = std::fs::read_to_string(partition_graphs.join(name))?;
        let map: PartitionMap = serde_json::from_str(&text)?;
        for axis in map
            .members
            .iter()
            .filter(|m| m.path_name == component)
        {
            axis_rows.push((map.partition, axis.start, axis.end));
        }
        maps.insert(map.partition, map);
    }
    axis_rows.sort_by_key(|&(partition, start, end)| (start, end, partition));
    let panel_component: String = if global_frame {
        let frame = territory_frame::TerritoryFrame::from_local_names(
            maps
                .values()
                .flat_map(|map| map.members.iter().map(|row| row.path_name.as_str())),
            record_origin,
        )?;
        let panel_component = frame.global_name(component)?;
        for map in maps.values_mut() {
            for row in &mut map.members {
                let translated =
                    territory_frame::translate_row(&row.path_name, row.start, row.end, record_origin)?;
                row.path_name = translated.path_name;
                row.start = translated.start;
                row.end = translated.end;
            }
        }
        for row in axis_rows.iter_mut() {
            let translated =
                territory_frame::translate_row(component, row.1, row.2, record_origin)?;
            row.1 = translated.start;
            row.2 = translated.end;
        }
        panel_component
    } else {
        component.to_string()
    };
    let window_partition: Vec<u32> =
        axis_rows.iter().map(|&(partition, _, _)| partition).collect();
    let window_axis: Vec<(u64, u64)> =
        axis_rows.iter().map(|&(_, start, end)| (start, end)).collect();
    // (The stored-frame convention, read the SAME way the call reads
    // it — the walks the asset serves depend on it.)
    let stored_frame = std::env::var("IMPG_REALIGN_STORED_WALK_SKELETON").is_ok();
    let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq, k, stored_frame);
    let rows_started = Instant::now();
    for &partition in &window_partition {
        let map = maps
            .get(&partition)
            .ok_or_else(|| invalid("axis partition map missing"))?;
        store.ensure_partition(&map.members, partition)?;
    }
    let rows_seconds = rows_started.elapsed().as_secs_f64();
    // THE ROW CONTEXT, serialized (the sequences replaced by their
    // FNV-64 CONTENT KEY — the bytes stay in the sources).
    let rows: Vec<PrecomputeRowAsset> = store
        .folds
        .iter()
        .map(|row| PrecomputeRowAsset {
            path_name: panel.name_map.path_to_name[row.path].clone(),
            start: row.start,
            end: row.end,
            seq_fnv: fnv1a64(&row.seq),
            stored_walk: row.stored_walk.clone(),
            walk: row.walk.clone(),
            frames: row.frames.clone(),
            members: row.members.clone(),
        })
        .collect();
    // THE FOLD LAYOUTS, per window (the live grouping, unmasked — the
    // masking filters are per-member predicates the serve applies at
    // its own call time; the grouping is the pure function recorded).
    let windows_started = Instant::now();
    let mut windows: Vec<PrecomputeWindowAsset> = Vec::with_capacity(window_partition.len());
    let mut folds_total: u64 = 0;
    for (window, &partition) in window_partition.iter().enumerate() {
        let map = maps
            .get(&partition)
            .ok_or_else(|| invalid("axis partition map missing"))?;
        let (axis_start, axis_end) = window_axis[window];
        let (_folds, fold_rep, member_fold, _, _) = build_window_folds_live(
            &store.folds,
            &store.dedup,
            &map.members,
            &panel_component,
            axis_start,
            axis_end,
            None,
            &[],
            &path_of_name,
            &|path_idx| panel.name_map.path_to_length[path_idx],
            k,
        )?;
        folds_total += fold_rep.len() as u64;
        windows.push(PrecomputeWindowAsset {
            partition,
            axis_start,
            axis_end,
            fold_rep,
            member_fold,
        });
    }
    let windows_seconds = windows_started.elapsed().as_secs_f64();
    // THE BIND SCAFFOLDING (SS10.6.25 slice 2 — the once-per-panel
    // price's new section): the canonical-scheme steps over the
    // territory's own member-row extents (the panel-knowable
    // coordinate set), merged and padded by the syncmer window's own
    // reach, derived through the SAME `canonical_steps_over_range`
    // the live bind runs — the builder IS the live derivation.
    let bind_started = Instant::now();
    let member_rows: Vec<(usize, u64, u64)> = maps
        .values()
        .flat_map(|map| map.members.iter())
        .map(|row| {
            let path = *path_of_name
                .get(&row.path_name)
                .ok_or_else(|| {
                    invalid(&format!(
                        "the partition maps name a member row absent from the syng: {:?}",
                        row.path_name
                    ))
                })?;
            Ok((path, row.start, row.end))
        })
        .collect::<io::Result<Vec<_>>>()?;
    let bind_paths = build_bind_scaffold(&panel, &fetch_seq, k, member_rows.into_iter())?;
    let bind_ranges: u64 = bind_paths.iter().map(|p| p.ranges.len() as u64).sum();
    let bind_steps: u64 = bind_paths.iter().map(|p| p.steps.len() as u64).sum();
    let bind_seconds = bind_started.elapsed().as_secs_f64();
    let manifest = PanelPrecomputeManifest {
        format: "impg-panel-precompute-loci/2".to_string(),
        panel: identity,
        k,
        component: panel_component,
        frame: if global_frame { "global" } else { "local" }.to_string(),
        record_origin,
        stored_frame,
        maps_digest,
        routes_digest,
        windows: windows.len() as u64,
        rows: rows.len() as u64,
        folds: folds_total,
        bind_paths: bind_paths.len() as u64,
        bind_ranges,
        bind_steps,
    };
    std::fs::create_dir_all(out_dir)?;
    let rows_bytes = encode_rows(&rows);
    let windows_bytes = encode_windows(&windows);
    let bind_bytes = encode_bind_scaffolding(&bind_paths);
    crate::genome_inference::atomic_write(&out_dir.join("rows.bin"), &rows_bytes)?;
    crate::genome_inference::atomic_write(&out_dir.join("windows.bin"), &windows_bytes)?;
    crate::genome_inference::atomic_write(
        &out_dir.join("bind-scaffolding.bin"),
        &bind_bytes,
    )?;
    crate::genome_inference::write_json(&out_dir.join("manifest.json"), &manifest)?;
    eprintln!(
        "[precompute] rows: {} distinct rows over {} partitions ({} windows) \
         [{rows_seconds:.1}s]",
        rows.len(),
        maps.len(),
        window_partition.len(),
    );
    eprintln!(
        "[precompute] windows: {} windows, {folds_total} folds [{windows_seconds:.1}s]",
        windows.len(),
    );
    eprintln!(
        "[precompute] bind scaffolding: {} paths, {} derivation ranges, {} steps \
         [{bind_seconds:.1}s]",
        bind_paths.len(),
        bind_ranges,
        bind_steps,
    );
    eprintln!(
        "[precompute] asset: {} (rows.bin {} bytes, windows.bin {} bytes, \
         bind-scaffolding.bin {} bytes, manifest.json; maps digest {:016x}, routes digest \
         {:016x}) [total {:.1}s]",
        out_dir.display(),
        rows_bytes.len(),
        windows_bytes.len(),
        bind_bytes.len(),
        maps_digest,
        routes_digest,
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

fn pack_edge(left: GraphNode, right: GraphNode) -> PackedEdge {
    (left as u64) << 32 | right as u64
}

fn multiset_overlap(left: &[(u64, u32)], right: &[(u64, u32)]) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                total += left[i].1 as u64 * right[j].1 as u64;
                i += 1;
                j += 1;
            }
        }
    }
    total as f64
}

fn merged_multiset(left: &[(u64, u32)], right: &[(u64, u32)]) -> Vec<(u64, u32)> {
    let (mut i, mut j) = (0usize, 0usize);
    let mut merged: Vec<(u64, u32)> = Vec::with_capacity(left.len() + right.len());
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                merged.push(left[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(right[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                merged.push((left[i].0, left[i].1 + right[j].1));
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&left[i..]);
    merged.extend_from_slice(&right[j..]);
    merged
}

fn differing_observed_mass(
    left: &[(u64, u32)],
    right: &[(u64, u32)],
    observed: &dyn Fn(u64) -> f64,
) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0.0f64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                total += observed(left[i].0) * left[i].1 as f64;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                total += observed(right[j].0) * right[j].1 as f64;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let delta = (left[i].1 as i64 - right[j].1 as i64).unsigned_abs();
                if delta > 0 {
                    total += observed(left[i].0) * delta as f64;
                }
                i += 1;
                j += 1;
            }
        }
    }
    while i < left.len() {
        total += observed(left[i].0) * left[i].1 as f64;
        i += 1;
    }
    while j < right.len() {
        total += observed(right[j].0) * right[j].1 as f64;
        j += 1;
    }
    total
}

/// THE SAME WALK WITH THE PER-KEY OBSERVED MASS READ FROM PARALLEL
/// ARRAYS — the committed certificate's dominant residual: one HashMap
/// probe per merged key per undecided test (~1,200 probes per test,
/// 13.6M tests at T00048 locus 2, ~1.0 us of each ~1.2 us test).
///
/// THE ARRAYS HOLD `observed(key)` — THE UNIT — AND NOT THE PRODUCT
/// `observed(key) * count as f64`. That is deliberate, and it must stay
/// that way: this walk multiplies by `count` on the less/greater/tail
/// branches but by `delta` (the count difference) on the equal branch,
/// so a stored product would reproduce three of the four products and
/// silently perturb the equal branch's bits. Holding the unit keeps
/// EVERY product the committed expression on identical operands, so
/// the total is the same f64 as `differing_observed_mass`'s. Do not
/// "optimise" this into a stored product.
///
/// The two key lists are walked in the committed order with the
/// committed `delta > 0` skip; only the probe is replaced (an array
/// read), so the accumulation sequence and every operand are the
/// committed ones.
fn differing_observed_mass_from_units(
    left: &[(u64, u32)],
    left_observed: &[f64],
    right: &[(u64, u32)],
    right_observed: &[f64],
) -> f64 {
    debug_assert_eq!(left.len(), left_observed.len());
    debug_assert_eq!(right.len(), right_observed.len());
    let (mut i, mut j, mut total) = (0usize, 0usize, 0.0f64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                total += left_observed[i] * left[i].1 as f64;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                total += right_observed[j] * right[j].1 as f64;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let delta = (left[i].1 as i64 - right[j].1 as i64).unsigned_abs();
                if delta > 0 {
                    total += left_observed[i] * delta as f64;
                }
                i += 1;
                j += 1;
            }
        }
    }
    while i < left.len() {
        total += left_observed[i] * left[i].1 as f64;
        i += 1;
    }
    while j < right.len() {
        total += right_observed[j] * right[j].1 as f64;
        j += 1;
    }
    total
}

// ---------------------------------------------------------------------------
// The alignment-connected fold relaxation (the segdup follow-up's
// stage-1 rule, env IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS).
//
// THE EXACT FOLD RULE coalesces rows only on identical spelled
// sequence AND identical contained relative walk — exact node+edge
// multiset equality at identical offsets. At near-identical segdup
// copies that rule is too strict: the truth's second haplotype
// exists as graph paths (its rows are present in the fold space)
// but the variant pockets — the node-multiset deltas between the
// near-identical copies — break the exact multiset, so the copies
// never coalesce and the exact correspondence refuses.
//
// THE RELAXED RULE. Rows fold by their connection in the locality's
// alignment structure — the pggb bubble: the connected alignment
// component containing the axis correspondence (two folds are
// alignment-adjacent iff their contained node multisets share a
// node; the bubble is the transitive closure around the axis fold)
// — rather than exact multiset equality alone, SUBJECT TO the
// safety property that preserves the fold's guarantee:
//
// A FOLD STILL MEANS NO CANDIDATE CAN DISTINGUISH THE MEMBERS. The
// members must be indistinguishable AT READ LEVEL — identical on
// the observed-material projection: every evidence unit's
// per-fold log-likelihood is bit-identical across the members
// (the operative, exactly-verified check; it is the definition of
// "no candidate distinguishes the members", computed from the
// scoring matrix itself), whose data-level form is: every node and
// edge in the symmetric difference of the members' contained
// multisets — the variant pockets — carries ZERO observed census
// mass at this locus. A pocket carrying read mass means some read's
// occurrence covers differing material, its placements vote against
// the members' different spellings, and its unit's likelihood
// differs — the fold does NOT happen (the merge is refused and the
// exact folds stand). This is the property that prevents folding
// genuinely-distinguishable material: an SNV pocket with even one
// observed read refuses the merge, whatever the node-sharing.
// (The conditions are transitive: pairwise massless symmetric
// differences union to a massless difference, and bit-identical
// columns are an equivalence — so a merged group is safe at every
// pair by construction.)
//
// The merged fold carries the representative's sequence, walk and
// pin skeleton (any member is a valid representative under the
// invariant — the asserted bit-identical columns ARE the proof),
// every member's member rows, and the additive multiset merge of
// the members' node/edge multisets (the pair machinery's own
// convention). Merging cannot change any class log-likelihood
// value (LL-identical columns map every class to a class of the
// same value), so the winner's likelihood and the truth pair's
// likelihood are preserved; class COUNTS, fold indices, called
// strain lists and QUAL spectra can change where merges land —
// the receipts name every merge group (the report's
// alignment_connected_folds field) so any changed call is
// attributable.
// ---------------------------------------------------------------------------

/// The per-locus merge anatomy, emitted in the report under the
/// relaxed rule (checker input: every merge named, every pocket
/// mass stated).
struct MergeGroupReceipt {
    /// The merged fold's index in the relaxed fold space.
    index: usize,
    /// The representative's index in the exact fold space (whose
    /// sequence/walk/skeleton the merged fold carries).
    representative: usize,
    /// The representative's first member row (the checker's
    /// representative-class identification input).
    representative_row: MemberRow,
    /// Every member's index in the exact fold space, ascending.
    members: Vec<usize>,
    /// The pairwise symmetric-difference node/edge observed mass
    /// (the variant pockets' read mass — zero at every admitted
    /// merge by the safety condition).
    pocket_node_mass: f64,
    pocket_edge_mass: f64,
}

/// The merge pass's honest counters (the receipt's refusal anatomy).
struct MergeCounters {
    /// Candidate pairs whose columns are bit-identical (the safety
    /// condition's input set).
    equal_column_pairs: u64,
    /// Refusals by reason (mutually exclusive, in check order).
    refused_outside_bubble: u64,
    refused_unshared: u64,
    refused_pocket_mass: u64,
    /// Admitted merges.
    merged_pairs: u64,
}

/// A minimal union-find over fold indices.
struct FoldUnionFind {
    parent: Vec<usize>,
}

impl FoldUnionFind {
    fn new(n: usize) -> Self {
        FoldUnionFind {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, index: usize) -> usize {
        let mut root = index;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut node = index;
        while self.parent[node] != root {
            let next = self.parent[node];
            self.parent[node] = root;
            node = next;
        }
        root
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

/// FNV-1a over the column's f64 bit patterns (the equal-column
/// bucket key; buckets are re-verified bit-exactly before use).
/// Whether two sorted contained multisets share at least one node
/// key (the alignment structure's adjacency: a shared interned node
/// IS an alignment edge between the two rows).
fn share_a_node(left: &[(u64, u32)], right: &[(u64, u32)]) -> bool {
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => return true,
        }
    }
    false
}

/// The multiset-intersection mass of two sorted node multisets (the
/// alignment correspondence measure: the shared material counted at
/// copy multiplicity, no thresholds).
fn shared_node_mass(left: &[(u64, u32)], right: &[(u64, u32)]) -> u64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                total += left[i].1.min(right[j].1) as u64;
                i += 1;
                j += 1;
            }
        }
    }
    total
}

/// THE MERGE PASS (pure; unit-tested). Consumes the exact-rule fold
/// space with its scoring matrix and territory images, returns the
/// relaxed fold space, the remapped axis fold, the merge anatomy and
/// the counters. With no axis fold or no eligible pair the output
/// is the input unchanged (the exact rule stands).
#[allow(clippy::type_complexity)]
fn alignment_connected_merge(
    folds: Vec<Fold>,
    matrix: ScoreMatrix,
    territory: Vec<(u64, u64)>,
    axis_fold: Option<usize>,
    territory_axis_fold: Option<usize>,
    observed_node: &HashMap<u64, f64>,
    observed_edge: &HashMap<u64, f64>,
) -> (
    Vec<Fold>,
    ScoreMatrix,
    Vec<(u64, u64)>,
    Option<usize>,
    Vec<MergeGroupReceipt>,
    MergeCounters,
) {
    let n = folds.len();
    let mut groups: Vec<MergeGroupReceipt> = Vec::new();
    let counters = MergeCounters {
        equal_column_pairs: 0,
        refused_outside_bubble: 0,
        refused_unshared: 0,
        refused_pocket_mass: 0,
        merged_pairs: 0,
    };
    let Some(axis) = axis_fold else {
        return (
            folds,
            matrix,
            territory,
            territory_axis_fold,
            groups,
            counters,
        );
    };
    if n < 2 {
        return (
            folds,
            matrix,
            territory,
            territory_axis_fold,
            groups,
            counters,
        );
    }
    // The alignment structure's connected components: one union per
    // shared node (the pggb bubble around the axis correspondence is
    // the axis fold's component).
    let mut node_carriers: HashMap<u64, Vec<usize>> = HashMap::new();
    for (index, fold) in folds.iter().enumerate() {
        for &(key, _) in &fold.nodes {
            node_carriers.entry(key).or_default().push(index);
        }
    }
    let mut components = FoldUnionFind::new(n);
    for carriers in node_carriers.values() {
        for &index in carriers.iter().skip(1) {
            components.union(carriers[0], index);
        }
    }
    let mut components = components;
    let axis_root = components.find(axis);
    // Component roots, hoisted (the parallel pair scan reads them).
    let roots: Vec<usize> = (0..n).map(|i| components.find(i)).collect();
    // The equal-column buckets (the safety condition's candidate
    // set: only bit-identical columns can be indistinguishable).
    // (Candidate C: the hash is the row's canonical VALUE hash — the
    // committed fnv-over-bits, read through the compressed form — so
    // two equal columns hash equal whichever encoding each row took.)
    let mut buckets: HashMap<u64, Vec<usize>> = HashMap::new();
    for index in 0..n {
        buckets.entry(matrix.ll_bits_hash(index)).or_default().push(index);
    }
    let node_mass = |key: u64| observed_node.get(&key).copied().unwrap_or(0.0);
    let edge_mass = |key: u64| observed_edge.get(&key).copied().unwrap_or(0.0);
    // The eligible-pair scan, parallel over buckets (folds, matrix
    // and roots are read-only here; the per-bucket counters sum to
    // the receipt's).
    struct BucketOutcome {
        pairs: Vec<(usize, usize)>,
        equal_column_pairs: u64,
        refused_outside_bubble: u64,
        refused_unshared: u64,
        refused_pocket_mass: u64,
    }
    let outcomes: Vec<BucketOutcome> = buckets
        .values()
        .collect::<Vec<_>>()
        .par_iter()
        .map(|bucket| {
            let mut out = BucketOutcome {
                pairs: Vec::new(),
                equal_column_pairs: 0,
                refused_outside_bubble: 0,
                refused_unshared: 0,
                refused_pocket_mass: 0,
            };
            for a in 0..bucket.len() {
                for b in (a + 1)..bucket.len() {
                    let (i, j) = (bucket[a].min(bucket[b]), bucket[a].max(bucket[b]));
                    // (Hash buckets can collide: verify the columns
                    // bit-exactly — the safety condition's exact
                    // predicate, never the hash.)
                    if !matrix.ll_bits_equal(i, j) {
                        continue;
                    }
                    out.equal_column_pairs += 1;
                    if roots[i] != axis_root || roots[j] != axis_root {
                        out.refused_outside_bubble += 1;
                        continue;
                    }
                    if !share_a_node(&folds[i].nodes, &folds[j].nodes) {
                        out.refused_unshared += 1;
                        continue;
                    }
                    // THE SAFETY PROPERTY (the derived condition): the
                    // variant pockets — the symmetric-difference node
                    // and edge material — must carry ZERO observed
                    // read mass, or the fold does NOT happen.
                    let pocket_nodes =
                        differing_observed_mass(&folds[i].nodes, &folds[j].nodes, &node_mass);
                    let pocket_edges =
                        differing_observed_mass(&folds[i].edges, &folds[j].edges, &edge_mass);
                    if pocket_nodes > 0.0 || pocket_edges > 0.0 {
                        out.refused_pocket_mass += 1;
                        continue;
                    }
                    out.pairs.push((i, j));
                }
            }
            out
        })
        .collect();
    let mut merges = FoldUnionFind::new(n);
    let mut counters = counters;
    for outcome in outcomes {
        counters.equal_column_pairs += outcome.equal_column_pairs;
        counters.refused_outside_bubble += outcome.refused_outside_bubble;
        counters.refused_unshared += outcome.refused_unshared;
        counters.refused_pocket_mass += outcome.refused_pocket_mass;
        counters.merged_pairs += outcome.pairs.len() as u64;
        for (i, j) in outcome.pairs {
            merges.union(i, j);
        }
    }
    // The merge groups: every fold's merge root, grouped and sorted
    // (the representative is the group's smallest old index —
    // deterministic).
    let mut groups_by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    for index in 0..n {
        let root = merges.find(index);
        groups_by_root.entry(root).or_default().push(index);
    }
    let mut sorted_groups: Vec<Vec<usize>> = groups_by_root
        .into_values()
        .filter(|members| members.len() > 1)
        .map(|mut members| {
            members.sort_unstable();
            members
        })
        .collect();
    sorted_groups.sort_unstable_by_key(|members| members[0]);
    let group_of_min: HashMap<usize, usize> = sorted_groups
        .iter()
        .enumerate()
        .map(|(group_index, members)| (members[0], group_index))
        .collect();
    let member_of: HashMap<usize, usize> = sorted_groups
        .iter()
        .enumerate()
        .flat_map(|(group_index, members)| {
            members
                .iter()
                .map(move |&member| (member, group_index))
        })
        .collect();
    // Rebuild the fold space: representatives (and solo folds) in old
    // index order; non-representative members fold into their group.
    let mut folds: Vec<Option<Fold>> = folds.into_iter().map(Some).collect();
    let mut new_folds: Vec<Fold> = Vec::with_capacity(n);
    // The new index of each old fold (group members share the merged
    // fold's index) and the representative of each new fold.
    let mut remap: Vec<usize> = vec![0; n];
    let mut new_representatives: Vec<usize> = Vec::with_capacity(n);
    for old in 0..n {
        if member_of.contains_key(&old) && !group_of_min.contains_key(&old) {
            // A non-representative group member: it folds into its
            // group's merged fold (remapped after the loop).
            continue;
        }
        remap[old] = new_folds.len();
        new_representatives.push(old);
        if let Some(&group_index) = group_of_min.get(&old) {
            let members = &sorted_groups[group_index];
            let representative =
                folds[old].take().expect("representative taken twice");
            let mut all_members: Vec<MemberRow> = Vec::new();
            let mut pocket_node_mass = 0.0f64;
            let mut pocket_edge_mass = 0.0f64;
            for &member in members.iter().skip(1) {
                // (The admitted-pair transitivity proof, asserted: every
                // member's column is bit-identical to the representative's.)
                assert!(
                    matrix.ll_bits_equal(member, old),
                    "merged group with distinguishable members"
                );
                let member_fold =
                    folds[member].take().expect("group member taken twice");
                all_members.extend(member_fold.members.iter().cloned());
                pocket_node_mass += differing_observed_mass(
                    &representative.nodes,
                    &member_fold.nodes,
                    &node_mass,
                );
                pocket_edge_mass += differing_observed_mass(
                    &representative.edges,
                    &member_fold.edges,
                    &edge_mass,
                );
            }
            all_members.extend(representative.members.iter().cloned());
            groups.push(MergeGroupReceipt {
                index: new_folds.len(),
                representative: old,
                representative_row: representative.members[0].clone(),
                members: members.clone(),
                pocket_node_mass,
                pocket_edge_mass,
            });
            // (The fold's material stays ONE hypothesis's — the
            // representative's multisets, the committed convention:
            // under the exact rule a fold's members do not multiply
            // its material either. The other members are PROVEN
            // equivalent, not additional copies.)
            new_folds.push(Fold {
                seq: representative.seq,
                len: representative.len,
                walk: representative.walk,
                node_positions: representative.node_positions,
                members: all_members,
                nodes: representative.nodes,
                edges: representative.edges,
            });
        } else {
            new_folds.push(folds[old].take().expect("solo fold taken twice"));
        }
    }
    // Every group member's remap lands on the merged fold's index.
    for members in &sorted_groups {
        let merged_index = remap[members[0]];
        for &member in members {
            remap[member] = merged_index;
        }
    }
    // The per-fold rows and territory images regroup onto each new
    // fold's representative (the merged fold carries the
    // representative's column — asserted bit-identical to every
    // member's above).
    // The per-fold rows regroup onto each new fold's representative
    // (the merged fold carries the representative's column — asserted
    // bit-identical to every member's above), at the row's own
    // encoding (the committed clone, compressed).
    let new_rows: Vec<ScoreRow> = new_representatives
        .iter()
        .map(|&old| matrix.row_clone(old))
        .collect();
    let new_matrix = ScoreMatrix {
        n_units: matrix.n_units,
        elsewhere: matrix.elsewhere,
        rows: new_rows,
        unit_fold: None,
    };
    let new_territory: Vec<(u64, u64)> =
        new_representatives.iter().map(|&old| territory[old]).collect();
    let new_territory_axis_fold = territory_axis_fold.map(|index| remap[index]);
    (
        new_folds,
        new_matrix,
        new_territory,
        new_territory_axis_fold,
        groups,
        counters,
    )
}

/// THE ALIGNMENT-CONNECTED TRUTH-SECOND CORRESPONDENCE (the
/// relaxation's attribution clause). The committed exact rule derives
/// the truth-second fold only at single-window partitions and only
/// when the haplotype's rows form exactly one fold — at the segdup
/// block the partitions are multi-window and the haplotype's tiling
/// is fragmented into several folds, so the exact rule refuses and
/// the locus is scored with the truth pair INEXPRESSIBLE even though
/// the truth's second haplotype exists as graph paths. Under the
/// relaxation the correspondence is read from the alignment structure
/// itself: among the folds carrying the second haplotype's member
/// rows, the fold connected to the axis correspondence by the
/// MAXIMUM shared contained-node multiset mass (the alignment's own
/// correspondence measure — the same connection the fold relation
/// reads; an argmax, not a threshold). Refused (None) when the axis
/// fold is absent, when no fold carries the haplotype, when the best
/// candidate shares NO node with the axis fold (not connected to the
/// bubble), or when the maximum is TIED (an ambiguous correspondence
/// — the no-fabrication discipline preserved). The receipt value
/// names every candidate and its shared mass (checker input).
fn connected_truth_second(
    folds: &[Fold],
    truth_first: Option<usize>,
    second_name: &str,
) -> (Option<usize>, Option<serde_json::Value>) {
    let Some(axis) = truth_first else {
        return (None, None);
    };
    let mut candidates: Vec<(usize, u64)> = folds
        .iter()
        .enumerate()
        .filter(|(_, fold)| {
            fold.members.iter().any(|m| m.path_name == second_name)
        })
        .map(|(index, fold)| {
            (index, shared_node_mass(&folds[axis].nodes, &fold.nodes))
        })
        .collect();
    candidates.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let receipt = Some(json!({
        "rule": "alignment-connected-argmax",
        "axis_fold": axis,
        "candidates": candidates.iter().map(|&(index, shared)| json!({
            "fold": index,
            "shared_node_mass": shared,
        })).collect::<Vec<_>>(),
    }));
    if candidates.is_empty() {
        return (None, receipt);
    }
    let best = candidates[0].1;
    if best == 0 {
        // Not connected to the axis correspondence at all: the
        // haplotype's folds lie outside the bubble.
        return (None, receipt);
    }
    let ties = candidates.iter().filter(|&&(_, shared)| shared == best).count();
    if ties > 1 {
        // An ambiguous correspondence: the no-fabrication rule.
        return (None, receipt);
    }
    (Some(candidates[0].0), receipt)
}

fn signature_cosine_distance(
    nodes_a: &[(u64, u32)],
    edges_a: &[(u64, u32)],
    nodes_b: &[(u64, u32)],
    edges_b: &[(u64, u32)],
) -> f64 {
    let norm = |nodes: &[(u64, u32)], edges: &[(u64, u32)]| {
        nodes.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
            + edges.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
    };
    let (norm_a, norm_b) = (norm(nodes_a, edges_a), norm(nodes_b, edges_b));
    if norm_a == 0.0 && norm_b == 0.0 {
        return 0.0;
    }
    if norm_a == 0.0 || norm_b == 0.0 {
        return 1.0;
    }
    let overlap = multiset_overlap(nodes_a, nodes_b) + multiset_overlap(edges_a, edges_b);
    1.0 - overlap / (norm_a * norm_b).sqrt()
}

fn class_physical_pairs(same_row: bool, members_a: usize, members_b: usize) -> u64 {
    if same_row {
        members_a as u64 * (members_a as u64 + 1) / 2
    } else {
        members_a as u64 * members_b as u64
    }
}

fn median_of(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).expect("finite relative jumps"));
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    }
}

struct KneeCut {
    has_knee: bool,
    cut: f64,
}

fn spectrum_knee(sorted: &[f64]) -> KneeCut {
    let positives: Vec<f64> = sorted.iter().copied().filter(|&d| d > 0.0).collect();
    if positives.len() < 2 {
        return KneeCut { has_knee: false, cut: 0.0 };
    }
    let jumps: Vec<f64> = positives
        .windows(2)
        .map(|pair| (pair[1] - pair[0]) / pair[1])
        .collect();
    let mut center = jumps.clone();
    let typical = median_of(&mut center);
    let (mut best, mut best_index) = (f64::NEG_INFINITY, 0usize);
    for (index, &jump) in jumps.iter().enumerate() {
        if jump >= best {
            best = jump;
            best_index = index;
        }
    }
    if best > typical {
        KneeCut { has_knee: true, cut: positives[best_index] }
    } else {
        KneeCut { has_knee: false, cut: 0.0 }
    }
}

fn qual_p(s_win: f64, called_classes: usize, alternative: Option<f64>) -> Option<f64> {
    if called_classes == 0 {
        return None;
    }
    let alternative = alternative.unwrap_or(0.0);
    let denominator = called_classes as f64 * s_win + alternative;
    (denominator > 0.0).then(|| s_win / denominator)
}

fn qual_from_p(p: f64) -> Option<f64> {
    (p < 1.0).then(|| -10.0 * (1.0 - p).log10())
}

struct ClusterQual {
    knee: Option<f64>,
    cut: f64,
    shape: &'static str,
    cluster_size: usize,
    k: usize,
    excluded: Vec<bool>,
    alternative: Option<f64>,
    p: Option<f64>,
}

/// The tie certificate's exact input: per tied called class (beyond
/// the winner), either the FULL band (its observed-mass distance to
/// every spectrum entry — the emitted form at every non-degenerate
/// locus) or, at a fully degenerate locus (every class ties at the
/// pure-E score), the compact CUT-INDEX LIST — the band's exact
/// <=cut predicate set. The two forms are EXCHANGEABLE BY
/// CONSTRUCTION: cluster_form_qual consumes from each band only the
/// predicate {i : band[i] <= cut} (the excluded-union and the
/// union-find edges — band values above the cut are never read), and
/// the cut-index list is exactly that predicate set, so k, excluded,
/// alternative, cluster_size, shape and p are identical under both.
enum TieCertificate<'a> {
    /// The naive full-band form — the ORACLE's own shape (the unit
    /// tests construct it verbatim; the product path never does).
    #[allow(dead_code)]
    Bands(&'a [(usize, Vec<f64>)]),
    CutIndices(&'a [(usize, CutSet)]),
}

impl TieCertificate<'_> {
    fn len(&self) -> usize {
        match self {
            TieCertificate::Bands(entries) => entries.len(),
            TieCertificate::CutIndices(entries) => entries.len(),
        }
    }
}

fn cluster_form_qual(
    s_win: f64,
    spectrum: &[(f64, f64)],
    tied: &TieCertificate,
) -> ClusterQual {
    let mut sorted: Vec<f64> = spectrum.iter().map(|&(distance, _)| distance).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
    let knee_cut = spectrum_knee(&sorted);
    let cut = knee_cut.cut;
    let mut excluded = vec![false; spectrum.len()];
    let mut band = 0usize;
    for (index, &(distance, _)) in spectrum.iter().enumerate() {
        if distance <= cut {
            excluded[index] = true;
            band += 1;
        }
    }
    let mut parent: Vec<usize> = (0..=tied.len()).collect();
    let root = |parent: &mut Vec<usize>, node: usize| {
        let mut current = node;
        while parent[current] != current {
            current = parent[current];
        }
        let mut walked = node;
        while parent[walked] != current {
            let next = parent[walked];
            parent[walked] = current;
            walked = next;
        }
        current
    };
    match tied {
        TieCertificate::Bands(entries) => {
            for (_, band_distances) in entries.iter() {
                for (index, &distance) in band_distances.iter().enumerate() {
                    if distance <= cut {
                        excluded[index] = true;
                    }
                }
            }
            for (left, (index_left, _)) in entries.iter().enumerate() {
                if spectrum[*index_left].0 <= cut {
                    let (a, b) = (root(&mut parent, 0), root(&mut parent, left + 1));
                    if a != b {
                        parent[a] = b;
                    }
                }
                for (right, (index_right, _)) in entries.iter().enumerate().skip(left + 1) {
                    if entries[left].1[*index_right] <= cut {
                        let (a, b) =
                            (root(&mut parent, left + 1), root(&mut parent, right + 1));
                        if a != b {
                            parent[a] = b;
                        }
                    }
                }
            }
        }
        TieCertificate::CutIndices(entries) => {
            // The compact form: each list IS the band's <=cut predicate
            // set, so the excluded-union is list membership and a
            // union-find edge (left, right) exists iff tied class
            // right's spectrum index appears in left's list. The
            // edge set is identical to the band form's (the band
            // path tests band_left[index_right] <= cut), so the
            // components — and k — are identical; the iteration is
            // list-driven, linear in the total predicate count
            // instead of quadratic in the tied-class count (the
            // degenerate loci tie ~35,510 classes — the quadratic
            // pairwise loop is the checker-side wall).
            let mut ordinal_of: HashMap<usize, usize> = HashMap::with_capacity(entries.len());
            for (ordinal, &(index, _)) in entries.iter().enumerate() {
                ordinal_of.insert(index, ordinal);
            }
            for (left, &(index, ref list)) in entries.iter().enumerate() {
                for member in list.iter() {
                    excluded[member] = true;
                }
                if spectrum[index].0 <= cut {
                    let (a, b) = (root(&mut parent, 0), root(&mut parent, left + 1));
                    if a != b {
                        parent[a] = b;
                    }
                }
                for member in list.iter() {
                    if let Some(&right) = ordinal_of.get(&member) {
                        let (a, b) =
                            (root(&mut parent, left + 1), root(&mut parent, right + 1));
                        if a != b {
                            parent[a] = b;
                        }
                    }
                }
            }
        }
    }
    let mut roots: Vec<usize> = (0..=tied.len()).map(|node| root(&mut parent, node)).collect();
    roots.sort_unstable();
    roots.dedup();
    let k = roots.len();
    let mut alternative: Option<f64> = None;
    for (index, &(_, score)) in spectrum.iter().enumerate() {
        if !excluded[index] {
            alternative = Some(alternative.map_or(score, |value| value.max(score)));
        }
    }
    let shape = if spectrum.is_empty() {
        "single_class"
    } else if knee_cut.has_knee {
        "knee"
    } else {
        "no_knee"
    };
    let p = qual_p(s_win, k, alternative);
    ClusterQual {
        knee: knee_cut.has_knee.then_some(cut),
        cut,
        shape,
        cluster_size: 1 + band,
        k,
        excluded,
        alternative,
        p,
    }
}

/// THE BAND'S `<=cut` PREDICATE SET, PACKED.
///
/// One BIT per spectrum position (8 B -> 1 bit per member, ~64x). At
/// the structurally degenerate shape this code already detects
/// (`called.len() == pair_count`) a locus carries `|called| - 1` sets
/// over a `pair_count`-long spectrum, so the unpacked form is quadratic
/// in the pair count with a 8-byte constant - the receipts-phase
/// audit's tens-of-GB class. Packed, the same content is 1/64th of it.
///
/// ITERATION IS ASCENDING, which is the order the producer pushes (it
/// walks the spectrum in position order) and therefore the exact order
/// the receipt lists the indices in - so the emitted bytes are the
/// unpacked form's own bytes, and `cluster_form_qual`'s list-driven
/// walk sees the same sequence.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct CutSet {
    /// Position bits, little-endian within each word.
    words: Vec<u64>,
    /// The number of SET bits (not the capacity), so the receipt's
    /// entry counter stays the predicate-set cardinality.
    len: usize,
}

impl CutSet {
    fn insert(&mut self, position: usize) {
        let word = position / 64;
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        let mask = 1u64 << (position % 64);
        if self.words[word] & mask == 0 {
            self.words[word] |= mask;
            self.len += 1;
        }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The set's members, ASCENDING (the producer's and the receipt's
    /// own order).
    fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.words
            .iter()
            .enumerate()
            .flat_map(|(word_index, &word)| {
                let mut remaining = word;
                std::iter::from_fn(move || {
                    if remaining == 0 {
                        return None;
                    }
                    let bit = remaining.trailing_zeros() as usize;
                    remaining &= remaining - 1;
                    Some(word_index * 64 + bit)
                })
            })
    }

    /// Build from ascending positions (the oracle/tests' own shape).
    fn from_ascending<I: IntoIterator<Item = usize>>(positions: I) -> CutSet {
        let mut set = CutSet::default();
        for position in positions {
            set.insert(position);
        }
        set
    }
}

/// THE TIED-CLASS cut-index lists (the wall fix's compute core,
/// extracted for the unit-test oracle; the per-pair constant's own fix):
/// per tied called class (beyond the winner, in called order), its
/// spectrum index and the exact <=cut predicate set of its band.
/// Each (tied class, spectrum entry) distance is decided by the
/// SOUND two-sided O(1) observed-mass bounds — D >= |Sm_n(t) -
/// Sm_n(e)| + |Sm_e(t) - Sm_e(e)| and D <= Sm_n(t) + Sm_n(e) +
/// Sm_e(t) + Sm_e(e), where Sm is the observed mass of the
/// (merged, additive) multiset, computed from the PER-FOLD mass
/// signatures (node_mass/edge_mass) built ONCE per locus — with the
/// exact differing_observed_mass as the fallback for the undecided
/// middle band, so the emitted predicate set is IDENTICAL to the
/// naive band matrix's by construction (the oracle tests below
/// prove it at BOTH the degenerate and the non-degenerate shapes).
/// The upper-bound-first ordering makes the common all-zero-mass
/// case (both multisets observed-empty, D = 0) decidable without
/// any merge. THE PER-PAIR CONSTANT'S FIX (the fold-merge open
/// term, closed): the naive certificate materialized |called| - 1
/// bands x |spectrum| entries, EVERY entry a fresh pair of merged
/// multisets (per-pair vector construction) plus two observed-mass
/// walks (per-pair full-column walks) — the ~6.4 ms-per-class-pair
/// constant the pre-merge fold-heavy shapes receipted (the ledger's
/// T03226 l15 anatomy) — while the cluster machinery consumes from
/// each band ONLY its <=cut predicate set (cluster_form_qual's
/// own contract), so this form computes exactly that at O(1) per
/// entry with the exact walk as the fallback. It now serves BOTH
/// the degenerate path (every class ties at the pure-E score) and
/// the non-degenerate path (a strict tie subset) — the answers are
/// identical to the naive band matrix's at both, and the emitted
/// certificate is the compact cut-index list everywhere (the
/// qual_tied_bands receipt member is replaced by
/// qual_tied_cut_indices at non-degenerate tie loci, mirroring the
/// committed degenerate precedent).
fn tied_class_cut_lists<NF, EF>(
    folds: &[Fold],
    called: &[(usize, usize)],
    spectrum: &[(f64, f64, f64, usize, usize, usize)],
    observed_node: &NF,
    observed_edge: &EF,
    cut: f64,
) -> Vec<(usize, CutSet)>
where
    NF: Fn(u64) -> f64 + Sync,
    EF: Fn(u64) -> f64 + Sync,
{
    let mut position_of: HashMap<(usize, usize), usize> = HashMap::with_capacity(spectrum.len());
    for (position, &(_, _, _, first, second, _)) in spectrum.iter().enumerate() {
        position_of.insert((first, second), position);
    }
    let node_mass: Vec<f64> = folds
        .iter()
        .map(|fold| {
            fold.nodes
                .iter()
                .map(|&(key, count)| observed_node(key) * count as f64)
                .sum()
        })
        .collect();
    let edge_mass: Vec<f64> = folds
        .iter()
        .map(|fold| {
            fold.edges
                .iter()
                .map(|&(key, count)| observed_edge(key) * count as f64)
                .sum()
        })
        .collect();
    let entries: Vec<(usize, usize, f64, f64)> = spectrum
        .iter()
        .map(|&(_, _, _, first, second, _)| {
            (
                first,
                second,
                node_mass[first] + node_mass[second],
                edge_mass[first] + edge_mass[second],
            )
        })
        .collect();
    // THE CERTIFICATE IS COMPUTED ONCE PER DISTINCT CLASS SHAPE (the
    // receipts-phase census's finding): a tied class's predicate set is a
    // pure function of its MERGED multisets and mass signatures, so two
    // tied classes with identical merged (nodes, edges) produce identical
    // lists BY CONSTRUCTION. The degenerate locus makes this the whole
    // cost: 7,749 tied classes over 7,749 spectrum entries = 46.5M
    // entries computed 7,749 times, 17.55 s, when ONE list would do
    // (measured: T00048 locus 3, distinct_rows=1, folds=124). Group by a
    // canonical hash of the merged multisets and verify the group on
    // match, so a hash collision can never merge two different shapes.
    let canonical = |first: usize, second: usize| -> u64 {
        let mut hash = 0xcbf29ce484222325u64;
        for &(key, count) in merged_multiset(&folds[first].nodes, &folds[second].nodes).iter() {
            for byte in key.to_le_bytes().iter().chain(count.to_le_bytes().iter()) {
                hash ^= *byte as u64;
                hash = hash.wrapping_mul(0x100000001b3);
            }
        }
        for &(key, count) in merged_multiset(&folds[first].edges, &folds[second].edges).iter() {
            for byte in key.to_le_bytes().iter().chain(count.to_le_bytes().iter()) {
                hash ^= *byte as u64;
                hash = hash.wrapping_mul(0x100000001b3);
            }
        }
        hash
    };
    // (representative per distinct shape, in called order)
    let mut representative_of: HashMap<u64, usize> = HashMap::new();
    let mut shape_of: Vec<usize> = Vec::with_capacity(called.len().saturating_sub(1));
    for (ordinal, &(first, second)) in called[1..].iter().enumerate() {
        let key = canonical(first, second);
        let slot = *representative_of.entry(key).or_insert(ordinal);
        shape_of.push(slot);
    }
    let distinct: Vec<bool> = {
        let mut is_first = vec![false; called.len().saturating_sub(1)];
        for (ordinal, &slot) in shape_of.iter().enumerate() {
            if slot == ordinal {
                is_first[ordinal] = true;
            }
        }
        is_first
    };
    // THE SPECTRUM'S SHAPES ARE BUILT ONCE, BEFORE THE LOOP: the inner
    // loop below needs an entry's merged multisets for every (class,
    // entry) test the two O(1) bounds cannot decide, and each such test
    // used to rebuild them from the entry's fold pair — two
    // `merged_multiset` allocations plus the merge walk, per test. An
    // entry's merged multisets are a function of its (first, second)
    // fold pair ALONE, so they are built once per spectrum index here
    // and INDEXED below. A pure cache: the bound tests, their order, the
    // exact walk and every f64 operation on the same operands are
    // untouched, so the predicate sets are bit-identical.
    //
    // MEMORY: merged(i, j) is at most len(i) + len(j) entries of 12 B,
    // plus one 8 B unit observed mass per entry (see the Shape fields),
    // so the cache is at most 20 * (sum of the locus's fold multiset
    // entries) * folds bytes — at the degenerate locus (T00048 locus 2:
    // 74,381 fold multiset entries, 124 folds, 7,749 spectrum entries)
    // ~184 MB. It is built only when there is a tied class to evaluate
    // (`called[1..]` is the loop's domain), so a locus with no tie
    // never pays for it.
    //
    // THE UNIT OBSERVED MASS RIDES WITH THE MULTISETS (the residual
    // fix): the undecided tests' walk probed the observed map once per
    // merged key. `node_observed`/`edge_observed` are index-parallel to
    // `nodes`/`edges` and hold `observed(key)` — the unit, NOT the
    // product; the walk's equal branch multiplies by `delta` rather
    // than by `count` (the reason is spelled out at
    // `differing_observed_mass_from_units`). The probes still happen
    // here, ONCE per spectrum entry, and the walk that ran 13.6M times
    // reads arrays instead.
    struct Shape {
        nodes: Vec<(u64, u32)>,
        edges: Vec<(u64, u32)>,
        node_observed: Vec<f64>,
        edge_observed: Vec<f64>,
    }
    let shape_of_pair = |first: usize, second: usize| -> Shape {
        let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
        let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
        let node_observed: Vec<f64> =
            nodes.iter().map(|&(key, _)| observed_node(key)).collect();
        let edge_observed: Vec<f64> =
            edges.iter().map(|&(key, _)| observed_edge(key)).collect();
        Shape {
            nodes,
            edges,
            node_observed,
            edge_observed,
        }
    };
    let spectrum_shapes: Vec<Shape> = if called.len() > 1 {
        spectrum
            .iter()
            .map(|&(_, _, _, first, second, _)| shape_of_pair(first, second))
            .collect()
    } else {
        Vec::new()
    };
    let computed: Vec<Option<(usize, CutSet)>> = called[1..]
        .par_iter()
        .enumerate()
        .map(|(ordinal, &(first, second))| {
            if !distinct[ordinal] {
                return None;
            }
            let index = *position_of
                .get(&(first, second))
                .expect("called class missing from the winner's spectrum");
            let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
            let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
            // the class's own unit observed mass, one probe per merged
            // key ONCE (the walk below then reads these arrays)
            let node_observed: Vec<f64> =
                nodes.iter().map(|&(key, _)| observed_node(key)).collect();
            let edge_observed: Vec<f64> =
                edges.iter().map(|&(key, _)| observed_edge(key)).collect();
            let (tn, te) = (
                node_mass[first] + node_mass[second],
                edge_mass[first] + edge_mass[second],
            );
            let mut list = CutSet::default();
            for (position, &(_, _, gn, ge)) in entries.iter().enumerate() {
                let lower = (tn - gn).abs() + (te - ge).abs();
                if lower > cut {
                    continue;
                }
                let upper = tn + gn + te + ge;
                if upper <= cut {
                    list.insert(position);
                    continue;
                }
                let other = &spectrum_shapes[position];
                let distance = differing_observed_mass_from_units(
                    &nodes,
                    &node_observed,
                    &other.nodes,
                    &other.node_observed,
                ) + differing_observed_mass_from_units(
                    &edges,
                    &edge_observed,
                    &other.edges,
                    &other.edge_observed,
                );
                if distance <= cut {
                    list.insert(position);
                }
            }
            Some((index, list))
        })
        .collect();
    // replicate each representative's predicate set across its shape class;
    // each tied class keeps its OWN spectrum index (only the set is shared)
    let mut by_slot: HashMap<usize, CutSet> = HashMap::new();
    for (ordinal, value) in computed.into_iter().enumerate() {
        if let Some((_, list)) = value {
            by_slot.insert(ordinal, list);
        }
    }
    called[1..]
        .iter()
        .enumerate()
        .map(|(ordinal, &(first, second))| {
            let index = *position_of
                .get(&(first, second))
                .expect("called class missing from the winner's spectrum");
            let slot = shape_of[ordinal];
            let list = by_slot
                .get(&slot)
                .expect("the class shape's representative was not computed")
                .clone();
            (index, list)
        })
        .collect()
}

/// THE TIE CLUSTER WITHOUT MATERIALIZING THE BANDS.
///
/// The committed certificate builds, for every tied class, its full
/// predicate set over the spectrum — `|called| x |spectrum|` entries (46.5
/// million at the degenerate locus, 14.8 s, and O(f^4) in the fold count
/// because `called` is itself f(f+1)/2 there). But the cluster consumes
/// that structure through exactly two questions:
///
///   k         = the number of CONNECTED COMPONENTS of the graph on
///               {winner} + tied classes, with an edge wherever
///               d(L, R) <= cut        (the lists encode exactly this
///               relation: R's spectrum index appears in L's list iff
///               d(L,R) <= cut, and `differing_observed_mass` is symmetric)
///   excluded  = the union of the banks — winner band OR any tied
///               class's band
///
/// Both are answered with the predicates the certificate already uses (the
/// SOUND two-sided observed-mass bounds, the exact walk as the undecided
/// fallback) — no threshold, no approximation — and neither needs the lists
/// to exist. The component sweep skips any pair already connected (a
/// union-find whose second question is free), and `excluded` is filled from
/// a SHRINKING candidate set, so the fan-out collapses as banks cover the
/// spectrum.
///
/// NOTE the deliberate omission: the `qual_tied_cut_indices` receipt member
/// is the certificate's own enumeration and is no longer produced by this
/// path. k, the excluded union, the alternative and the QUAL — everything
/// the answer is a function of — are computed identically.
fn tie_cluster_fast<NF, EF>(
    folds: &[Fold],
    called: &[(usize, usize)],
    spectrum: &[(f64, f64, f64, usize, usize, usize)],
    observed_node: &NF,
    observed_edge: &EF,
) -> ClusterQual
where
    NF: Fn(u64) -> f64 + Sync,
    EF: Fn(u64) -> f64 + Sync,
{
    let mut sorted: Vec<f64> = spectrum.iter().map(|e| e.0).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
    let knee_cut = spectrum_knee(&sorted);
    let cut = knee_cut.cut;
    let mut excluded = vec![false; spectrum.len()];
    let mut band = 0usize;
    for (index, entry) in spectrum.iter().enumerate() {
        if entry.0 <= cut {
            excluded[index] = true;
            band += 1;
        }
    }
    let node_mass: Vec<f64> = folds
        .iter()
        .map(|fold| {
            fold.nodes
                .iter()
                .map(|&(key, count)| observed_node(key) * count as f64)
                .sum()
        })
        .collect();
    let edge_mass: Vec<f64> = folds
        .iter()
        .map(|fold| {
            fold.edges
                .iter()
                .map(|&(key, count)| observed_edge(key) * count as f64)
                .sum()
        })
        .collect();
    // the spectrum's own per-pair mass signatures (no multisets yet)
    let signatures: Vec<(f64, f64)> = spectrum
        .iter()
        .map(|&(_, _, _, first, second, _)| {
            (
                node_mass[first] + node_mass[second],
                edge_mass[first] + edge_mass[second],
            )
        })
        .collect();
    let n_tied = called.len().saturating_sub(1);
    let mut parent: Vec<usize> = (0..=n_tied).collect();
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            node = parent[node];
        }
        let mut walked = node;
        while parent[walked] != node {
            let next = parent[walked];
            parent[walked] = node;
            walked = next;
        }
        node
    }
    if n_tied > 0 {
        struct Shape {
            nodes: Vec<(u64, u32)>,
            edges: Vec<(u64, u32)>,
            tn: f64,
            te: f64,
        }
        let shapes: Vec<Shape> = (0..n_tied)
            .map(|ordinal| {
                let (first, second) = called[ordinal + 1];
                Shape {
                    nodes: merged_multiset(&folds[first].nodes, &folds[second].nodes),
                    edges: merged_multiset(&folds[first].edges, &folds[second].edges),
                    tn: node_mass[first] + node_mass[second],
                    te: edge_mass[first] + edge_mass[second],
                }
            })
            .collect();
        // THE BAND PREDICATE, by the certificate's own bounds: reject on the
        // lower bound, admit on the upper bound, walk exactly in between.
        let in_band = |a: &Shape, b: &Shape| -> bool {
            let lower = (a.tn - b.tn).abs() + (a.te - b.te).abs();
            if lower > cut {
                return false;
            }
            let upper = a.tn + b.tn + a.te + b.te;
            if upper <= cut {
                return true;
            }
            let node_lookup = |key: u64| observed_node(key);
            let edge_lookup = |key: u64| observed_edge(key);
            let distance = differing_observed_mass(&a.nodes, &b.nodes, &node_lookup)
                + differing_observed_mass(&a.edges, &b.edges, &edge_lookup);
            distance <= cut
        };
        let in_band_with_signature = |a: &Shape, signature: (f64, f64), index: usize| -> bool {
            let lower = (a.tn - signature.0).abs() + (a.te - signature.1).abs();
            if lower > cut {
                return false;
            }
            let upper = a.tn + signature.0 + a.te + signature.1;
            if upper <= cut {
                return true;
            }
            let (_, _, _, first, second, _) = spectrum[index];
            let other = Shape {
                nodes: merged_multiset(&folds[first].nodes, &folds[second].nodes),
                edges: merged_multiset(&folds[first].edges, &folds[second].edges),
                tn: signature.0,
                te: signature.1,
            };
            in_band(a, &other)
        };
        // (1) the winner's edges: a tied class joins the winner when ITS OWN
        // spectrum distance is at or inside the cut (the committed rule)
        let mut position_of: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::with_capacity(spectrum.len());
        for (position, &(_, _, _, first, second, _)) in spectrum.iter().enumerate() {
            position_of.insert((first, second), position);
        }
        for ordinal in 0..n_tied {
            let (first, second) = called[ordinal + 1];
            let position = position_of[&(first, second)];
            if spectrum[position].0 <= cut {
                let (a, b) = (
                    root(&mut parent, 0),
                    root(&mut parent, ordinal + 1),
                );
                if a != b {
                    parent[a] = b;
                }
            }
        }
        // (2) the class-class edges: skip any pair already connected - the
        // second question costs a find, not a distance
        for left in 0..n_tied {
            for right in (left + 1)..n_tied {
                let (a, b) = (
                    root(&mut parent, left + 1),
                    root(&mut parent, right + 1),
                );
                if a == b {
                    continue;
                }
                if in_band(&shapes[left], &shapes[right]) {
                    parent[a] = b;
                }
            }
        }
        // (3) the excluded union, from a SHRINKING candidate set: once a
        // position is covered, no later class needs to look at it
        let mut remaining: Vec<usize> = (0..spectrum.len())
            .filter(|&index| !excluded[index])
            .collect();
        for shape in &shapes {
            if remaining.is_empty() {
                break;
            }
            remaining.retain(|&index| {
                let covered = in_band_with_signature(shape, signatures[index], index);
                if covered {
                    excluded[index] = true;
                }
                !covered
            });
        }
    }
    let mut roots: Vec<usize> = (0..=n_tied).map(|node| root(&mut parent, node)).collect();
    roots.sort_unstable();
    roots.dedup();
    let k = roots.len();
    let mut alternative: Option<f64> = None;
    for (index, entry) in spectrum.iter().enumerate() {
        if !excluded[index] {
            alternative = Some(alternative.map_or(entry.1, |value| value.max(entry.1)));
        }
    }
    let shape = if spectrum.is_empty() {
        "single_class"
    } else if knee_cut.has_knee {
        "knee"
    } else {
        "no_knee"
    };
    let p = qual_p(1.0, k, alternative);
    ClusterQual {
        knee: knee_cut.has_knee.then_some(cut),
        cut,
        shape,
        cluster_size: 1 + band,
        k,
        excluded,
        alternative,
        p,
    }
}

/// The canonical diploid pair mixture: P(read | generated from the pair)
/// with the equal 1/2 homolog prior (the balanced 15x/15x sample).
#[inline]
fn mix_logsumexp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let m = a.max(b);
    m + (0.5 * (a - m).exp() + 0.5 * (b - m).exp()).ln()
}

/// Inverse of flat = j*(j+1)/2 + i (i <= j).
fn unrank_pair(flat: usize) -> (usize, usize) {
    let mut j = (((8.0 * flat as f64 + 1.0).sqrt() - 1.0) / 2.0).floor() as usize;
    while (j + 1) * (j + 2) / 2 <= flat {
        j += 1;
    }
    while j * (j + 1) / 2 > flat {
        j -= 1;
    }
    let i = flat - j * (j + 1) / 2;
    (i, j)
}

// ---------------------------------------------------------------------------
// The scoring core (the new likelihood layer).
// ---------------------------------------------------------------------------

/// The uniform substitution scoring constants derived from the reads'
/// own quality bytes (epsilon and phred are carried for the receipts).
#[derive(Clone, Copy)]
#[allow(dead_code)]
struct Scoring {
    a: f64,
    b: f64,
    epsilon: f64,
    phred: u8,
}

impl Scoring {
    /// The canonical per-placement closed form: m matching bases and
    /// c mismatching bases. Both the factorized path and the direct
    /// per-base walk produce the integers (m, c); this one expression
    /// turns them into the placement score, so integer equality
    /// implies bit-exact score equality.
    #[inline]
    fn placement(&self, m: u32, c: u32) -> f64 {
        m as f64 * self.a + c as f64 * self.b
    }

    /// The model's derived per-placement MINIMUM likelihood: every
    /// base of the read mismatching (read_len() * ln(eps/3)) — the
    /// rigorous lower bound of the read's likelihood under ANY
    /// placement. It remains the stated lower bound of the local
    /// branch's per-placement terms; as a logsumexp term it
    /// underflows to +0.0 exactly wherever an anchored placement
    /// exists (150 * ln(eps/3) ~ -1546 against placement scores ~ -1).
    /// (The read length is PER-CALL state now — passed by the caller,
    /// never a process-global, so concurrent calls in one process
    /// cannot read each other's length.)
    #[inline]
    fn floor(&self, read_len: usize) -> f64 {
        read_len as f64 * self.b
    }
}

/// The plain two-term logsumexp.
#[inline]
fn logsumexp2(a: f64, b: f64) -> f64 {
    let max = a.max(b);
    if max == f64::NEG_INFINITY {
        return max;
    }
    max + ((a - max).exp() + (b - max).exp()).ln()
}

/// THE ELSEWHERE BRANCH E(read) (slice C, derived — no tuning
/// constants): each read was generated at ONE origin, uniform
/// (max-entropy) over the donor's diploid genome; a candidate models
/// the locale only, so the read's likelihood under a candidate
/// marginalizes generated-here (the candidate's local spell) against
/// generated-elsewhere (the genome outside this locale).
///
/// The elsewhere branch: the read matches its origin at the MEASURED
/// per-base rate (the census measured zero unmatched placement bp
/// over 1,063,522 verified occurrences — reads are exact segments of
/// their origins, at the reads' own Phred-derived rate A), and the
/// origin is uniform over the genome's 2*(G-149) both-strand
/// full-read placements — the same max-entropy spirit as the old
/// Poisson model's beta = M/|U|: the unit of match likelihood spread
/// uniformly over the universe of placements, the universe measured
/// from the panel (the instrument's genome universe, truth-free: the
/// mean per-strain total path length, doubled for the diploid).
///
///   E(read) = len(read) * A - ln(2 * (G_diploid - len(read) + 1))
///
/// With uniform L150 reads and uniform quality bytes this is ONE
/// derived constant per run. THE FLOOR for a read a candidate cannot
/// place is E — not the all-mismatch 150*B (which survives only as
/// the per-placement lower bound inside the local branch, a bit-exact
/// no-op wherever an anchor pins).
#[inline]
fn elsewhere_log_prob(read_len: usize, match_log_prob: f64, genome_diploid_bp: f64) -> f64 {
    read_len as f64 * match_log_prob
        - (2.0 * (genome_diploid_bp - read_len as f64 + 1.0)).ln()
}

/// One anchored placement of a record's pattern onto one candidate
/// fold: the census orientation that discovered it and the
/// FOLD-RELATIVE pinned window starts per canonical anchor (None =
/// absent from the fold's walk or ambiguous under the monotone rule).
struct Pin {
    orientation: u8,
    pinned: Vec<Option<i64>>,
    /// The merged bp of the pinned anchor windows in hull coordinates
    /// (mirror-invariant: the mirrored window set is the reflection of
    /// the forward set, so the merged coverage length is identical).
    backbone_bp: u64,
}

/// The spelled fold sequence as a JSON-safe string (the sequence
/// emission's material: the fold's own spelled walk from the local
/// graph — the panel fetch of the member row [start, end) window,
/// the same convention-(a) material the truth-QV alignments use).
/// Fail-closed on non-ASCII bytes (the panel is DNA ASCII).
pub fn fold_sequence_string(seq: &[u8]) -> io::Result<String> {
    if seq.iter().any(|b| !b.is_ascii()) {
        return Err(invalid(
            "a called fold's spelled sequence carries non-ASCII bytes",
        ));
    }
    Ok(seq.iter().map(|&b| b as char).collect())
}

/// A minimal fnmatch-style glob matcher (`*` any run, `?` one char;
/// no escapes) for the call-locus --exclude-paths operator (the
/// held-out masking experiment's stranger filter, folded onto the
/// call path: candidate rows whose member path names match any glob
/// leave the candidate universe, loudly counted).
pub fn path_matches_glob(path_name: &str, pattern: &str) -> bool {
    fn inner(name: &[u8], pattern: &[u8]) -> bool {
        match pattern.first() {
            None => name.is_empty(),
            Some(b'*') => (1..=name.len())
                .rev()
                .any(|k| inner(&name[k..], &pattern[1..]))
                || inner(name, &pattern[1..]),
            Some(b'?') => !name.is_empty() && inner(&name[1..], &pattern[1..]),
            Some(&b) => name.first() == Some(&b) && inner(&name[1..], &pattern[1..]),
        }
    }
    inner(path_name.as_bytes(), pattern.as_bytes())
}

/// The component's axis rows (start, end), in axis order, from the
/// locality maps (windowN.gfa.map.json) - the same rows the scoring
/// path derives from --partition-graphs (the axis intervals of the
/// component lane, sorted by start). One spelling of the axis for
/// every region-resolution site.
pub fn component_axis_rows(
    partition_graphs: &std::path::Path,
    component: &str,
) -> io::Result<Vec<(u64, u64)>> {
    let mut rows: Vec<(u64, u64, u32)> = Vec::new();
    for entry in std::fs::read_dir(partition_graphs)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.ends_with(".gfa.map.json") {
            continue;
        }
        let text = std::fs::read_to_string(&path)?;
        let map: PartitionMap = serde_json::from_str(&text)?;
        for axis in map.members.iter().filter(|m| m.path_name == component) {
            rows.push((axis.start, axis.end, map.partition));
        }
    }
    ensure(
        !rows.is_empty(),
        &format!(
            "the component {component:?} carries no axis rows in {:?}",
            partition_graphs
        ),
    )?;
    rows.sort_by_key(|&(start, end, partition)| (start, end, partition));
    Ok(rows.into_iter().map(|(start, end, _)| (start, end)).collect())
}

/// THE REGION-AS-INTERFACE RESOLVER (the owner's ruling: the CLI's
/// interface is COORDINATES - "detect any locus, any position, any
/// allele"): a region `chr:start-end` (0-based half-open, the panel
/// path-naming frame) resolves against the component's axis rows
/// (the coordinate frame) to the covered windows' locus ids - the
/// component-local indices of the rows that overlap the region in
/// the component's own path frame. The returned loci are in axis
/// order, distinct and sorted by construction - the --loci
/// convention verbatim.
pub fn resolve_region_loci(
    component: &str,
    axis_rows: &[(u64, u64)],
    region: &str,
) -> io::Result<Vec<u32>> {
    let invalid_region = || invalid(&format!(
        "the region must be chr:start-end (0-based half-open, the \
         panel path-naming frame): {region:?}"));
    let (chrom, bounds) = region.split_once(':').ok_or_else(invalid_region)?;
    let (lo, hi) = bounds.split_once('-').ok_or_else(invalid_region)?;
    let (lo, hi): (u64, u64) = (lo.parse().map_err(|_| invalid_region())?,
        hi.parse().map_err(|_| invalid_region())?);
    ensure(lo < hi, "the region is empty: {region:?}")?;
    // the component's own path frame: its name ends chr:start-end
    let (name, interval) = component
        .rsplit_once(':')
        .ok_or_else(|| invalid(&format!(
            "the component name carries no coordinate frame: {component:?}")))?;
    let component_chrom = name
        .rsplit_once('#')
        .map(|(_, chrom)| chrom)
        .ok_or_else(|| invalid(&format!(
            "the component name carries no chromosome: {component:?}")))?;
    ensure(component_chrom == chrom,
        "the region's chromosome {chrom:?} is not the axis frame {component_chrom:?}")?;
    let (offset, end): (u64, u64) = interval
        .split_once('-')
        .and_then(|(s, e)| Some((s.parse().ok()?, e.parse().ok()?)))
        .ok_or_else(|| invalid(&format!(
            "the component name carries no interval: {component:?}")))?;
    ensure(lo >= offset && hi <= end,
        "the region {region:?} is outside the axis frame {component:?}")?;
    let (lo, hi) = (lo - offset, hi - offset);
    let mut loci = Vec::new();
    for (index, (start, end)) in axis_rows.iter().enumerate() {
        if *start < hi && *end > lo {
            loci.push(index as u32);
        }
    }
    ensure(!loci.is_empty(),
        "the region {region:?} covers no axis window (frame {component:?})")?;
    Ok(loci)
}

/// True when the member row leaves the candidate universe under the
/// exclusion globs. THE AXIS ROW IS NEVER EXCLUDED (it is the
/// coordinate frame of the call, not a candidate; the full substrate
/// masking of the held-out experiment - axis included, panel rebuilt -
/// remains the experiment's own path).
pub fn row_is_excluded(path_name: &str, component: &str, globs: &[String]) -> bool {
    path_name != component && globs.iter().any(|g| path_matches_glob(path_name, g))
}

/// THE FLANK-OVERLAP CANDIDATE FILTER (the owner's design, the
/// presence semantics ruled at the fork): a panel member row is KEPT
/// in a window's candidate universe only if its row CARRIES ALIGNED
/// MATERIAL OVER BOTH FLANKS of the window's axis extent - the first
/// N bp and the last N bp, in the FOCAL REFERENCE's own coordinate
/// frame (the partition maps' member rows live in exactly that
/// frame: the axis row spans the window's extent and every member
/// row's [start, end) is its aligned extent in the same frame). THE
/// HULL-EDGE GEOMETRY, measured before the semantics was ruled: the
/// member rows are alignment-correspondence hulls and structurally
/// stop 0-50bp SHORT of the region's edges (a syncmer window cannot
/// start in the last k+w-1 bases - the truth haplotypes' own rows
/// end ~25bp short), so a literal containment reading of "covers the
/// flank" removes 91-100% of every window's universe INCLUDING the
/// truth folds at ~every window; the ruled semantics is PRESENCE
/// (the owner's own degenerate-case wording): the row must OVERLAP
/// both flank intervals. N is an absolute bp count
/// (`--flank-overlap 10000`) or a percent of the window length
/// (`--flank-overlap 10%` -> ceil(pct x length)). THE DEGENERATE
/// CASE: a window shorter than N reduces the rule to simple
/// presence in the window. OFF BY DEFAULT (unset: the committed
/// candidate universe verbatim).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FlankOverlap {
    Bp(u64),
    Percent(f64),
}

impl std::str::FromStr for FlankOverlap {
    type Err = String;
    fn from_str(text: &str) -> Result<Self, String> {
        let invalid = || format!(
            "invalid --flank-overlap {text:?}: N in bp (e.g. 10000) or a percent (e.g. 10%)"
        );
        if let Some(pct) = text.strip_suffix('%') {
            let pct: f64 = pct.parse().map_err(|_| invalid())?;
            if !(pct > 0.0 && pct <= 100.0) {
                return Err(format!(
                    "invalid --flank-overlap {text:?}: the percent must be in (0, 100]"
                ));
            }
            Ok(FlankOverlap::Percent(pct))
        } else {
            let bp: u64 = text.parse().map_err(|_| invalid())?;
            if bp == 0 {
                return Err(format!(
                    "invalid --flank-overlap {text:?}: N must be positive (0 is the default-off)"
                ));
            }
            Ok(FlankOverlap::Bp(bp))
        }
    }
}

impl FlankOverlap {
    /// The flank size N in bp for a window of the given length.
    pub fn flank_bp(&self, window_len: u64) -> u64 {
        match *self {
            FlankOverlap::Bp(n) => n,
            FlankOverlap::Percent(pct) => ((pct / 100.0) * window_len as f64).ceil() as u64,
        }
    }

    /// The CLI spelling of the setting (the receipt field).
    pub fn setting(&self) -> String {
        match *self {
            FlankOverlap::Bp(n) => n.to_string(),
            FlankOverlap::Percent(pct) => format!("{pct}%"),
        }
    }
}

/// True when the member row [member_start, member_end) stays in the
/// window's candidate universe under the flank-overlap rule (THE
/// PRESENCE SEMANTICS, ruled at the measured fork): the row must
/// CARRY ALIGNED MATERIAL OVER BOTH FLANK INTERVALS of the window's
/// axis extent [axis_start, axis_end) - it must OVERLAP the first N
/// bp [axis_start, axis_start+N) AND the last N bp
/// [axis_end-N, axis_end). A window shorter than N degenerates to
/// simple presence in the window. THE AXIS ROW IS NEVER EXCLUDED (it
/// is the coordinate frame of the call - the --exclude-paths
/// precedent); the caller enforces that, this predicate is pure.
pub fn row_covers_flanks(
    member_start: u64,
    member_end: u64,
    axis_start: u64,
    axis_end: u64,
    flank: FlankOverlap,
) -> bool {
    let window_len = axis_end.saturating_sub(axis_start);
    let n = flank.flank_bp(window_len);
    if window_len < n {
        // THE DEGENERATE CASE (the owner's rule): the window is
        // shorter than N - simple presence in the region.
        return member_start < axis_end && member_end > axis_start;
    }
    // PRESENCE in the first N bp [axis_start, axis_start + n) and
    // the last N bp [axis_end - n, axis_end): the row's aligned
    // hull must reach into both flanks (overlap, not contain - the
    // hull-edge geometry: rows structally stop short of the edges).
    member_start < axis_start + n
        && member_end > axis_start
        && member_start < axis_end
        && member_end > axis_end - n
}

/// The COALESCED CALLED SET (the owner's ambiguity ruling: a tied
/// class call emits every called fold's sequence, not just the
/// representative pair): the distinct folds of all called classes,
/// sorted by fold index. An untied call's set is exactly the
/// representative pair.
pub fn coalesced_called_set(called: &[(usize, usize)]) -> Vec<usize> {
    let mut folds: Vec<usize> = called
        .iter()
        .flat_map(|&(first, second)| [first, second])
        .collect();
    folds.sort_unstable();
    folds.dedup();
    folds
}

/// One folded candidate: identical sequence AND walk (relative) —
/// coalesced copies are ONE hypothesis.
struct Fold {
    seq: Vec<u8>,
    len: u64,
    /// The representative row's CONTAINED walk steps (windows fully
    /// inside the row extent), relative to its start — the pinnable
    /// anchor skeleton; the committed partition GFAs' P lines spell
    /// these steps among their edge-overlapping and gap steps.
    walk: Vec<(u64, i32)>,
    /// node abs id -> relative (bp, sign) occurrences.
    node_positions: HashMap<u32, Vec<(u64, i32)>>,
    members: Vec<MemberRow>,
    /// Contained-step node usage multiset (the QUAL spectrum's
    /// material signature), sorted by key.
    nodes: Vec<(u64, u32)>,
    edges: Vec<(u64, u32)>,
}

impl Fold {
    fn positions_of(&self, node_abs: u32) -> Vec<(u64, i32)> {
        self.node_positions
            .get(&node_abs)
            .cloned()
            .unwrap_or_default()
    }
}

/// One evidence unit: a distinct read variant of one record with its
/// instance count (the record-once discipline: the record's instances
/// are its evidence, each instance voting once per row).
struct Unit {
    record: usize,
    variant: usize,
    mirror: bool,
    read: usize,
    w_lo: u64,
    count: u64,
}

/// The bracketing pinned anchor for one hull-relative skipped
/// position d (slice A's committed serving rule: the pinned anchor
/// with the greatest own position <= d, else the least).
fn serving_anchor(pinned: &[Option<i64>], own_hull: &[u64], d: i64) -> Option<(usize, i64)> {
    let mut best: Option<(usize, i64)> = None;
    for (index, &r) in pinned.iter().enumerate() {
        let Some(_r) = r else { continue };
        let own_pos = own_hull[index] as i64;
        let better = match best {
            None => true,
            Some((_, best_pos)) => {
                let (own_left, best_left) = (own_pos <= d, best_pos <= d);
                if own_left && !best_left {
                    true
                } else if own_left == best_left {
                    if own_left {
                        own_pos > best_pos
                    } else {
                        own_pos < best_pos
                    }
                } else {
                    false
                }
            }
        };
        if better {
            best = Some((index, own_pos));
        }
    }
    best
}

/// The FACTORIZED per-placement score for one read variant: the
/// backbone (merged pinned-window bp, matches by node identity and
/// sign — asserted per base by the direct check) plus the skipped-base
/// votes against the fold's spelled sequence. Returns None when the
/// placement is invalid for this variant (a voted base projects
/// outside the comparison domain — the fold's extent, or under the
/// normalized convention the fold's TERRITORY IMAGE [img_lo, img_hi),
/// the uniform placement prior's support requires the read fully
/// inside it).
#[allow(clippy::too_many_arguments)]
fn placement_score(
    unit: &Unit,
    canonical: &[(i32, u64)],
    pin: &Pin,
    fold_seq: &[u8],
    _fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
    img_lo: i64,
    img_hi: i64,
) -> Option<(u32, u32)> {
    let read = &reads[unit.read];
    let forward = physical_forward(pin.orientation, unit.mirror);
    let own_positions = own_positions_of(canonical, unit.mirror, unit.w_lo);
    let own_hull = own_positions_of(canonical, unit.mirror, 0);
    let skipped = skipped_positions(&own_positions, k, read.len() as u64);
    let mut matches = pin.backbone_bp as u32;
    let mut mismatches = 0u32;
    for &(lo, hi) in &skipped {
        for i in lo..hi {
            let d = i as i64 - unit.w_lo as i64;
            let (index, own_pos) = serving_anchor(&pin.pinned, &own_hull, d)?;
            let r = pin.pinned[index].unwrap();
            let coord = if forward {
                r + (d - own_pos)
            } else {
                r + (own_pos + k as i64 - 1 - d)
            };
            if coord < img_lo || coord >= img_hi {
                return None;
            }
            let base = fold_seq[coord as usize];
            let read_base = if forward {
                read[i as usize]
            } else {
                complement(read[i as usize])
            };
            if read_base == base {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    Some((matches, mismatches))
}

/// The DIRECT per-base recomputation of one placement's score for one
/// variant (the brute-force side of the exactness proof, computed
/// in-process over the whole pilot domain): the backbone is verified
/// base-by-base against the fold's sequence instead of trusted from
/// node identity, and the skipped votes are recomputed. Integer
/// equality with the factorized (m, c) is the exactness assertion.
#[allow(clippy::too_many_arguments)]
fn direct_score(
    unit: &Unit,
    canonical: &[(i32, u64)],
    pin: &Pin,
    fold_seq: &[u8],
    fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
) -> io::Result<(u32, u32)> {
    let read = &reads[unit.read];
    let forward = physical_forward(pin.orientation, unit.mirror);
    let own_positions = own_positions_of(canonical, unit.mirror, unit.w_lo);
    let own_hull = own_positions_of(canonical, unit.mirror, 0);
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    // The backbone, per base over the MERGED pinned-window coverage:
    // syncmer windows overlap heavily (k=63, spacing <= w), so the
    // per-anchor loops would count the overlap bases several times.
    // Each base is verified through ONE covering pinned anchor (the
    // greatest own position at or before the base whose window covers
    // it) against the fold's sequence by direct comparison.
    let backbone_windows: Vec<(u64, u64)> = pin
        .pinned
        .iter()
        .enumerate()
        .filter_map(|(j, &r)| r.map(|_| (own_hull[j], own_hull[j] + k)))
        .collect();
    for (lo, hi) in merge_windows(backbone_windows) {
        for d in lo..hi {
            let mut cover: Option<usize> = None;
            for (j, &r) in pin.pinned.iter().enumerate() {
                if r.is_none() {
                    continue;
                }
                let p = own_hull[j];
                if p <= d && d < p + k {
                    match cover {
                        None => cover = Some(j),
                        Some(current) if p > own_hull[current] => cover = Some(j),
                        _ => {}
                    }
                }
            }
            let j = cover.expect("merged coverage base without a covering anchor");
            let r = pin.pinned[j].unwrap();
            let t = d - own_hull[j];
            let read_i = unit.w_lo as usize + d as usize;
            let read_base = if forward {
                read[read_i]
            } else {
                complement(read[read_i])
            };
            let coord = if forward {
                r + t as i64
            } else {
                r + k as i64 - 1 - t as i64
            };
            ensure(
                coord >= 0 && (coord as usize) < fold_seq.len(),
                "direct backbone outside the fold",
            )?;
            if read_base == fold_seq[coord as usize] {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    // The skipped votes, per base.
    let skipped = skipped_positions(&own_positions, k, read.len() as u64);
    for &(lo, hi) in &skipped {
        for i in lo..hi {
            let d = i as i64 - unit.w_lo as i64;
            let (index, own_pos) = serving_anchor(&pin.pinned, &own_hull, d)
                .ok_or_else(|| invalid("direct vote without a serving anchor"))?;
            let r = pin.pinned[index].unwrap();
            let coord = if forward {
                r + (d - own_pos)
            } else {
                r + (own_pos + k as i64 - 1 - d)
            };
            ensure(
                coord >= 0 && coord < fold_len as i64,
                "direct vote outside the fold",
            )?;
            let base = fold_seq[coord as usize];
            let read_base = if forward {
                read[i as usize]
            } else {
                complement(read[i as usize])
            };
            if read_base == base {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    let _ = own_positions;
    Ok((matches, mismatches))
}

/// The per-unit per-fold log-likelihood: the uniform placement prior
/// (both strands, 2 * (len - 149) full-read placements) plus the
/// logsumexp over the valid anchored placements. NEG_INFINITY when no
/// placement is valid (the fold shares no anchor, or every pin leaves
/// the read partially outside the fold — the anchored evidence model
/// has no placement there; the exactness receipt measures what the
/// unanchored floor would be on the sample).
fn unit_log_likelihood(
    scores: &[f64],
    fold_len: u64,
    read_len: usize,
) -> f64 {
    if scores.is_empty() || fold_len < read_len as u64 {
        return f64::NEG_INFINITY;
    }
    let prior = -(2.0 * (fold_len - read_len as u64 + 1) as f64).ln();
    let best = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = scores.iter().map(|&s| (s - best).exp()).sum();
    prior + best + sum.ln()
}

// ---------------------------------------------------------------------------
// THE ROW-COMPRESSED SCORING MATRIX (candidate C — the marathon
// profile's dense-matrix memory tail, the guard-kill class's fix).
//
// THE MEASURED DISEASE (SS10.6.3's receipts): the per-locus scoring
// matrix holds a dense (folds x units) ll/pinned/mismatches triple at
// 13 B/pair — ~33GB at the worst measured marathon locus (T03411:
// 490 x 5,160,612 = 2.53e9 pairs) — the 64GiB guard-kill class's
// direct cause (134 kills; completions 54-71GB). At the >=98%-no-pin
// loci the mass is the CONSTANT no-pin triple (ll = the elsewhere
// branch, pinned = false, mismatches = 0) for nearly every pair.
//
// THE EXACTNESS ARGUMENT (bit-for-bit, NO approximation): the pair
// loop's pinned path can never produce the no-pin triple — a pinned
// pair's local branch is finite (scores nonempty), so its
// marginalization logsumexp2(local, E) is STRICTLY greater than E;
// and an unpinned pair's mismatches stay 0 (nothing accumulated) with
// ll exactly E (the empty-scores local branch is NEG_INFINITY and
// logsumexp2(-inf, E) = E). The triple therefore PARTITIONS: an entry
// equals the no-pin constant iff the unit is unpinned on that fold —
// a value-identical compression. Every consumer reads the SAME f64
// bits, the SAME bool, the SAME u32 it read from the dense form, in
// the SAME iteration order (the class scan's accumulation order is
// unit order, unchanged), so every downstream float operation is
// unchanged: the outputs are byte-identical by construction.
//
// THE ROW ENCODING (a DERIVED choice per row, never tuned): a row
// keeps the committed dense form iff the sparse form would be no
// smaller — the sparse form costs 16 B/pinned entry (unit + ll +
// mismatches) against the dense 13 B/unit — so the pinned-majority
// marathon loci (T03641 at 84% pinned) keep the committed size
// EXACTLY (never larger), and the no-pin-majority loci (T03803 at 2%
// pinned) collapse ~39x. The bound is per row: min(13 x n_units,
// 16 x pinned_in_row) bytes.
// ---------------------------------------------------------------------------

/// The SPARSE row: the pinned entries only (unit, ll, mismatches),
/// sorted by unit, structure-of-arrays. An absent unit reads the
/// no-pin constant (elsewhere, false, 0) — exactly the dense form's
/// value there, per the exactness argument above. (A pinned=false
/// entry always equals the constant — proven above — so the sparse
/// row IS the pinned set.)
struct ScoreRowSparse {
    units: Vec<u32>,
    lls: Vec<f64>,
    mismatches: Vec<u32>,
}

/// The DENSE row: the committed layout verbatim (13 B/unit).
struct ScoreRowDense {
    ll: Vec<f64>,
    pinned: Vec<bool>,
    mismatches: Vec<u32>,
}

enum ScoreRow {
    Dense(ScoreRowDense),
    Sparse(ScoreRowSparse),
}

impl ScoreRow {
    /// Build the row from the pair loop's dense accumulators, taking
    /// the SMALLER of the two encodings (the derived bound, never a
    /// tuned crossover): dense bytes are 13/unit, sparse bytes are
    /// 16 per entry that differs from the no-pin constant.
    fn from_pair_loop(
        ll: Vec<f64>,
        pinned: Vec<bool>,
        mismatches: Vec<u32>,
        elsewhere: f64,
    ) -> ScoreRow {
        let default_bits = elsewhere.to_bits();
        let sparse_len = ll
            .iter()
            .zip(&pinned)
            .zip(&mismatches)
            .filter(|((&value, &is_pinned), &m)| {
                is_pinned || value.to_bits() != default_bits || m != 0
            })
            .count();
        if 16 * sparse_len as u64 >= 13 * ll.len() as u64 {
            return ScoreRow::Dense(ScoreRowDense {
                ll,
                pinned,
                mismatches,
            });
        }
        let mut units = Vec::with_capacity(sparse_len);
        let mut lls = Vec::with_capacity(sparse_len);
        let mut mism = Vec::with_capacity(sparse_len);
        for (unit, ((&value, &is_pinned), &m)) in
            ll.iter().zip(&pinned).zip(&mismatches).enumerate()
        {
            if is_pinned || value.to_bits() != default_bits || m != 0 {
                units.push(unit as u32);
                lls.push(value);
                mism.push(m);
            }
        }
        ScoreRow::Sparse(ScoreRowSparse {
            units,
            lls,
            mismatches: mism,
        })
    }

    /// ll at `unit` (a random-access probe; the monotone-scan form is
    /// the cursor below). Absent sparse entries read the no-pin
    /// constant.
    fn ll_probe(&self, unit: usize, elsewhere: f64) -> f64 {
        match self {
            ScoreRow::Dense(dense) => dense.ll[unit],
            ScoreRow::Sparse(sparse) => match sparse
                .units
                .binary_search(&(unit as u32))
            {
                Ok(index) => sparse.lls[index],
                Err(_) => elsewhere,
            },
        }
    }

    /// The pinned flag at `unit`.
    fn pinned_probe(&self, unit: usize) -> bool {
        match self {
            ScoreRow::Dense(dense) => dense.pinned[unit],
            ScoreRow::Sparse(sparse) => {
                sparse.units.binary_search(&(unit as u32)).is_ok()
            }
        }
    }

    /// The mismatches at `unit` (the committed u32 value).
    fn mismatches_probe(&self, unit: usize) -> u32 {
        match self {
            ScoreRow::Dense(dense) => dense.mismatches[unit],
            ScoreRow::Sparse(sparse) => {
                match sparse.units.binary_search(&(unit as u32)) {
                    Ok(index) => sparse.mismatches[index],
                    Err(_) => 0,
                }
            }
        }
    }

    /// The pinned pairs of this row in UNIT ORDER: (unit,
    /// mismatches) — the enumeration the exactness sample and the
    /// named-fold strata walk (the dense form's `(0..n_units)
    /// .filter(pinned)` in the same order).
    fn pinned_units(&self) -> impl Iterator<Item = (usize, u32)> + '_ {
        match self {
            ScoreRow::Dense(dense) => dense
                .pinned
                .iter()
                .enumerate()
                .filter(|(_, &is_pinned)| is_pinned)
                .map(|(unit, _)| (unit, dense.mismatches[unit]))
                .collect::<Vec<_>>()
                .into_iter(),
            ScoreRow::Sparse(sparse) => sparse
                .units
                .iter()
                .zip(&sparse.mismatches)
                .map(|(&unit, &m)| (unit as usize, m))
                .collect::<Vec<_>>()
                .into_iter(),
        }
    }

    /// Bit-exact equality of the two rows' LL COLUMNS (the merge's
    /// own safety predicate — the committed comparison was over the
    /// ll matrix alone, and this one is too). Both rows walk units in
    /// order, so every compared pair of values is the pair the dense
    /// comparison compared.
    fn ll_bits_equal(&self, other: &ScoreRow, n_units: usize, elsewhere: f64) -> bool {
        let mut a = RowCursor::new(self);
        let mut b = RowCursor::new(other);
        for unit in 0..n_units {
            if a.ll_at(unit, elsewhere).to_bits()
                != b.ll_at(unit, elsewhere).to_bits()
            {
                return false;
            }
        }
        true
    }

    /// The canonical LL-column hash (the merge's bucketing): the
    /// value bits in unit order — variant-independent, so two equal
    /// rows hash equal whichever encoding each took. (The hash only
    /// buckets candidates; the exact predicate above decides.)
    fn ll_bits_hash(&self, n_units: usize, elsewhere: f64) -> u64 {
        let mut cursor = RowCursor::new(self);
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        for unit in 0..n_units {
            for byte in cursor.ll_at(unit, elsewhere).to_bits().to_le_bytes() {
                hash ^= byte as u64;
                hash = hash.wrapping_mul(0x100_0000_01b3);
            }
        }
        hash
    }

    /// Serialize the row's LL column in the DENSE form's exact bytes
    /// (serde_json's `Vec<Vec<f64>>` element shape: `[` + the values
    /// joined by `,` + `]`), reconstructing the dense values from the
    /// compressed form one at a time — the caller streams this row to
    /// the ingredients spill, so the dense column never materializes.
    fn write_ll_row_json(
        &self,
        writer: &mut impl std::io::Write,
        n_units: usize,
        elsewhere: f64,
    ) -> std::io::Result<()> {
        let mut cursor = RowCursor::new(self);
        writer.write_all(b"[")?;
        for unit in 0..n_units {
            if unit > 0 {
                writer.write_all(b",")?;
            }
            serde_json::to_writer(&mut *writer, &cursor.ll_at(unit, elsewhere))?;
        }
        writer.write_all(b"]")?;
        Ok(())
    }

    /// The row's live bytes (the memory bound's own accounting).
    fn live_bytes(&self) -> u64 {
        match self {
            ScoreRow::Dense(dense) => {
                13 * dense.ll.len() as u64
            }
            ScoreRow::Sparse(sparse) => {
                16 * sparse.units.len() as u64
            }
        }
    }
}

/// A monotone cursor over one row (the class scan's per-unit walk):
/// O(1) amortized per unit, sparse or dense — no binary search on the
/// scan path.
struct RowCursor<'a> {
    row: &'a ScoreRow,
    index: usize,
}
impl<'a> RowCursor<'a> {
    fn new(row: &'a ScoreRow) -> RowCursor<'a> {
        RowCursor { row, index: 0 }
    }
    /// The ll at `unit`, assuming monotone non-decreasing access.
    fn ll_at(&mut self, unit: usize, elsewhere: f64) -> f64 {
        match self.row {
            ScoreRow::Dense(dense) => dense.ll[unit],
            ScoreRow::Sparse(sparse) => {
                while self.index < sparse.units.len()
                    && (sparse.units[self.index] as usize) < unit
                {
                    self.index += 1;
                }
                if self.index < sparse.units.len()
                    && sparse.units[self.index] as usize == unit
                {
                    sparse.lls[self.index]
                } else {
                    elsewhere
                }
            }
        }
    }
}

/// The per-locus scoring matrix: one row per fold, the committed
/// triple (ll, pinned, mismatches) in the per-row smaller encoding.
struct ScoreMatrix {
    n_units: usize,
    elsewhere: f64,
    rows: Vec<ScoreRow>,
    /// THE UNIT-FOLD REDUCTION's plan (the class scan's unit-dimension
    /// twin of the bit-equal-row memo; env IMPG_UNIT_FOLD, DEFAULT
    /// OFF — see `UnitFoldPlan`): `None` IS the committed unit-order
    /// accumulation, verbatim. Attached once per locus by `fold_units`
    /// at the class scan site, never built with the flag unset.
    unit_fold: Option<UnitFoldPlan>,
}

impl ScoreMatrix {
    fn n_folds(&self) -> usize {
        self.rows.len()
    }

    /// ll at (fold, unit) — the dense form's `matrix[fold][unit]`.
    fn ll(&self, fold: usize, unit: usize) -> f64 {
        self.rows[fold].ll_probe(unit, self.elsewhere)
    }

    /// The pinned flag — the dense form's `pinned[fold][unit]`.
    fn pinned(&self, fold: usize, unit: usize) -> bool {
        self.rows[fold].pinned_probe(unit)
    }

    /// The mismatches — the dense form's `mismatches[fold][unit]`.
    fn mismatches(&self, fold: usize, unit: usize) -> u32 {
        self.rows[fold].mismatches_probe(unit)
    }

    /// The pinned pairs of a fold in unit order.
    fn pinned_units(&self, fold: usize) -> impl Iterator<Item = (usize, u32)> + '_ {
        self.rows[fold].pinned_units()
    }

    /// THE CLASS PAIR SUM — the exhaustive classes' own accumulation,
    /// verbatim in order and arithmetic: `total += count *
    /// mix_logsumexp(ll_a[u], ll_b[u])` over units ascending. The
    /// cursors read the same values the dense rows held (the
    /// compressed form's absent entries ARE the no-pin constant), so
    /// the sum is the committed sum bit for bit.
    ///
    /// With a `unit_fold` plan attached (env IMPG_UNIT_FOLD) the walk is
    /// over the plan's SURVIVING COLUMN GROUPS instead: the same
    /// `count * mix_logsumexp(...)` term per group at its smallest unit
    /// index, ascending — the original unit order with the redundant
    /// terms removed — plus the constant columns' ONE scalar
    /// `count_total * c` as the sum's final term (mix(c,c) = c exactly,
    /// so it is pair-independent; it stays IN the value, so every
    /// reported class LL still carries the whole mass). A different
    /// association of the same addends: exact in real arithmetic, NOT
    /// bit-identical (see `UnitFoldPlan`).
    fn class_pair_sum(
        &self,
        first: usize,
        second: usize,
        unit_counts: &[f64],
    ) -> f64 {
        if let Some(plan) = self.unit_fold.as_ref() {
            let mut a = RowCursor::new(&self.rows[first]);
            let mut b = RowCursor::new(&self.rows[second]);
            let mut total = 0.0f64;
            for &(unit, count) in plan.groups() {
                let unit = unit as usize;
                total += count
                    * mix_logsumexp(
                        a.ll_at(unit, self.elsewhere),
                        b.ll_at(unit, self.elsewhere),
                    );
            }
            return match plan.constant_scalar {
                Some(scalar) => total + scalar,
                None => total,
            };
        }
        let mut a = RowCursor::new(&self.rows[first]);
        let mut b = RowCursor::new(&self.rows[second]);
        let mut total = 0.0f64;
        for (unit, &count) in unit_counts.iter().enumerate() {
            total += count
                * mix_logsumexp(
                    a.ll_at(unit, self.elsewhere),
                    b.ll_at(unit, self.elsewhere),
                );
        }
        total
    }

    /// Build and attach the unit-fold plan (THE UNIT-FOLD REDUCTION,
    /// env IMPG_UNIT_FOLD — see `UnitFoldPlan`). `None` when the
    /// reduction declines (one fold: every column is trivially constant
    /// and the whole sum would become a single product — a shape
    /// measured nowhere), leaving the committed walk in place.
    fn fold_units(&mut self, unit_counts: &[f64]) -> Option<&UnitFoldPlan> {
        self.unit_fold = UnitFoldPlan::build(self, unit_counts);
        self.unit_fold.as_ref()
    }

    /// The attached plan (the class-scan site's receipt reader).
    fn unit_fold(&self) -> Option<&UnitFoldPlan> {
        self.unit_fold.as_ref()
    }

    /// The merge's bit-exact column predicate.
    fn ll_bits_equal(&self, a: usize, b: usize) -> bool {
        self.rows[a].ll_bits_equal(&self.rows[b], self.n_units, self.elsewhere)
    }

    /// The merge's bucketing hash.
    fn ll_bits_hash(&self, fold: usize) -> u64 {
        self.rows[fold].ll_bits_hash(self.n_units, self.elsewhere)
    }

    /// The row's canonical COLUMN-CONTENT key hash (the bit-equal-row
    /// memoization's bucketing): the row's pinned entries whose ll
    /// bits differ from the elsewhere constant's, as (unit, bits)
    /// pairs in unit order, length-fed first — variant-independent
    /// and encoding-independent by construction, because the column
    /// IS its non-E entries (an absent entry reads E, and a pinned
    /// entry whose value is exactly E reads E too, so both contribute
    /// nothing to the column; `ll_bits_equal`'s own predicate is
    /// exactly the equality of these entries' sets). Two bit-equal
    /// columns therefore hash equal whichever encoding each row
    /// took; a hash collision is decided by the exact predicate on
    /// match — the merge's own discipline, never the hash alone.
    /// The walk is over the row's PINNED entries only (both encodings
    /// keep them ascending), so the key costs O(pinned), never the
    /// full O(units) column walk.
    fn column_key_hash(&self, fold: usize) -> u64 {
        let elsewhere_bits = self.elsewhere.to_bits();
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |hash: &mut u64, bytes: &[u8]| {
            for &byte in bytes {
                *hash ^= byte as u64;
                *hash = hash.wrapping_mul(0x100_0000_01b3);
            }
        };
        let mut entries: u64 = 0;
        match &self.rows[fold] {
            ScoreRow::Sparse(sparse) => {
                for (&unit, &ll) in sparse.units.iter().zip(&sparse.lls) {
                    let bits = ll.to_bits();
                    if bits != elsewhere_bits {
                        entries += 1;
                        feed(&mut hash, &unit.to_le_bytes());
                        feed(&mut hash, &bits.to_le_bytes());
                    }
                }
            }
            ScoreRow::Dense(dense) => {
                for (unit, &is_pinned) in dense.pinned.iter().enumerate() {
                    if is_pinned {
                        let bits = dense.ll[unit].to_bits();
                        if bits != elsewhere_bits {
                            entries += 1;
                            feed(&mut hash, &(unit as u32).to_le_bytes());
                            feed(&mut hash, &bits.to_le_bytes());
                        }
                    }
                }
            }
        }
        feed(&mut hash, &entries.to_le_bytes());
        hash
    }

    /// The representative's row for a merged fold (the committed
    /// clone, at the row's own encoding).
    fn row_clone(&self, fold: usize) -> ScoreRow {
        match &self.rows[fold] {
            ScoreRow::Dense(dense) => ScoreRow::Dense(ScoreRowDense {
                ll: dense.ll.clone(),
                pinned: dense.pinned.clone(),
                mismatches: dense.mismatches.clone(),
            }),
            ScoreRow::Sparse(sparse) => ScoreRow::Sparse(ScoreRowSparse {
                units: sparse.units.clone(),
                lls: sparse.lls.clone(),
                mismatches: sparse.mismatches.clone(),
            }),
        }
    }

    /// The matrix's live bytes (the memory bound receipt).
    fn live_bytes(&self) -> u64 {
        self.rows.iter().map(|row| row.live_bytes()).sum()
    }

    /// Serialize the WHOLE ll matrix in the dense form's exact bytes
    /// (serde_json's `Vec<Vec<f64>>` shape), one row at a time — the
    /// ingredients receipt's ll_matrix, streamed so the dense
    /// materialization never exists (candidate C's second pole: the
    /// committed form buffered the whole matrix as JSON text in RAM,
    /// ~20 B/pair at the worst marathon loci — a second matrix,
    /// larger than the first, held while the first was still alive).
    fn write_ll_matrix_json(
        &self,
        writer: &mut impl std::io::Write,
    ) -> std::io::Result<()> {
        writer.write_all(b"[")?;
        for fold in 0..self.rows.len() {
            if fold > 0 {
                writer.write_all(b",")?;
            }
            self.rows[fold]
                .write_ll_row_json(&mut *writer, self.n_units, self.elsewhere)?;
        }
        writer.write_all(b"]")?;
        Ok(())
    }

    /// The committed dense construction (the tests' and the merge
    /// tests' own shape).
    fn from_dense(
        matrix: Vec<Vec<f64>>,
        pinned: Vec<Vec<bool>>,
        mismatches: Vec<Vec<u32>>,
        elsewhere: f64,
    ) -> ScoreMatrix {
        let n_units = matrix.first().map(|row| row.len()).unwrap_or(0);
        let rows = matrix
            .into_iter()
            .zip(pinned)
            .zip(mismatches)
            .map(|((ll, pinned), mismatches)| {
                ScoreRow::from_pair_loop(ll, pinned, mismatches, elsewhere)
            })
            .collect();
        ScoreMatrix {
            n_units,
            elsewhere,
            rows,
            unit_fold: None,
        }
    }
}

/// The class scan's accumulation-order toggle (the owner's Ruling
/// 3): Exact = the committed unit-order accumulation (the default,
/// bit-reproducible); SparseSkip = the fast mode (the constant
/// entries skipped, the float accumulation regrouped — deltas
/// measured and named).
#[derive(Clone, Copy)]
enum ClassScanOrder {
    /// The committed exhaustive unit-order accumulation, bit-
    /// reproducible — the receipts' form (the toggle-doctrine restore:
    /// IMPG_CLASS_SCAN_ORDER=exact).
    Exact,
    /// The legacy sparse-skip accumulation order (IMPG_CLASS_SCAN-
    /// ORDER=sparse-skip), unchanged.
    SparseSkip,
    /// THE SANCTIONED BOUNDED CLASS SCAN (the owner's ruling on
    /// SS10.6.21's frontier-1 design; the PRODUCTION DEFAULT, env
    /// unset): the per-fold excess-mass bound E*C + A_f + A_g with
    /// the derived exp-underflow skip threshold — a branch-and-bound
    /// over the class-pair space whose kept pairs are summed by the
    /// EXACT committed unit-order accumulation, so every computed
    /// value is bit-identical to the exhaustive scan's and every
    /// answer (the MAP pair, the called set, the QUAL terms, the
    /// truth fields, the best-outside rival) is preserved BY THE
    /// BOUND'S OWN ARITHMETIC. The skipped pairs' instrument entries
    /// become the labeled bounded receipt (the null skip
    /// certificate + the per-pair bounds), so no bounded output can
    /// be mistaken for an exhaustive one.
    Bounded,
}

/// The ulp distance between two f64s on the monotone bit axis
/// (0 for equal bits) — the sparse-skip compare receipt's named
/// delta unit.
fn ulp_distance(a: f64, b: f64) -> u64 {
    let key = |x: f64| {
        let bits = x.to_bits();
        if bits >> 63 == 0 {
            bits | 0x8000_0000_0000_0000
        } else {
            !bits
        }
    };
    key(a).abs_diff(key(b))
}

/// One row of the fast class-scan index: the PINNED units in
/// ascending order with their ll values (a sparse row's own arrays
/// borrowed verbatim; a dense row's list built once per locus).
enum FastRow<'a> {
    Borrowed { units: &'a [u32], lls: &'a [f64] },
    Owned { units: Vec<u32>, lls: Vec<f64> },
}

impl<'a> FastRow<'a> {
    fn units(&self) -> &[u32] {
        match self {
            FastRow::Borrowed { units, .. } => units,
            FastRow::Owned { units, .. } => units,
        }
    }
    fn lls(&self) -> &[f64] {
        match self {
            FastRow::Borrowed { lls, .. } => lls,
            FastRow::Owned { lls, .. } => lls,
        }
    }
}

/// THE FAST CLASS-SCAN INDEX (the sparse-skip form — the previously
/// REFUSED accumulation-order change, §10.6.8, now the owner's
/// Ruling 3 performance path behind IMPG_CLASS_SCAN_ORDER=1). Per
/// fold: the pinned units ascending with their ll values, plus the
/// locus's total unit-count mass. A class pair's sum then walks only
/// the UNION of the two rows' pinned units — every other (fold, unit)
/// entry is the no-pin constant on BOTH folds (mix_logsumexp(E, E) =
/// E exactly, proven per the §10.6.8 exactness argument), so its
/// whole contribution returns as ONE `elsewhere * constant_mass`
/// term where constant_mass = total_count − (the union's counts).
/// That regrouping CHANGES the float accumulation order relative to
/// the committed unit-order walk — the exact mode's byte-identity
/// invariant — so the fast mode's sums are measured and named
/// (IMPG_CLASS_SCAN_COMPARE), never hidden, and the model string is
/// suffixed so no fast-mode output can be mistaken for an exact one.
struct ClassScanFastIndex<'a> {
    n_units: usize,
    elsewhere: f64,
    total_count: f64,
    rows: Vec<FastRow<'a>>,
}

impl<'a> ClassScanFastIndex<'a> {
    fn build(matrix: &'a ScoreMatrix, unit_counts: &[f64]) -> ClassScanFastIndex<'a> {
        let mut rows = Vec::with_capacity(matrix.rows.len());
        for row in &matrix.rows {
            rows.push(match row {
                ScoreRow::Sparse(sparse) => FastRow::Borrowed {
                    units: &sparse.units,
                    lls: &sparse.lls,
                },
                ScoreRow::Dense(dense) => {
                    let mut units = Vec::new();
                    let mut lls = Vec::new();
                    for (unit, &is_pinned) in dense.pinned.iter().enumerate() {
                        if is_pinned {
                            units.push(unit as u32);
                            lls.push(dense.ll[unit]);
                        }
                    }
                    FastRow::Owned { units, lls }
                }
            });
        }
        ClassScanFastIndex {
            n_units: matrix.n_units,
            elsewhere: matrix.elsewhere,
            total_count: unit_counts.iter().sum(),
            rows,
        }
    }

    /// The sparse-skip class pair sum (the fast mode's accumulation,
    /// deterministic: the union walk ascending, the constant mass as
    /// the final term).
    fn class_pair_sum(&self, first: usize, second: usize, unit_counts: &[f64]) -> f64 {
        let (a_units, a_lls) = (self.rows[first].units(), self.rows[first].lls());
        let (b_units, b_lls) = (self.rows[second].units(), self.rows[second].lls());
        let (mut ai, mut bi) = (0usize, 0usize);
        let mut total = 0.0f64;
        let mut union_count = 0.0f64;
        while ai < a_units.len() || bi < b_units.len() {
            let unit = if bi >= b_units.len()
                || (ai < a_units.len() && a_units[ai] <= b_units[bi])
            {
                a_units[ai]
            } else {
                b_units[bi]
            };
            let la = if ai < a_units.len() && a_units[ai] == unit {
                let value = a_lls[ai];
                ai += 1;
                value
            } else {
                self.elsewhere
            };
            let lb = if bi < b_units.len() && b_units[bi] == unit {
                let value = b_lls[bi];
                bi += 1;
                value
            } else {
                self.elsewhere
            };
            let count = unit_counts[unit as usize];
            total += count * mix_logsumexp(la, lb);
            union_count += count;
        }
        let constant_mass = self.total_count - union_count;
        total + self.elsewhere * constant_mass
    }

    /// The constant (fold, unit) entries this pair's union walk
    /// skips — the compare receipt's honest measurement of what the
    /// sparse-skip form actually skips at this locus.
    fn constant_units_skipped(&self, first: usize, second: usize) -> u64 {
        let (a_units, b_units) = (self.rows[first].units(), self.rows[second].units());
        let (mut ai, mut bi) = (0usize, 0usize);
        let mut union_len = 0usize;
        while ai < a_units.len() || bi < b_units.len() {
            match (ai < a_units.len(), bi < b_units.len()) {
                (true, false) => ai += 1,
                (false, true) => bi += 1,
                (true, true) => {
                    if a_units[ai] < b_units[bi] {
                        ai += 1;
                    } else if b_units[bi] < a_units[ai] {
                        bi += 1;
                    } else {
                        ai += 1;
                        bi += 1;
                    }
                }
                (false, false) => unreachable!(),
            }
            union_len += 1;
        }
        (self.n_units - union_len) as u64
    }
}

/// The f64 exp underflow point, DERIVED BY EVALUATION (never a tuned
/// constant): the largest f64 t with exp(t) == +0.0, found by a
/// bit-axis binary search between a known-nonzero anchor (-744 —
/// e^-744 ~ 2^-1073.4, far above the boundary) and the
/// real-arithmetic boundary -1075*ln 2 (known to underflow). Every
/// x <= t underflows to exactly +0.0, so a class pair whose exact
/// sum sits at or below best + t contributes EXACTLY +0.0 to every
/// QUAL posterior term the exhaustive run computes — the derived
/// skip threshold of the sanctioned bounded class scan (the same
/// derivation the committed test exp_underflow_point_derived_by_
/// evaluation proves, and the probe's own constant).
fn exp_underflow_point() -> f64 {
    let mut lo = -744.0f64;
    let mut hi = -1075.0f64 * 2.0f64.ln();
    assert!(lo.exp() > 0.0);
    assert_eq!(hi.exp(), 0.0);
    while hi.to_bits() - lo.to_bits() > 1 {
        // the midpoint on the bit axis WITHOUT the u64 sum overflow
        // (both anchors are negative doubles: their bit patterns sum
        // past u64::MAX).
        let mid = f64::from_bits(lo.to_bits() + (hi.to_bits() - lo.to_bits()) / 2);
        if mid == lo || mid == hi {
            break;
        }
        if mid.exp() == 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// THE PER-FOLD EXCESS-MASS BOUND INDEX (the sanctioned frontier-1
/// design, SS10.6.21 — the owner's "take the better, faster one"):
/// every per-(fold, unit) log-likelihood is either exactly the
/// elsewhere constant E (the no-pin branch — the committed
/// row-compression exactness argument: a pinned pair's local branch
/// is finite, so its marginalization logsumexp2(local, E) rounds to
/// >= E, with float equality only by rounding; an unpinned entry is
/// exactly E) or E + alpha with alpha >= 0, so with C = sum_u c_u
/// the class pair sum decomposes exactly (real arithmetic) as
/// S(f,g) = E*C + sum_u c_u * mix0(alpha_f[u], alpha_g[u]) and, since
/// mix0(a, b) <= a + b for a, b >= 0 (WLOG a >= b: mix0 = a +
/// ln(0.5 + 0.5e^{b-a}) <= a),
///
/// ```text
///     S(f,g) <= E*C + A_f + A_g,   A_f = sum_u c_u * alpha_f[u]
/// ```
///
/// — ONE SCALAR PER FOLD, one walk over the fold's pinned entries:
/// the whole bound index costs a fraction of ONE class-pair walk per
/// fold pair (the probe's receipted 0.01-1.2% of the scan it
/// bounds). If E is not finite the bounds are NaN and NOTHING is
/// skipped (NaN fails every skip comparison) — the conservative
/// exhaustive fallback, never a silent wrong prune.
struct ClassScanBoundIndex {
    /// E * C (one multiply; C = the locus's total unit-count mass).
    e_times_c: f64,
    /// A_f per fold (sum over the fold's pinned units ascending —
    /// deterministic, one spelling).
    excess: Vec<f64>,
    /// THE DERIVED ROUNDING-ERROR SLACK's own inputs (SS10.6.22's
    /// fix: the bound's f64 evaluation can sit strictly BELOW the
    /// exact unit-order sum — the receipted 12-ulp shortfall at
    /// T03226 locus 16 — so every skip comparison consumes the
    /// bound TOGETHER WITH its own forward-error term, derived in
    /// `pair_slack` below from the two accumulations' op counts and
    /// magnitudes; no tuned constants anywhere).
    /// C itself (integer-valued, exact in f64 below 2^53).
    total_count: f64,
    /// |E| — the mix evaluation's error is bounded by it (the max
    /// term m sits in [E, 0]).
    abs_elsewhere: f64,
    /// max over every pinned entry of (ll - E) >= 0 — the largest
    /// spread the mix's exp/ln terms can carry at this locus.
    max_spread: f64,
    /// n_units — the exact accumulation's own length.
    n_units: usize,
    /// pinned entries per fold (the bound chain's own length).
    pinned_counts: Vec<usize>,
}

/// The per-locus bounded-scan state carried between the scan's two
/// stages: STAGE ONE (the deterministic search for best) at the
/// class-scan site, STAGE TWO (the derived kept set + the truth
/// window) once the truth pair is resolved, and the best-outside
/// fallback afterwards. Every set it carries is a deterministic
/// function of the matrix alone, so two bounded runs are
/// byte-identical.
struct BoundedScanState {
    /// E*C + A_f + A_g per flat pair — the skip certificate's own
    /// numbers, emitted in the labeled bounded receipt.
    pair_bounds: Vec<f64>,
    /// The DERIVED rounding-error slack per flat pair (the bound's
    /// own forward-error term, `ClassScanBoundIndex::pair_slack` —
    /// every skip comparison consumes bound + slack; SS10.6.22's
    /// fix for the bound's inadmissible-by-rounding shapes).
    pair_slack: Vec<f64>,
    /// The exact sums computed so far (None = skipped so far; filled
    /// by stage two's kept set and the best-outside fallback — all
    /// deterministic).
    computed: Vec<Option<f64>>,
    /// Stage one's running best — the exhaustive run's own value,
    /// bit for bit.
    best: f64,
    /// Stage one's own wall (the bound index + the search),
    /// receipted separately from stage two's kept-set pass.
    stage_one_seconds: f64,
    /// THE BIT-EQUAL-ROW MEMOIZATION's per-locus state (the fourth
    /// reduction): the canonical row classes (with the receipt
    /// counters) and the persistent memo table shared by stage two
    /// and the best-outside fallback — None is the direct form
    /// (IMPG_CLASS_SCAN_MEMO=off).
    scan_memo: Option<ClassScanMemo>,
    /// The memo table: one slot per canonical flat (None = not yet
    /// computed). Sized pair_count (the bounded path's own vectors
    /// are that size already).
    memo_values: Vec<Option<f64>>,
}

impl ClassScanBoundIndex {
    fn build(matrix: &ScoreMatrix, unit_counts: &[f64]) -> ClassScanBoundIndex {
        let total_count: f64 = unit_counts.iter().sum();
        let mut excess = Vec::with_capacity(matrix.rows.len());
        let mut pinned_counts = Vec::with_capacity(matrix.rows.len());
        let mut max_spread = 0.0f64;
        for row in &matrix.rows {
            let mut a = 0.0f64;
            let mut pinned = 0usize;
            let mut spread = 0.0f64;
            match row {
                ScoreRow::Sparse(sparse) => {
                    for (&unit, &ll) in sparse.units.iter().zip(&sparse.lls) {
                        a += unit_counts[unit as usize] * (ll - matrix.elsewhere);
                        pinned += 1;
                        // total_cmp propagates a NaN spread upward
                        // (the not-finite-E case) — the conservative
                        // never-skip direction.
                        if ll.total_cmp(&matrix.elsewhere) == std::cmp::Ordering::Greater
                            && (ll - matrix.elsewhere).total_cmp(&spread)
                                == std::cmp::Ordering::Greater
                        {
                            spread = ll - matrix.elsewhere;
                        }
                    }
                }
                ScoreRow::Dense(dense) => {
                    for (unit, &is_pinned) in dense.pinned.iter().enumerate() {
                        if is_pinned {
                            a += unit_counts[unit] * (dense.ll[unit] - matrix.elsewhere);
                            pinned += 1;
                            if (dense.ll[unit] - matrix.elsewhere).total_cmp(&spread)
                                == std::cmp::Ordering::Greater
                            {
                                spread = dense.ll[unit] - matrix.elsewhere;
                            }
                        }
                    }
                }
            }
            if spread.total_cmp(&max_spread) == std::cmp::Ordering::Greater {
                max_spread = spread;
            }
            excess.push(a);
            pinned_counts.push(pinned);
        }
        ClassScanBoundIndex {
            e_times_c: matrix.elsewhere * total_count,
            excess,
            total_count,
            abs_elsewhere: matrix.elsewhere.abs(),
            max_spread,
            n_units: unit_counts.len(),
            pinned_counts,
        }
    }

    /// The admissible upper bound E*C + A_f + A_g (deterministic:
    /// one spelling, the folds' flat indices in order).
    fn pair_bound(&self, first: usize, second: usize) -> f64 {
        self.e_times_c + self.excess[first] + self.excess[second]
    }

    /// THE DERIVED ROUNDING-ERROR SLACK — a rigorous upper bound, in
    /// the standard IEEE-754 round-to-nearest model (unit roundoff
    /// u = 2^-53; gamma(n) = n*u/(1 - n*u); the standard calculus
    /// gamma(a)*(1 + gamma(b)) <= gamma(a + b)), on
    /// (exact_f64_sum - bound_f64) for the pair:
    ///
    /// ```text
    ///     exact  <= real_sum + err_X      (the unit-order accumulation)
    ///     real_sum <= real_bound          (the mix0 <= alpha+beta lemma)
    ///     real_bound <= bound + err_B     (the bound chain)
    ///     =>  exact <= bound + err_X + err_B.
    /// ```
    ///
    /// err_X — the exact accumulation's own drift ABOVE the real
    /// sum: every unit's term is c * mix with mix in [E, 0]
    /// (logsumexp >= max >= E), so -real_sum <= |E|*C; the n_units
    /// adds plus the n_units count-products contribute at most
    /// gamma(2n)*|E*C| (the op headroom is doubled to carry the
    /// |E*C| magnitude's own read off the rounded e_times_c); mix's
    /// OWN evaluation (m + ln(0.5*e^{a-m} + 0.5*e^{b-m}): two
    /// subtractions bounded by the spread, two exps of values in
    /// (0, 1], one inner sum in [0.5, 1], one ln in [-ln 2, 0], one
    /// final add of m with |m| <= |E|) drifts by at most
    /// u*(|E| + 12*spread + 16) per unit, and the count product's
    /// own rounding one more u*|E| per unit — together
    /// u*(2|E| + 12*spread + 16)*C summed over the counts.
    /// err_B — the bound chain's own drift BELOW the real bound:
    /// the e_times_c product, the two excess accumulations at three
    /// roundings per pinned entry (the ll-E subtract, the count
    /// product, the accumulation add), and the two final adds —
    /// the standard mixed-sign bound gamma(ops)*(sum of term
    /// magnitudes), with ops = 3*(k_f + k_g) + 3, the term
    /// magnitudes read off the rounded index values, and the factor
    /// 2 carrying those reads' own rounding.
    /// A non-finite result fails every skip comparison exactly the
    /// way a NaN bound does — the conservative exhaustive fallback,
    /// never a silent wrong prune.
    fn pair_slack(&self, first: usize, second: usize) -> f64 {
        let u = f64::EPSILON / 2.0;
        let gamma = |ops: f64| {
            if ops * u < 1.0 {
                ops * u / (1.0 - ops * u)
            } else {
                f64::INFINITY
            }
        };
        let n = self.n_units as f64;
        let err_x = gamma(4.0 * n + 4.0) * self.e_times_c.abs()
            + u * (2.0 * self.abs_elsewhere + 12.0 * self.max_spread + 16.0)
                * self.total_count;
        let ops_b =
            3.0 * (self.pinned_counts[first] + self.pinned_counts[second]) as f64 + 3.0;
        let err_b =
            gamma(2.0 * ops_b) * (self.e_times_c.abs() + self.excess[first] + self.excess[second]);
        err_x + err_b
    }
}

/// The unit-fold column hash's FNV-1a 64 prime (the class scan's own
/// hashing discipline: the hash only BUCKETS candidates, the exact
/// column-key comparison decides every merge).
const UNIT_FOLD_HASH_PRIME: u64 = 0x0000_0100_0000_01b3;

/// THE UNIT-FOLD REDUCTION (the class scan's unit-dimension twin of the
/// bit-equal-row memo, behind env `IMPG_UNIT_FOLD` — DEFAULT OFF).
/// The class scan's inner term walks UNITS, and at the deep loci most of
/// those units are redundant (docs/syng-gmem-bwt/unit-redundancy.md:
/// T00048 locus 8 = 45,423 units in 842 bit-distinct ll columns, 40,777
/// of them one constant column). Two reductions, both grouped by the
/// EXACT f64 bit pattern of the WHOLE column (`ll[f][u]` over every fold
/// — nothing rounded, nothing hashed-loose):
///
///   1. CONSTANT columns (every fold bit-equal): the pair's term is
///      `count[u] * mix(c,c)`, and mix(c,c) = c BIT FOR BIT (exp(0) = 1
///      and ln(1) = +0.0 are exact, so `m + 0.0 = m`; the identity is
///      asserted by the `mix_logsumexp` unit tests), so such columns
///      contribute `count_total * c` to EVERY pair identically — ONE
///      scalar, computed once per locus and added as each pair sum's
///      final term. It is never dropped from the reported value: the
///      class LL (and so every QUAL term derived from it) still carries
///      the constant mass, which is what keeps the numbers comparable
///      with the committed run's.
///   2. IDENTICAL columns: units whose whole column is bit-identical are
///      ONE term with their summed count mass (integer-exact in f64),
///      evaluated by the SAME `count * mix_logsumexp(...)` arithmetic as
///      today, in the surviving groups' ORIGINAL unit order (each group
///      at its smallest unit index) — the only change is FEWER TERMS.
///
/// THE EXACTNESS STATEMENT (never over-claimed): the reduction is exact
/// in REAL arithmetic — every dropped addend is replaced by a real
/// quantity equal to it, so the sum's real value is unchanged — but it
/// is NOT bit-identical: it removes addends from the middle of a fixed
/// left-to-right IEEE-754 accumulation and merges duplicate products, so
/// round-to-nearest sees a different association. Measured at T00048
/// locus 8 (unit-redundancy.md §6.1, the same V3 form): 65,520 of the
/// 65,535 kept pairs change bits, max 1,362 ulp (1.59e-7 absolute at the
/// sums' ~-7.1e5 magnitude); the argmax pair is unchanged with a 1.5e7x
/// margin and that locus carries no class ties. This is therefore a
/// LABELLED accumulation change, never the exact path's default: with
/// the env unset the committed unit-order walk runs verbatim.
///
/// THE BOUNDED SCAN'S CAVEAT (stated, not hidden): `pair_bounds` is
/// untouched (it is built from per-fold aggregates, and the row totals
/// are unchanged), but its derived slack prices the BOUND's own
/// evaluation error, not the accumulation's reassociation — so a
/// skipped pair could in principle sit a few hundred ulp above the
/// running best. The measured margin is 1.5e7x and the gate is
/// calls.jsonl, which the receipts' decision fields compare.
struct UnitFoldPlan {
    /// the surviving column groups, ASCENDING by their smallest unit
    /// index: (representative unit, the group's summed unit count) —
    /// the original summation order with the redundant terms removed.
    groups: Vec<(u32, f64)>,
    /// the constant columns' pair-independent contribution: the sum
    /// over the constant groups of `count_total * c`, computed ONCE per
    /// locus; `None` when the locus has no constant column.
    constant_scalar: Option<f64>,
    /// the receipt numbers (the `[unit-fold]` line's own fields).
    units: usize,
    distinct_columns: usize,
    constant_columns: usize,
    constant_units: usize,
    pair_terms: usize,
    build_seconds: f64,
}

impl UnitFoldPlan {
    /// The plan for one locus's FINAL (post-merge) matrix — `None`
    /// when the reduction declines: with fewer than two folds every
    /// column is trivially constant and the whole class sum would
    /// collapse to a single product, a shape measured nowhere, so the
    /// committed walk stays.
    fn build(matrix: &ScoreMatrix, unit_counts: &[f64]) -> Option<UnitFoldPlan> {
        let started = Instant::now();
        let n_folds = matrix.n_folds();
        let n_units = matrix.n_units;
        if n_folds < 2 || n_units == 0 {
            return None;
        }
        let elsewhere = matrix.elsewhere;
        let elsewhere_bits = elsewhere.to_bits();
        // THE COLUMN KEY: per unit, its NON-CONSTANT entries (fold, value
        // bits) in fold order. An entry equal to the elsewhere constant
        // reads the constant whichever way it is stored (absent from a
        // sparse row, or a pinned/unpinned entry whose value IS E), so
        // the non-constant entries ARE the column — the row-compression
        // exactness argument, transposed — and two units' columns are
        // bit-identical iff their key lists are. The walk is over the
        // rows' non-constant entries only (604,977 of 19,486,467 cells at
        // T00048 locus 8), never the full F x U matrix.
        let mut keys: Vec<Vec<(u32, u64)>> = vec![Vec::new(); n_units];
        for fold in 0..n_folds {
            match &matrix.rows[fold] {
                ScoreRow::Sparse(sparse) => {
                    for (&unit, &ll) in sparse.units.iter().zip(&sparse.lls) {
                        let bits = ll.to_bits();
                        if bits != elsewhere_bits {
                            keys[unit as usize].push((fold as u32, bits));
                        }
                    }
                }
                ScoreRow::Dense(dense) => {
                    for (unit, &ll) in dense.ll.iter().enumerate() {
                        let bits = ll.to_bits();
                        if bits != elsewhere_bits {
                            keys[unit].push((fold as u32, bits));
                        }
                    }
                }
            }
        }
        // the bucketing hash (the merge's own key discipline: buckets
        // only — the exact key-list comparison below decides)
        let mut buckets: HashMap<u64, Vec<u32>> = HashMap::new();
        for (unit, key) in keys.iter().enumerate() {
            let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
            for &(fold, bits) in key {
                hash = (hash ^ bits).wrapping_mul(UNIT_FOLD_HASH_PRIME);
                hash = (hash ^ fold as u64).wrapping_mul(UNIT_FOLD_HASH_PRIME);
            }
            hash = (hash ^ key.len() as u64).wrapping_mul(UNIT_FOLD_HASH_PRIME);
            buckets.entry(hash).or_default().push(unit as u32);
        }
        // THE COLUMN CLASSES: within a bucket every unit joins the first
        // class whose key list is EQUAL (a hash collision can never
        // merge). Bucket members are pushed in unit order, so each class
        // keeps its SMALLEST unit index as the representative — the
        // merged group's own summation position. The mass is summed in
        // ascending unit order and the counts are integers, so it is
        // exact in f64.
        let mut classes: Vec<(u32, f64, usize)> = Vec::new();
        for bucket in buckets.values() {
            let mut bucket_classes: Vec<(u32, f64, usize)> = Vec::new();
            for &unit in bucket {
                let key = &keys[unit as usize];
                match bucket_classes
                    .iter_mut()
                    .find(|class| keys[class.0 as usize] == *key)
                {
                    Some(class) => {
                        class.1 += unit_counts[unit as usize];
                        class.2 += 1;
                    }
                    None => bucket_classes.push((
                        unit,
                        unit_counts[unit as usize],
                        1,
                    )),
                }
            }
            classes.extend(bucket_classes);
        }
        classes.sort_unstable_by_key(|class| class.0);
        let mut groups: Vec<(u32, f64)> = Vec::new();
        let mut constant_columns = 0usize;
        let mut constant_units = 0usize;
        let mut constant_scalar: Option<f64> = None;
        for &(rep, mass, members) in &classes {
            let key = &keys[rep as usize];
            // a CONSTANT column: every fold holds the same value bits,
            // i.e. the key is empty (all folds read the elsewhere
            // constant) or it names every fold once with one repeated
            // value.
            let constant_value = if key.is_empty() {
                Some(elsewhere)
            } else if key.len() == n_folds
                && key.iter().all(|&(_, bits)| bits == key[0].1)
            {
                Some(f64::from_bits(key[0].1))
            } else {
                None
            };
            match constant_value {
                Some(value) => {
                    constant_columns += 1;
                    constant_units += members;
                    constant_scalar =
                        Some(constant_scalar.unwrap_or(0.0) + mass * value);
                }
                None => groups.push((rep, mass)),
            }
        }
        Some(UnitFoldPlan {
            pair_terms: groups.len() + usize::from(constant_scalar.is_some()),
            groups,
            constant_scalar,
            units: n_units,
            distinct_columns: classes.len(),
            constant_columns,
            constant_units,
            build_seconds: started.elapsed().as_secs_f64(),
        })
    }

    /// The surviving groups (the pair sum's own walk).
    fn groups(&self) -> &[(u32, f64)] {
        &self.groups
    }

    /// The units folded into another unit's column (the receipt's
    /// "merged" mass = units - distinct columns).
    fn merged_units(&self) -> usize {
        self.units - self.distinct_columns
    }

    /// The distinct columns that are NOT constant (each formed by
    /// merging its units, possibly a class of one).
    fn merged_columns(&self) -> usize {
        self.distinct_columns - self.constant_columns
    }
}

/// THE BIT-EQUAL-ROW MEMOIZATION (the fourth expressibility
/// reduction — the plateau lever, the named Frontier-1 residual,
/// now implemented): the class scan's per-(class-pair, unit) work is
/// a PURE FUNCTION of the pair's evidence-row bit patterns — the
/// committed accumulation `total += count * mix_logsumexp(ll_a[u],
/// ll_b[u])` reads, at every unit, the two rows' column values and
/// nothing else — so two folds whose ll columns are bit-equal (the
/// merge's own `ll_bits_equal` predicate, whichever encoding each
/// row took) produce BIT-IDENTICAL class pair sums against every
/// third fold, and their hom pairs are identical too. The scan
/// therefore computes ONE sum per unordered (row-class, row-class)
/// pair and shares it with every bit-equal pair: the memo SHARES
/// VALUES, it never reorders — every pair's own flat slot is
/// filled with the representative pair's value, so class_lls, best,
/// the called set, the QUAL terms and every receipt field are the
/// direct scan's own bits.
/// THE EXACTNESS PROOF, BY CONSTRUCTION: for folds f, g with
/// representatives r = canon[f], s = canon[g] (rows bit-equal by
/// the exact predicate), class_pair_sum(r, s) executes the IDENTICAL
/// f64 op sequence as class_pair_sum(f, g) — at each unit u the
/// cursor reads the same bits (column equality is exactly the
/// `ll_bits_equal` predicate), `count` is per-unit, and the
/// accumulation walks units ascending in both — and IEEE-754
/// round-to-nearest is a deterministic function of the input bits
/// and the op sequence, so the results are equal bit for bit. Proven per-shape by
/// `bit_equal_rows_share_one_class_sum` and per-fixture by
/// `memoized_scan_is_bit_identical_to_direct_on_randomized_rows`.
/// THE KEY DISCIPLINE (the fold-signature doctrine): bucketing by
/// the canonical column-content hash (`column_key_hash` —
/// encoding-independent, collision-decided by the exact predicate
/// `ll_bits_equal` on match, never the hash alone). The build is
/// O(total pinned entries) — the key walks only the row's pinned
/// non-E entries — plus the within-bucket exact compares (buckets
/// are the merge's own: tiny).
/// THE TOGGLE (the measurement restore path, the bounded scan's own
/// doctrine): IMPG_CLASS_SCAN_MEMO=off restores the direct
/// per-pair form verbatim; the default is ON. The memo changes no
/// emitted byte either way (the values are the direct scan's own),
/// so no output label is warranted; the phase census (env
/// IMPG_LOCUS_PHASE_CENSUS) receipts the hit rates per shape.
struct ClassScanMemo {
    /// fold -> the representative fold of its bit-equal row class
    /// (the class's smallest fold index — deterministic; a fold
    /// whose row is unique represents itself).
    canon: Vec<usize>,
    /// the row classes' count (the receipt's "distinct rows").
    distinct_rows: usize,
    /// the distinct canonical PAIR sums (over the row classes):
    /// distinct_rows * (distinct_rows + 1) / 2.
    canonical_pairs: u64,
    /// the build's own wall (the overhead receipt; counted always,
    /// printed under the census).
    build_seconds: f64,
    /// the receipt counters (always counted, printed under the
    /// census): how many canonical sums the scan sites actually
    /// computed, and how many pair lookups were served from the
    /// memo (the hit rate's numerator is the pairs not computed).
    computed_sums: AtomicU64,
    memo_hits: AtomicU64,
}

impl ClassScanMemo {
    fn build(matrix: &ScoreMatrix) -> ClassScanMemo {
        let started = Instant::now();
        let n = matrix.rows.len();
        let mut canon: Vec<usize> = (0..n).collect();
        // the canonical-content buckets (collision-decided below)
        let mut buckets: HashMap<u64, Vec<usize>> = HashMap::new();
        for fold in 0..n {
            buckets
                .entry(matrix.column_key_hash(fold))
                .or_default()
                .push(fold);
        }
        for bucket in buckets.values() {
            if bucket.len() < 2 {
                continue;
            }
            let mut sorted = bucket.clone();
            sorted.sort_unstable();
            // within a bucket, each fold joins the FIRST earlier fold
            // whose column is bit-equal (the exact predicate decides
            // every candidate — a hash collision can never merge).
            // Bit-equality of columns is an equivalence relation
            // (it IS equality of the column function), so the first
            // match's canon is already the class's smallest index.
            for a in 1..sorted.len() {
                let fold = sorted[a];
                for b in (0..a).rev() {
                    let other = sorted[b];
                    if matrix.ll_bits_equal(fold, other) {
                        canon[fold] = canon[other];
                        break;
                    }
                }
            }
        }
        let distinct_rows = canon.iter().copied().collect::<HashSet<_>>().len();
        let canonical_pairs = distinct_rows as u64 * (distinct_rows as u64 + 1) / 2;
        ClassScanMemo {
            canon,
            distinct_rows,
            canonical_pairs,
            build_seconds: started.elapsed().as_secs_f64(),
            computed_sums: AtomicU64::new(0),
            memo_hits: AtomicU64::new(0),
        }
    }

    /// The pair's CANONICAL flat index: the flat index of its
    /// representatives' pair (canon preserves nothing about ORDER
    /// between the two folds, so the min/max spelling of `unrank_pair`
    /// is restated here; the value at that slot is the pair's own).
    fn canonical_flat(&self, first: usize, second: usize) -> usize {
        let (a, b) = (self.canon[first], self.canon[second]);
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        hi * (hi + 1) / 2 + lo
    }

    /// The memo's per-shape receipt numbers (the census line's own
    /// fields; every counter counted always, printed under the env).
    fn note_computed(&self, count: u64) {
        self.computed_sums.fetch_add(count, Ordering::Relaxed);
    }

    fn note_hits(&self, count: u64) {
        self.memo_hits.fetch_add(count, Ordering::Relaxed);
    }
}

/// THE BOUNDED CLASS SCAN's STAGE ONE — the deterministic
/// branch-and-bound search for best (the strict running-best form
/// the receipts name): pairs are walked in bound-DESCENDING order
/// (ties by flat index — one spelling, so the plan and every
/// computed value are functions of the matrix alone, never of the
/// thread schedule), in rounds of geometrically doubling width (a
/// derived schedule — 1, 2, 4, ... pairs per round, no tuned
/// constant); each round's exact sums are computed in parallel by
/// the committed unit-order accumulation (bit-identical to the
/// exhaustive scan's own value for every pair it touches), the
/// running best is the max over everything computed so far, and
/// the search stops when every remaining pair's bound PLUS ITS OWN
/// DERIVED ROUNDING-ERROR SLACK is strictly below the running best
/// — such a pair's exact sum is <= its bound + slack < best, so it
/// can neither win nor tie the emitted argmax. THE WINNER IS NEVER
/// SKIPPED, BY THE SLACK-CORRECTED ARITHMETIC: the winner's bound +
/// slack >= its exact sum = best >= the running best at every
/// point, and the skip requires bound + slack < running best.
/// THE SLACK (SS10.6.22's fix) is the bound's own forward-error
/// term, DERIVED from the two accumulations' op counts and
/// magnitudes — the bound's f64 evaluation can otherwise sit
/// strictly BELOW the exact unit-order sum (the receipted 12-ulp
/// shortfall at T03226 locus 16), and the stop rule would fire
/// legitimately between bit-identical ties. best is the exhaustive
/// best bit for bit; the provisional called set is the computed
/// ties in flat order (the exhaustive form's own order) — the
/// emitted set is REBUILT from the final class_lls after stage
/// two (see `called_set_from_class_lls`).
/// The pinned-union census (process-wide; a receipt instrument, no
/// answer change). Reset points: none needed — the closing line reports
/// the running totals for the run.
static UNION_CENSUS_PAIRS: AtomicU64 = AtomicU64::new(0);
static UNION_CENSUS_DIFF: AtomicU64 = AtomicU64::new(0);
static UNION_CENSUS_MAXULP: AtomicU64 = AtomicU64::new(0);
static UNION_CENSUS_SKIPPED: AtomicU64 = AtomicU64::new(0);
static UNION_CENSUS_RECHECKS: AtomicU64 = AtomicU64::new(0);

fn bounded_stage_one(
    matrix: &ScoreMatrix,
    unit_counts: &[f64],
    memo: Option<&ClassScanMemo>,
) -> (Vec<f64>, Vec<f64>, Vec<Option<f64>>, f64, Vec<Option<f64>>) {
    let pair_count = matrix.rows.len() * (matrix.rows.len() + 1) / 2;
    let index = ClassScanBoundIndex::build(matrix, unit_counts);
    let pair_bounds: Vec<f64> = (0..pair_count)
        .into_par_iter()
        .map(|flat| {
            let (first, second) = unrank_pair(flat);
            index.pair_bound(first, second)
        })
        .collect();
    // the derived rounding-error slack per pair (never tuned: the
    // bound's own forward-error term — see pair_slack)
    let pair_slack: Vec<f64> = (0..pair_count)
        .into_par_iter()
        .map(|flat| {
            let (first, second) = unrank_pair(flat);
            index.pair_slack(first, second)
        })
        .collect();
    // the deterministic search plan: bound descending (f64::total_cmp
    // orders NaN above everything, so a NaN bound — the not-finite-E
    // case — is walked early and never skips anything), ties by flat
    // index ascending.
    let mut order: Vec<u32> = (0..pair_count as u32).collect();
    order.sort_by(|&a, &b| {
        pair_bounds[b as usize]
            .total_cmp(&pair_bounds[a as usize])
            .then(a.cmp(&b))
    });
    // the stop certificate's suffix max: the largest bound + slack
    // over every remaining pair of the plan (a per-pair slack does
    // not preserve the plan's bound order, so the whole remaining
    // suffix is folded once; NaN and +inf propagate upward through
    // total_cmp — the conservative never-stop cases).
    let mut dominated_suffix: Vec<f64> = vec![f64::NEG_INFINITY; order.len() + 1];
    for position in (0..order.len()).rev() {
        let flat = order[position] as usize;
        let dominated = pair_bounds[flat] + pair_slack[flat];
        dominated_suffix[position] = if dominated
            .total_cmp(&dominated_suffix[position + 1])
            == std::cmp::Ordering::Greater
        {
            dominated
        } else {
            dominated_suffix[position + 1]
        };
    }
    let mut computed: Vec<Option<f64>> = vec![None; pair_count];
    // THE BIT-EQUAL-ROW MEMO TABLE (the fourth reduction at the
    // bounded path): one slot per canonical flat, persistent across
    // the rounds — each distinct row-class pair sum is computed
    // exactly once (by the committed accumulation on the
    // representative rows, bit-identical to the direct per-pair
    // form by construction), and every later pair — in this round or
    // a later one — reads the shared value. With the memo off the
    // table stays empty and the rounds take the direct form.
    let mut memo_values: Vec<Option<f64>> =
        memo.map(|_| vec![None; pair_count]).unwrap_or_default();
    // THE PINNED-UNION OPTION. The ll matrix is mostly the elsewhere
    // constant (2.9% pinned at the tail locus), so in principle a pair
    // needs ~2,700 corrections, not 45,423 — 17x on that locus, 11.1x
    // territory-wide IN UNIT-PAIR COUNT.
    //
    // THE REALIZED SPEEDUP IS 2.1x, MEASURED, and the reason is a
    // tension worth knowing: the excess-mass bound keeps exactly the
    // pairs with HIGH excess mass (the pinned-heavy folds), so the pairs
    // it evaluates are the ones where sparsity buys least. The bound and
    // the union walk work against each other. The measured phase census
    // at T00048 (24 loci): classes_seconds 4.96 -> 2.33 (2.1x), while
    // qual_spectrum (16.0s) and qual_certificate (14.8s) are untouched —
    // and those, not the class scan, are the model's poles at an
    // ordinary territory.
    // These are THE exact-sum sites of the production (bounded) scan, so
    // this is where the work actually is. With IMPG_CLASS_SCAN_UNION the
    // union walk IS the value; with IMPG_CLASS_SCAN_COMPARE both forms are
    // computed, the exact one is kept, and the delta is receipted.
    let union_enabled = std::env::var("IMPG_CLASS_SCAN_UNION").is_ok();
    let union_compare = std::env::var("IMPG_CLASS_SCAN_COMPARE").is_ok();
    let union_fast = if union_enabled || union_compare {
        Some(ClassScanFastIndex::build(matrix, unit_counts))
    } else {
        None
    };
    let sum_of = |first: usize, second: usize| -> f64 {
        match union_fast.as_ref() {
            None => matrix.class_pair_sum(first, second, unit_counts),
            Some(fast) => {
                let fast_sum = fast.class_pair_sum(first, second, unit_counts);
                if union_enabled && !union_compare {
                    return fast_sum;
                }
                let exact = matrix.class_pair_sum(first, second, unit_counts);
                UNION_CENSUS_PAIRS.fetch_add(1, Ordering::Relaxed);
                UNION_CENSUS_SKIPPED
                    .fetch_add(fast.constant_units_skipped(first, second), Ordering::Relaxed);
                if exact.to_bits() != fast_sum.to_bits() {
                    UNION_CENSUS_DIFF.fetch_add(1, Ordering::Relaxed);
                    let ulp = ulp_distance(exact, fast_sum);
                    UNION_CENSUS_MAXULP.fetch_max(ulp, Ordering::Relaxed);
                }
                // COMPARE mode keeps the EXACT value (its whole point is a
                // no-answer-change census); union-only emits the fast one.
                if union_enabled && !union_compare {
                    fast_sum
                } else {
                    exact
                }
            }
        }
    };
    let mut running_best = f64::NEG_INFINITY;
    let mut next = 0usize;
    let mut chunk = 1usize;
    while next < order.len() {
        let end = (next + chunk).min(order.len());
        let sums: Vec<(u32, f64)> = match memo {
            None => order[next..end]
                .par_iter()
                .map(|&flat| {
                    let (first, second) = unrank_pair(flat as usize);
                    (flat, sum_of(first, second))
                })
                .collect(),
            Some(memo) => {
                // the round's pairs, keyed by their canonical flats
                let keyed: Vec<(u32, usize)> = order[next..end]
                    .iter()
                    .map(|&flat| {
                        let (first, second) = unrank_pair(flat as usize);
                        (flat, memo.canonical_flat(first, second))
                    })
                    .collect();
                // the distinct canonical sums this round still
                // misses — each computed exactly once, in parallel,
                // by the committed unit-order accumulation on the
                // representative rows
                let mut missing: Vec<usize> =
                    keyed.iter().map(|&(_, canon)| canon).collect();
                missing.sort_unstable();
                missing.dedup();
                missing.retain(|canon| memo_values[*canon].is_none());
                let fresh: Vec<(usize, f64)> = missing
                    .par_iter()
                    .map(|&canon| {
                        let (first, second) = unrank_pair(canon);
                        (canon, sum_of(first, second))
                    })
                    .collect();
                memo.note_computed(fresh.len() as u64);
                memo.note_hits(
                    keyed.len() as u64 - fresh.len() as u64,
                );
                for (canon, sum) in fresh {
                    memo_values[canon] = Some(sum);
                }
                // the round's own pairs read the shared values (the
                // flat slots keep the pair order — values shared,
                // never reordered)
                keyed.into_iter()
                    .map(|(flat, canon)| {
                        (flat, memo_values[canon].expect("the canonical sum is filled"))
                    })
                    .collect()
            }
        };
        for (flat, sum) in sums {
            computed[flat as usize] = Some(sum);
            if sum > running_best {
                running_best = sum;
            }
        }
        // strict stop, slack-corrected: every remaining pair's exact
        // sum is <= its bound + slack (the slack's own derivation),
        // so when the largest remaining bound + slack is strictly
        // below the running best, none can win or tie it.
        if end == order.len()
            || dominated_suffix[end].partial_cmp(&running_best)
                == Some(std::cmp::Ordering::Less)
        {
            break;
        }
        next = end;
        chunk *= 2;
    }
    (pair_bounds, pair_slack, computed, running_best, memo_values)
}

/// THE BOUNDED CLASS SCAN's STAGE-TWO KEEP THRESHOLD, derived from
/// the emitted consumers' own arithmetic (no tuned constants):
/// (a) the QUAL posterior terms need every pair whose exp(ll -
/// best) can be nonzero, i.e. every pair with ll - best above the
/// exp-underflow point — a pair with bound + slack < best + t has
/// exact < best + t, so its every QUAL contribution is exactly +0.0 (the
/// committed bit-identity argument); (b) the truth fields (the
/// receipt's truth_rank / truth_tied_classes / truth_log_likelihood)
/// need every pair that can be strictly above or tied with the
/// truth pair — a pair with bound + slack < truth_ll has exact < truth_ll,
/// so it is neither; where the truth pair sits INSIDE the QUAL
/// window (the receipted shape: the log gaps are tens-to-hundreds
/// of nats, under the 745-nat window) the truth threshold adds
/// nothing. The kept set is therefore {bound + the derived
/// rounding-error slack >= min(best + t, truth_ll)} — the minimum,
/// deterministic, sound for every emitted consumer (the slack is
/// SS10.6.22's fix: without it a pair with bound < best + t can
/// still carry exact >= best + t — the bound's f64 evaluation
/// rounds below the exact unit-order sum at some shapes, the
/// receipted 12-ulp shortfall at T03226 locus 16).
fn bounded_keep_threshold(best: f64, truth_ll: Option<f64>) -> f64 {
    let qual = best + exp_underflow_point();
    match truth_ll {
        Some(truth) => qual.min(truth),
        None => qual,
    }
}

/// THE FINAL-VALUE TIE-SET CONSTRUCTION — best and the called set
/// as pure functions of the class_lls vector, the exhaustive
/// form's own definition (every pair at best, in flat order). The
/// bounded scan emits through this AFTER stage two's kept window
/// (the T03226 locus-16 fix, SS10.6.22: the emitted set is a
/// function of the final exact values — which stage two has
/// already recovered bit-identical — never of the bound's
/// rounding); the exhaustive branch uses the same one spelling.
fn called_set_from_class_lls(class_lls: &[f64]) -> (f64, Vec<(usize, usize)>) {
    let best = class_lls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let called = (0..class_lls.len())
        .filter(|&flat| class_lls[flat] == best)
        .map(unrank_pair)
        .collect();
    (best, called)
}

/// The COLLINEARITY verdict of one placement's pins (the slice C
/// anatomy): the pinned anchors must sit at ONE offset for the
/// placement to be ONE CONTINUOUS path through the fold's spell (a
/// substitution-only generation spells the read along a single
/// offset; the piecewise serving of slice B can collect flanks at
/// mutually incompatible offsets — the flank-indel class). For
/// physically-forward placements the invariant is r_j - own_hull_j;
/// for physically-reverse placements the reflection r_j + own_hull_j
/// (the reverse geometry). A TRUE placement is always collinear: the
/// read is a contiguous segment of the row's source genome, the true
/// assignment is among the monotone assignments, so every pinned
/// anchor sits at its true position.
fn placement_offset(
    pinned: &[Option<i64>],
    own_hull: &[u64],
    forward: bool,
) -> (bool, Option<i64>) {
    let mut offsets: Vec<i64> = Vec::new();
    for (j, &r) in pinned.iter().enumerate() {
        if let Some(r) = r {
            let own = own_hull[j] as i64;
            offsets.push(if forward { r - own } else { r + own });
        }
    }
    match offsets.first() {
        None => (false, None),
        Some(&first) => (offsets.iter().all(|&o| o == first), Some(first)),
    }
}

/// The ONE-CONTINUOUS-PATH score: every base of the read projected
/// at the placement's single offset — the whole-read form of the
/// backbone + serving-anchor votes, INCLUDING the bases slice B
/// abstained inside unpinned anchor windows (at a collinear offset
/// the whole read projects, pinned or not). Requires collinear pins.
/// The third return is the out-of-extent base count (a measurement
/// for the anatomy; the placement-validity rule itself is
/// unchanged).
#[allow(clippy::too_many_arguments)]
fn single_offset_score(
    unit: &Unit,
    offset_const: i64,
    forward: bool,
    fold_seq: &[u8],
    fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
) -> (u32, u32, u32) {
    let read = &reads[unit.read];
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut out_of_extent = 0u32;
    for i in 0..read.len() as i64 {
        let d = i - unit.w_lo as i64;
        let coord = if forward {
            offset_const + d
        } else {
            offset_const + k as i64 - 1 - d
        };
        if coord < 0 || coord >= fold_len as i64 {
            out_of_extent += 1;
            continue;
        }
        let base = fold_seq[coord as usize];
        let read_base = if forward {
            read[i as usize]
        } else {
            complement(read[i as usize])
        };
        if read_base == base {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }
    (matches, mismatches, out_of_extent)
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

// (Phase 3 — the parallelism pilot: THE SEAM MAP. The two measured-
// dominant seams (the per-record binding loop and the per-locus loop)
// run on DEDICATED std::threads at a configurable width, NOT on the
// rayon pool. MEASURED REASON, both effects caught by the first 2-wide
// pilot run and its control: (1) the committed serial runs were never
// single-threaded in the loci loop — the per-locus scoring/classes/
// spectrum par_iters ride rayon's default pool (the whole box), and
// bounding RAYON_NUM_THREADS to the rung width silently serialized
// THOSE phases too (chrMT loci 1.8s serial-default -> 8.2s at a
// 2-wide pool); (2) nesting the seam par_iter in the SAME pool adds
// contention on top (18.5s: two concurrent loci each enqueueing
// thousands of inner jobs onto the same 2 workers). The seam map
// instead spawns `width` plain scoped threads, stride-assigned over
// the items; from a non-rayon thread every inner par_iter installs
// rayon's GLOBAL pool — the committed inner behavior, unchanged. The
// seam bodies are pure per-item functions over immutable shared state
// (compiler-enforced: the map takes `&(dyn Fn + Sync)`); each result
// lands in its own pre-sized slot and is read out in item order, so
// the outputs are interleaving-independent by construction — the
// byte-identity gate against the serial receipts is the race
// detector.)
/// The per-locus ingredients spill path (candidate C's streaming
/// emission): the receipts phase streams the locus's ingredients
/// line here, and the ordered emission below copies it into the
/// ingredients file and removes the spill — same bytes as the
/// committed buffered form, bounded RAM.
fn ingredients_spill_path(out: &Path, locus: u32) -> PathBuf {
    out.with_file_name(format!(
        "{}.ingredients.spill-{locus}.tmp",
        out.file_name().and_then(|name| name.to_str()).unwrap_or("")
    ))
}

fn seam_map<T: Sync, R: Send>(
    items: &[T],
    width: usize,
    f: &(dyn Fn(&T) -> io::Result<R> + Sync),
) -> io::Result<Vec<R>> {
    ensure(width >= 1, "the seam width must be at least 1")?;
    let n = items.len();
    if n == 0 || width == 1 {
        return items.iter().map(f).collect();
    }
    let out: std::sync::Mutex<Vec<Option<R>>> =
        std::sync::Mutex::new((0..n).map(|_| None).collect());
    let errors: std::sync::Mutex<Vec<(usize, io::Error)>> =
        std::sync::Mutex::new(Vec::new());
    let out_ref = &out;
    let errors_ref = &errors;
    std::thread::scope(|scope| {
        for w in 0..width {
            scope.spawn(move || {
                let mut i = w;
                while i < n {
                    match f(&items[i]) {
                        Ok(value) => out_ref.lock().unwrap()[i] = Some(value),
                        Err(error) => errors_ref.lock().unwrap().push((i, error)),
                    }
                    i += width;
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if let Some((_, error)) = errors.into_iter().min_by_key(|(i, _)| *i) {
        return Err(error);
    }
    let out = out.into_inner().unwrap();
    Ok(out
        .into_iter()
        .map(|slot| slot.expect("the seam map left an unset slot"))
        .collect())
}

// ---------------------------------------------------------------------------
// THE SPAN-CACHED LANE SEQUENCE SOURCE (phase 1, lever 2's call-bound fix,
// at module scope so the cache's own byte-identity contract is unit-tested
// directly). The load policy: whole-lane (the committed per-territory frame
// — the lanes are the locus.fa record slices) or span-scoped (the GLOBAL
// frame — the lanes are whole panel paths, so the cache holds the touched
// SPAN per lane, grown by union on demand; the derived bound is the crops
// actually requested, never a constant).
//
// THE DELTA-EXTEND GROWTH RULE (the union-regrowth re-decompression fix):
// on a miss the cached span is extended by ONLY the missing sub-ranges —
// fetch [new_lo, lo) and/or [hi, new_hi) — and stitched to the cached
// interior; the covered interior is never re-decompressed. The committed
// union-regrowth form re-fetched the WHOLE grown span on every miss (the
// measured disease: 277,796 loads fetching 56.7GB for L1's ~117.6MB unique
// territory envelope — ~482x re-decompression; SS10.5.1). BYTE-IDENTITY,
// the served-crop argument restated for the stitch: a crop of a
// Sources::fetch answer equals the direct fetch of the crop (the property
// the covered branch serves by), so fetch(new_lo, lo) is fetch(new_lo,
// new_hi)'s own [0, lo - new_lo) crop and fetch(hi, new_hi) is its
// [hi - new_lo, new_hi - new_lo) crop; the cached span is byte-equal to
// fetch(lo, hi) by the maintained invariant (base case: a direct fetch;
// inductive step: this argument), so the concatenation IS fetch(new_lo,
// new_hi) byte for byte — the buffer the union-regrowth form held. Every
// served crop is therefore a slice of the same Sources::fetch answer
// covering it, and the outputs are invariant; the identity gate re-proves
// it end-to-end.
// ---------------------------------------------------------------------------
struct LaneCache {
    /// The route sources, shared with the chunk's other territories when
    /// the caller derived them once per process (`ScoreInputsShared`):
    /// `fetch` takes `&self` (its memo is a Mutex), so a shared fetch
    /// returns the same bytes a private one would.
    sources: std::sync::Arc<routes::Sources>,
    lane_lengths: Vec<u64>,
    whole_lane: bool,
    seqs: Vec<
        std::sync::Mutex<
            Option<(u64, u64, std::sync::Arc<Vec<u8>>)>,
        >,
    >,
    loads: AtomicU64,
    bytes: AtomicU64,
    nanos: AtomicU64,
    /// THE PER-LOAD SEEK-PATTERN TRACE (SS10.6.13 duty 2's receipt
    /// instrument; default OFF — env IMPG_LANE_TRACE=<path> appends
    /// one JSONL line per REAL sources.fetch (the AGC seek): source,
    /// the fetched range, the delta kind, and the fetch's own wall.
    /// The bind seek-floor anatomy's measured half: the offline
    /// analysis classifies the loads (size histogram, per-lane run
    /// counts, adjacency) against the ~0.4ms per-fetch fixed floor
    /// the bind_seek_floor_anatomy_receipt measures at the AGC.
    lane_trace: std::sync::OnceLock<Option<std::sync::Mutex<std::io::BufWriter<std::fs::File>>>>,
}
impl LaneCache {
    fn trace_load(&self, source: usize, start: u64, end: u64, kind: &str, nanos: u128) {
        use std::io::Write as _;
        let trace = self
            .lane_trace
            .get_or_init(|| {
                std::env::var("IMPG_LANE_TRACE")
                    .ok()
                    .and_then(|path| {
                        std::fs::OpenOptions::new()
                            .create(true)
                            .append(true)
                            .open(path)
                            .ok()
                    })
                    .map(|file| std::sync::Mutex::new(std::io::BufWriter::new(file)))
            })
            .as_ref();
        if let Some(trace) = trace {
            let line = format!(
                "{{\"source\":{source},\"start\":{start},\"end\":{end},\
                 \"bytes\":{},\"kind\":\"{kind}\",\"nanos\":{nanos}}}\n",
                end.saturating_sub(start),
            );
            if let Ok(mut writer) = trace.lock() {
                let _ = writer.write_all(line.as_bytes());
            }
        }
    }
    fn lane(&self, source: usize) -> io::Result<std::sync::Arc<Vec<u8>>> {
        let started = Instant::now();
        let mut slot = self.seqs[source].lock().unwrap();
        if let Some((_, _, seq)) = slot.as_ref() {
            return Ok(std::sync::Arc::clone(seq));
        }
        let len = self.lane_lengths[source];
        let seq = std::sync::Arc::new(self.sources.fetch(source, 0, len)?);
        self.loads.fetch_add(1, Ordering::Relaxed);
        self.bytes.fetch_add(len, Ordering::Relaxed);
        self.nanos
            .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        *slot = Some((0, len, std::sync::Arc::clone(&seq)));
        Ok(seq)
    }
    fn span(
        &self,
        source: usize,
        start: u64,
        end: u64,
    ) -> io::Result<(u64, std::sync::Arc<Vec<u8>>)> {
        // The GLOBAL-frame policy: the cached span grows by union to
        // cover the requested crop; a covering span serves the crop
        // as a slice (the same bytes the direct fetch returns). The
        // span's OWN [lo, hi) bounds are returned WITH it: the span is
        // span-RELATIVE memory (a span grown from lo > 0 is not
        // addressable by absolute path coordinates — the fetch below
        // slices at (start - lo), and a covered crop equals the
        // direct fetch exactly).
        let started = Instant::now();
        let mut slot = self.seqs[source].lock().unwrap();
        let len = self.lane_lengths[source];
        if let Some((lo, hi, seq)) = slot.as_ref() {
            if *lo <= start && end <= *hi {
                return Ok((*lo, std::sync::Arc::clone(seq)));
            }
            let seq = std::sync::Arc::clone(seq);
            let (new_lo, new_hi) = ((*lo).min(start), (*hi).max(end));
            // SEGMENT-ALIGNED GROWTH (the seek-floor fix, SS10.6.13
            // duty 2): the AGC's per-fetch cost is a FIXED floor
            // (~0.4ms per touched segment at the panel AGC — the pack
            // re-read + whole-pack zstd decompress + LZ reference
            // prepare get_contig_range pays EVERY fetch, measured by
            // bind_seek_floor_anatomy_receipt: a 900B fetch costs the
            // same 0.4ms as a 6KB fetch), so growing the cached span
            // only to the requested edge re-pays the whole segment's
            // decompression on the next adjacent request. The fetch
            // below therefore grows to the touched SEGMENTS' own
            // boundaries — the archive's layout unit, from its public
            // descriptors, a derived bound (identity for
            // non-segment-addressed sources) — and every later
            // request inside the covered segments is a cache hit.
            // Byte-identity: the covered span is still a contiguous
            // [lo, hi) whose buffer is exactly fetch(lo, hi) by the
            // delta-extend induction (a crop of a fetch answer IS the
            // fetch of the crop — the Sources contract), and every
            // served crop is a slice of it; the alignment changes
            // only WHICH superset range is fetched, never a served
            // byte. The counters receipt the aligned volume honestly
            // (bytes rises toward the segment-hull of the touched
            // set, loads collapse — the fixed-cost floor is what the
            // load count multiplies).
            let (alo, ahi) = self.sources.aligned_crop(source, new_lo, new_hi);
            let mut grown: Vec<u8> = Vec::with_capacity((ahi - alo) as usize);
            if alo < *lo {
                let fetch_started = Instant::now();
                grown.extend_from_slice(&self.sources.fetch(source, alo, *lo)?);
                let fetch_nanos = fetch_started.elapsed().as_nanos();
                self.trace_load(source, alo, *lo, "left", fetch_nanos);
                self.loads.fetch_add(1, Ordering::Relaxed);
                self.bytes.fetch_add(*lo - alo, Ordering::Relaxed);
            }
            grown.extend_from_slice(&seq[..]);
            if ahi > *hi {
                let fetch_started = Instant::now();
                grown.extend_from_slice(&self.sources.fetch(source, *hi, ahi)?);
                let fetch_nanos = fetch_started.elapsed().as_nanos();
                self.trace_load(source, *hi, ahi, "right", fetch_nanos);
                self.loads.fetch_add(1, Ordering::Relaxed);
                self.bytes.fetch_add(ahi - *hi, Ordering::Relaxed);
            }
            self.nanos
                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            let grown = std::sync::Arc::new(grown);
            *slot = Some((alo, ahi, std::sync::Arc::clone(&grown)));
            return Ok((alo, grown));
        }
        // (The initial fetch takes the same alignment: the first
        // request in a segment pays the segment's decompression once
        // and every later request inside it is a slice serve.)
        let (alo, ahi) = self.sources.aligned_crop(source, start, end.min(len));
        let fetch_started = Instant::now();
        let seq = std::sync::Arc::new(self.sources.fetch(source, alo, ahi)?);
        let fetch_nanos = fetch_started.elapsed().as_nanos();
        self.trace_load(source, alo, ahi, "initial", fetch_nanos);
        self.loads.fetch_add(1, Ordering::Relaxed);
        self.bytes.fetch_add(ahi - alo, Ordering::Relaxed);
        self.nanos
            .fetch_add(fetch_nanos as u64, Ordering::Relaxed);
        *slot = Some((alo, ahi, std::sync::Arc::clone(&seq)));
        Ok((alo, seq))
    }
    fn fetch(&self, source: usize, start: u64, end: u64) -> io::Result<Vec<u8>> {
        // The Sources::fetch contract, mirrored exactly.
        ensure(
            source < self.lane_lengths.len()
                && start < end
                && end <= self.lane_lengths[source],
            "invalid source crop",
        )?;
        if self.whole_lane {
            let lane = self.lane(source)?;
            let lo = start as usize;
            let hi = end as usize;
            return Ok(lane[lo..hi].to_vec());
        }
        let (span_lo, span) = self.span(source, start, end)?;
        let lo = (start - span_lo) as usize;
        let hi = (end - span_lo) as usize;
        Ok(span[lo..hi].to_vec())
    }
}

/// THE INPUTS PHASE'S COMPONENT BREAKDOWN (IMPG_INPUTS_CENSUS=1,
/// receipt-side, DEFAULT OFF): sequential laps over the call's inputs
/// phase, emitted as one stderr line beside the phase's own total.
/// THE SHAREABILITY AUDIT'S INSTRUMENT: the parts are timed one after
/// another (never nested), so they sum to the phase total, and the
/// line states WHICH part is a pure function of (panel, routes) —
/// shareable by a whole chunk of territories — and which is
/// inherently per-territory. Default off: `lap` returns before any
/// Instant read, so no timing call executes and no receipt byte
/// changes.
pub struct InputsLaps {
    enabled: bool,
    started: Instant,
    laps: Vec<(&'static str, f64)>,
}
impl InputsLaps {
    pub fn start() -> Self {
        InputsLaps {
            enabled: std::env::var("IMPG_INPUTS_CENSUS").is_ok(),
            started: Instant::now(),
            laps: Vec::new(),
        }
    }
    pub fn lap(&mut self, component: &'static str) {
        if !self.enabled {
            return;
        }
        self.laps.push((component, self.started.elapsed().as_secs_f64()));
        self.started = Instant::now();
    }
    pub fn emit(&self, total_seconds: f64) {
        if !self.enabled {
            return;
        }
        let parts = self
            .laps
            .iter()
            .map(|(component, seconds)| format!("{component} {seconds:.2}s"))
            .collect::<Vec<_>>()
            .join(", ");
        eprintln!("[score] inputs breakdown: {parts} [total {total_seconds:.1}s]");
    }
}

/// THE CHUNK-INVARIANT SCORE INPUTS (the chunk's per-process derivation
/// share): the call's inputs phase splits along one line — the half that
/// is a pure function of (panel, routes) and the half that is the
/// territory's own. The first half is the panel's content identity (an
/// FNV over every sidecar: ~18.5GB of the HPRC panel), the routes graph,
/// the opened route sources, and the panel/route name maps; the SAME
/// values for every territory of a chunk, and a chunk runs its rows
/// (width 32 at the product form) inside ONE child process. The second
/// half — the territory's cut census, its locality partition maps, its
/// region resolution, its axis-node attribution — differs per territory
/// and stays where it is.
///
/// This struct is the first half, derived ONCE per process and captured
/// by every territory's call: the same shape as the fork-COW panel share
/// (`shared_panel` in the combined cut+call driver), where the chunk's
/// territories hold a reference instead of re-deriving. BYTE-IDENTITY:
/// every member is the value the per-territory derivation produced — the
/// same files hashed, the same graph parsed, the same name maps built —
/// and `sources` is shared behind an `Arc` because `Sources::fetch` takes
/// `&self` (its memo is a Mutex): a concurrent second territory's fetch
/// returns the same bytes as its own. No territory's product bit reads a
/// value this share did not compute.
pub struct ScoreInputsShared {
    /// The (panel, routes) pair this context was derived from: a call
    /// naming a different pair is refused loudly, never served another
    /// panel's identity.
    panel_prefix: String,
    routes_dir: PathBuf,
    identity: PanelIdentity,
    graph: routes::Graph,
    path_of_name: HashMap<String, usize>,
    source_of_name: HashMap<String, usize>,
    sources: std::sync::Arc<routes::Sources>,
}

impl ScoreInputsShared {
    /// Derive the chunk-invariant half over the panel the caller already
    /// holds loaded (the combined cut+call process's ONE panel load). The
    /// identity-versus-routes check the per-territory path performs lives
    /// here, once: it binds this panel to this routes inventory.
    pub fn load(
        panel: &SyngIndex,
        panel_prefix: &str,
        routes_dir: &Path,
        laps: &mut InputsLaps,
    ) -> io::Result<Self> {
        let identity = PanelIdentity::read(panel_prefix)?;
        laps.lap("panel-identity");
        let graph: routes::Graph = read_json(&routes_dir.join("graph.json"))?;
        if graph.panel != identity || !graph.generation_complete {
            return Err(invalid("incompatible route inventory"));
        }
        laps.lap("routes-graph");
        let lanes = graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect::<Vec<_>>();
        let sources = routes::Sources::open(&graph.source_paths, lanes)?;
        laps.lap("routes-sources");
        let path_of_name: HashMap<String, usize> = panel
            .name_map
            .name_to_path
            .iter()
            .map(|(name, &path)| (name.clone(), path as usize))
            .collect();
        let source_of_name: HashMap<String, usize> = graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.id))
            .collect();
        laps.lap("name-maps");
        Ok(ScoreInputsShared {
            panel_prefix: panel_prefix.to_string(),
            routes_dir: routes_dir.to_path_buf(),
            identity,
            graph,
            path_of_name,
            source_of_name,
            sources: std::sync::Arc::new(sources),
        })
    }

    /// The reuse gate: the shared inputs speak for exactly the (panel,
    /// routes) pair this call names — the values are held, never
    /// reinterpreted for another panel.
    fn ensure_names(&self, panel_prefix: &str, routes_dir: &Path) -> io::Result<()> {
        ensure(
            self.panel_prefix == panel_prefix && self.routes_dir == routes_dir,
            &format!(
                "the shared score inputs were derived from panel {:?} + routes {}; this \
                 call names panel {panel_prefix:?} + routes {}",
                self.panel_prefix,
                self.routes_dir.display(),
                routes_dir.display(),
            ),
        )
    }
}

/// THE INJECTED-INDEX SEAM (the shared-index cut+call lever): `run` is
/// the committed path — it loads the panel itself, exactly as before.
/// `run_injected` runs the SAME machinery over a panel the caller
/// already holds loaded (the combined cut+call process: the ONE
/// global-panel load serves the batch's cuts AND calls, deleting the
/// per-call inputs pole — the ~245s global load, 58-70% of the seam
/// rung's per-lane wall). The injected panel must be the panel named
/// by `options.panel` (the identity check below binds it to the
/// routes graph exactly as the loaded path does; a mismatched panel
/// fails loudly at that check).
pub fn run(options: Options) -> io::Result<()> {
    run_injected(options, None, None)
}
/// `shared_inputs` is the chunk-invariant half of the inputs phase, when
/// the caller (the combined cut+call driver) already derived it for this
/// process: `Some(context)` reuses it, `None` derives it here — the same
/// values either way (`ScoreInputsShared::ensure_names` binds the reuse
/// to the (panel, routes) pair this call names).
pub fn run_injected(
    options: Options,
    injected_panel: Option<&SyngIndex>,
    shared_inputs: Option<&ScoreInputsShared>,
) -> io::Result<()> {
    // (THE PER-CALL READ LENGTH — the concurrent-calls audit's
    // finding, the one mutable process-global the call path held:
    // every committed receipt's E is computed at the FRESH-PROCESS
    // value 150 BEFORE the derive cache stores the cache's own
    // declared length, and every later read (the skipped-positions
    // sites, the unit likelihood, the scoring receipt) reads the
    // DECLARED value. The value is therefore PER-CALL state — a
    // plain local, initialized fresh, stored once at the derive-cache
    // load — so the combined process's CONCURRENT calls cannot read
    // each other's length (the serial form's first-gate divergence,
    // the stale-148 class, is impossible by construction, not by
    // reset discipline).
    // (Deferred to the derive-cache load below, where the declared
    // length is stored — the only assignment; every read after it,
    // including the seam threads' closures, sees the declared
    // value.)
    let read_length;
    let started = Instant::now();
    let rss = RssGuard::new(options.rss_budget_gib);
    // (THE RECORD-ORIGIN GATE: the content-origin parameter belongs to
    // the GLOBAL frame alone — the committed local frame is a closed
    // frame that derives the substrate's content convention itself,
    // per record, against the AGC; a nonzero origin without the
    // global frame is a caller error, refused loudly.)
    ensure(
        options.record_origin == 0 || options.global_frame,
        "--record-origin requires --global-frame (the committed local frame \
         is closed: it never needed the absolute content origin)",
    )?;
    // (THE FAR-BIN PRUNE'S FRAME GATE: the prune is the global-frame
    // product default — the working form the §10.6.18 receipts gated
    // — and the committed LOCAL frame is byte-frozen; a default-on
    // prune reaching the local frame is refused LOUDLY, naming the
    // opt-out, never silently applied to the frozen path.)
    far_bin_prune_frame_gate(options.global_frame, options.far_bin_prune)?;
    // (Phase 3 — the parallelism pilot: THE SEAM SWITCH. The committed
    // default is the serial code path — byte-for-byte the phase-1
    // behavior. With IMPG_REALIGN_PARALLEL_SEAMS set (and the width in
    // IMPG_REALIGN_SEAM_WIDTH, default 2), the two measured-dominant
    // seams (the per-record binding loop and the per-locus loop) run
    // on dedicated seam threads — see `seam_map` for the measured
    // reasons the seams do NOT share the rayon pool with the per-locus
    // scoring/classes par_iters. The seam bodies are pure per-record /
    // per-locus computations whose results land in per-item slots read
    // out in item order, so the receipts are thread-interleaving-
    // independent BY CONSTRUCTION; the byte-identity gate against the
    // serial receipts is the race detector.)
    let seam_parallel = std::env::var("IMPG_REALIGN_PARALLEL_SEAMS").is_ok();
    let seam_width: usize = std::env::var("IMPG_REALIGN_SEAM_WIDTH")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(2);
    ensure(seam_width >= 1, "the seam width must be at least 1")?;
    // (Stage 2, the extent normalization and the alignment-induced
    // locality domain — rules 8 and 9 of the model doc above. The
    // conventions are pinned exactly once below: the selection flag
    // decides, with `Env` (the instrument default) deferring to the
    // committed environment switches bit-for-bit.)
    // (The alignment-induced locality domain, the owner's go: under
    // IMPG_REALIGN_LOCALITY_DOMAIN the maps' member lists ARE the
    // locality folds (the connected component of the axis partition's
    // pggb graph containing the axis row; the fold material = the
    // path-interval hull of the maximal monotone axis-correspondence,
    // built receipt-side by locality-domain-maps.py), and the locus's
    // records re-derive per locality: an occurrence votes iff its
    // interval overlaps a locality fold row on the occurrence's OWN
    // path — the census's own territory-touch convention over the
    // locality extents instead of the partition's extents. All the
    // standing machinery (the territory-image comparison rule, the
    // marginal realignment likelihood, the DP disciplines, the QUAL
    // bounds) is unchanged; the locality domain's folds make the
    // territory images exact (the fold IS the correspondence core),
    // so the locality mode implies the territory rule.)
    let (env_territory, env_locality) = (
        std::env::var("IMPG_REALIGN_TERRITORY_NORMALIZED").is_ok(),
        std::env::var("IMPG_REALIGN_LOCALITY_DOMAIN").is_ok(),
    );
    // (The alignment-connected fold relaxation — the segdup
    // follow-up's stage-1 rule, beside the exact fold rule: env
    // IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS. Rows fold by their
    // connection in the locality's alignment structure — the pggb
    // bubble containing the axis correspondence — subject to the
    // derived safety property (the members indistinguishable AT
    // READ LEVEL: bit-identical per-unit likelihood columns, the
    // variant pockets carrying zero observed read mass, or the fold
    // does NOT happen), and the truth-second correspondence is read
    // from the same alignment connection. The full derived rule is
    // documented at the merge pass. Env-unset: the exact rule's code
    // path verbatim — the committed receipts byte-identical.)
    let alignment_connected =
        std::env::var("IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS").is_ok();
    // (The call-core autopsy instrument, receipt-side, DEFAULT OFF:
    // IMPG_LOCUS_SCORING_CENSUS=1 wraps the per-locus scoring pair
    // loop's four CPU-work classes in per-call Instant timers,
    // atomically accumulated across the fold par_iter and emitted as
    // a "scoring_census" member of the per-locus receipt: the
    // per-orientation anchor/pin build (anchor_correspondence + the
    // merged backbone), placement_score (the pinned pair's own
    // alignment), the ALWAYS-ON direct_score exactness assertion
    // (candidate B's target), and the marginalization
    // (unit_log_likelihood + logsumexp2 + the row build). Default
    // off: no timing call executes, no receipt byte changes. The
    // honest rung is the SERIAL call (RAYON_NUM_THREADS=1, no seam),
    // where every measured second is a CPU second; at the timer's
    // own overhead (two Instant reads + one atomic add per class
    // call) the census run's walls are upper bounds, never used as
    // the committed receipt shape.)
    let locus_scoring_census = std::env::var("IMPG_LOCUS_SCORING_CENSUS").is_ok();
    // (THE LOCUS-PHASE CENSUS, receipt-side, DEFAULT OFF — the fold
    // anatomy's own instrument, the OPEN per-pair constant's receipt:
    // IMPG_LOCUS_PHASE_CENSUS=1 wraps the alignment-connected merge
    // pass, the class scan, and the QUAL cluster's own phases (the
    // winner's spectrum build, the tied-band certificate, the
    // degenerate cut-index lists, the cluster form, the best-outside
    // fallback) in per-phase Instant timers and emits them beside
    // the shape counters that price them (the certificate's band
    // entries, the per-fold multiset lengths, the spectrum entries)
    // as a "phase_census" member of the per-locus receipt — counters
    // computed at write time, never a raw dump. Default off: no
    // timing call executes, no receipt byte changes. The honest rung
    // is the SERIAL call (RAYON_NUM_THREADS=1), where every measured
    // second is a CPU second.)
    let locus_phase_census = std::env::var("IMPG_LOCUS_PHASE_CENSUS").is_ok();
    // (THE FOLD-SET DIGEST, gate-side, DEFAULT OFF: IMPG_FOLD_SET
    // _DIGEST=1 emits a per-locus "fold_set_digest" — an FNV-1a over
    // the FINAL fold space's ordered identities (each fold's ordered
    // member rows and length, in fold order — exactly the fields the
    // receipt's fold_identities carries, in that order, so the digest
    // is recomputable from any committed receipt) — so the fold-set
    // bit-identity gate compares DIGESTS computed in-process, never
    // two fold trees on disk. Default off: no receipt byte changes.)
    let fold_set_digest = std::env::var("IMPG_FOLD_SET_DIGEST").is_ok();
    // (THE OWNER'S RULING 1 — the exactness assertion OUT of the
    // production path: the ALWAYS-ON direct_score re-verification
    // (the §10.6.14 autopsy: 0.95 µs/pair, 32.9 core-s at L1 = 16%
    // of the ex-load call, 58% of the pinned alignment work) becomes
    // a receipt/gate/test instrument, env-gated:
    // IMPG_SCORE_EXACTNESS_ASSERT restores the committed always-on
    // form VERBATIM (every placement re-verified base-by-base, the
    // per-locus and whole-run FACTORIZATION gates enforced);
    // unset — the production default — the pair loop trusts the
    // factorized (m, c) exactly as the committed receipts'
    // assertion proved it, and the exactness SAMPLE receipt (the
    // independent checker's vote-level detail) plus the unit tests
    // keep proving the machinery. The assertion emits no product
    // bit: the byte-identity gate against the committed receipts is
    // the proof of the removal.)
    let exactness_assert = std::env::var("IMPG_SCORE_EXACTNESS_ASSERT").is_ok();
    // (THE OWNER'S RULING 3 — the accumulation-order TOGGLE at the
    // class scan — extended by the SS10.6.21 sanction: the BOUNDED
    // scan (the per-fold excess-mass branch-and-bound) is the
    // PRODUCTION DEFAULT (the owner: "take the better, faster
    // one"); IMPG_CLASS_SCAN_ORDER=exact restores the committed
    // exhaustive unit-order accumulation VERBATIM (the receipts'
    // form, the toggle doctrine); =sparse-skip is the legacy fast
    // accumulation (unchanged, superseded by the bounded default —
    // its receipted gain at the pruned form was -5.7%, under the
    // bounded scan's 98-99.8% at the receipted worst shapes);
    // =bounded is the explicit default form. Any other value is a
    // loud error, never a silent fallback. The --truth-qv
    // assessment mode keeps its EXHAUSTIVE form regardless (the
    // truth-rank threshold keeps nearly everything at the receipted
    // shapes, SS10.6.21's own caveat) — resolved at the scan site,
    // so the assessment receipts stay byte-identical.)
    let class_scan_order = match std::env::var("IMPG_CLASS_SCAN_ORDER") {
        Err(_) => ClassScanOrder::Bounded,
        Ok(value) if value == "exact" => ClassScanOrder::Exact,
        Ok(value) if value == "sparse-skip" => ClassScanOrder::SparseSkip,
        Ok(value) if value == "bounded" => ClassScanOrder::Bounded,
        Ok(other) => {
            return Err(invalid(&format!(
                "unknown IMPG_CLASS_SCAN_ORDER {other:?} \
                 (exact | sparse-skip | bounded)"
            )))
        }
    };
    // (The sparse-skip DELTA RECEIPT, receipt-side, DEFAULT OFF:
    // IMPG_CLASS_SCAN_COMPARE=1 computes BOTH the exact and the
    // sparse-skip sum for every class pair, receipts per locus how
    // many sums differ / the max ulp delta / the constant entries
    // skipped — and uses the EXACT value downstream, so the
    // compare form itself changes no output byte.)
    let class_scan_compare = std::env::var("IMPG_CLASS_SCAN_COMPARE").is_ok();
    // (THE BIT-EQUAL-ROW MEMOIZATION's toggle — the fourth
    // expressibility reduction, ON by default at every exact-sum
    // site of the class scan: IMPG_CLASS_SCAN_MEMO=off restores the
    // direct per-pair form verbatim (the bounded scan's own toggle
    // doctrine, the measurement restore path). The memo is
    // exact-preserving BY CONSTRUCTION — bit-equal rows share one
    // sum, the flat order and every emitted bit untouched — so no
    // output label is warranted; the per-shape hit rates are
    // receipted under IMPG_LOCUS_PHASE_CENSUS.)
    let scan_memo_enabled =
        std::env::var("IMPG_CLASS_SCAN_MEMO").map_or(true, |v| v != "off");
    // (THE UNIT-FOLD REDUCTION, model-side, DEFAULT OFF:
    // IMPG_UNIT_FOLD=1 folds the class scan's UNIT dimension — the
    // bit-equal-row memo's twin on the other axis. Two exactness-
    // preserving reductions, grouped by the EXACT f64 bit pattern of
    // each unit's whole ll column: a CONSTANT column (every fold
    // bit-equal) contributes count_total * c to every pair identically
    // (mix(c,c) = c exactly) so it leaves the per-pair loop as ONE
    // scalar, still counted in the reported class LL; units whose
    // columns are BIT-IDENTICAL share one term with their summed
    // count mass, in the surviving groups' original unit order.
    // Exact in REAL arithmetic (the same addends, nothing rounded), NOT
    // bit-identical (a different association of the fixed left-to-
    // right accumulation — measured at T00048 locus 8: 65,520 of
    // 65,535 kept pairs change bits, max 1,362 ulp, argmax unchanged by
    // 1.5e7x; docs/syng-gmem-bwt/unit-redundancy.md). Unset — the
    // committed unit-order walk, verbatim, and no receipt byte
    // changes; set — the labeled `unit_fold` receipt member and the
    // per-locus `[unit-fold]` stderr line state the mode and its
    // shape counters.)
    let unit_fold_enabled = std::env::var("IMPG_UNIT_FOLD").is_ok();
    // The selection flag pins the convention exactly: Env defers to
    // the committed switches verbatim (the instrument default);
    // Default pins the committed base convention (env ignored);
    // Territory/Locality pin the rules on (locality implies
    // territory — the locality fold IS the correspondence core).
    let locality_domain = match options.selection {
        Selection::Env => env_locality,
        Selection::Default => false,
        Selection::Territory | Selection::Locality => true,
    };
    let territory_normalized = match options.selection {
        Selection::Env => env_territory || env_locality,
        Selection::Default => false,
        Selection::Territory | Selection::Locality => true,
    };
    // THE SCORED LOCI - the request, resolved below once the
    // component's axis rows are loaded: the explicit --loci list,
    // or THE REGION-AS-INTERFACE FILTER (--region, the owner's
    // ruling: "detect any locus, any position, any allele"), a
    // chr:start-end range resolved against the component's own
    // frame to the covered windows. The two are mutually exclusive.
    enum LociRequest {
        Listed(Vec<u32>),
        Region(String),
    }
    let loci_request = if let Some(region) = &options.region {
        ensure(options.loci.is_empty(), "--region and --loci are mutually exclusive")?;
        LociRequest::Region(region.clone())
    } else {
        let loci: Vec<u32> = options
            .loci
            .split(',')
            .map(|s| s.trim().parse::<u32>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| invalid("invalid loci list"))?;
        ensure(!loci.is_empty(), "no loci requested")?;
        {
            let mut sorted = loci.clone();
            sorted.sort_unstable();
            sorted.dedup();
            ensure(sorted == loci, "loci must be distinct and sorted")?;
        }
        LociRequest::Listed(loci)
    };

    // ------------------------------------------------------------- inputs
    // (Phase instrumentation for the dominance measurement: timers and
    // counters only, all emitted to stderr; no receipt field changes,
    // no behavior change.)
    let inputs_started = Instant::now();
    let mut inputs_laps = InputsLaps::start();
    let fetch_calls = AtomicU64::new(0);
    let fetch_bytes = AtomicU64::new(0);
    let loaded_panel;
    let panel: &SyngIndex = match injected_panel {
        Some(panel) => panel,
        None => {
            loaded_panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
            &loaded_panel
        }
    };
    inputs_laps.lap("panel-load");
    // THE CHUNK-INVARIANT HALF (panel identity, routes graph, route
    // sources, name maps): reused when this process already derived it
    // for the chunk's other territories, derived here otherwise.
    let derived_inputs;
    let inputs: &ScoreInputsShared = match shared_inputs {
        Some(shared) => {
            shared.ensure_names(&options.panel, &options.routes)?;
            inputs_laps.lap("shared-inputs-reused");
            shared
        }
        None => {
            derived_inputs =
                ScoreInputsShared::load(panel, &options.panel, &options.routes, &mut inputs_laps)?;
            &derived_inputs
        }
    };
    let identity = &inputs.identity;
    let k = panel.syncmer_length_bp() as u64;
    let sources = &inputs.sources;
    let path_of_name = &inputs.path_of_name;
    let source_of_name = &inputs.source_of_name;
    // (Phase 1, lever 2 — the call-bound fix: a per-lane in-memory
    // sequence cache over the route sources. Phase 0 measured the
    // binding verify loop at ~90%/93% of the whole instrument wall,
    // CALL-bound at ~129-190us per ~72bp crop — the per-call
    // random-access seek/decompress of Sources::fetch into the AGC —
    // and measured the touched subset: one 93.3MB AGC file (3.34GB
    // uncompressed over 9,901 lanes), of which the exhaustive chrI
    // run touches 840 occurrence lanes (384.5MB) and the partition-
    // map member rows 1,286 lanes (514.7MB total). DERIVED CHOICE:
    // load each touched LANE once, whole, through the same validated
    // Sources::fetch path (one sequential decompression per lane —
    // the per-locality batch-prefetch option taken at its natural
    // locality, the contig; a window/block cache would still pay a
    // decompression call per 64KB window), then serve every crop as
    // a memory slice of the cached buffer. Byte-identity: the lane
    // buffer IS the 0..len crop Sources::fetch returns (length- and
    // alphabet-validated, uppercased); a crop of it is a crop of the
    // answer, so every fetch_seq call returns the same bytes as the
    // per-call AGC fetch it replaces — the identity gate proves it
    // end-to-end. The cache holds only touched lanes (<= ~515MB at
    // chrI), never the 3.34GB whole panel.)
    // (THE GLOBAL-PANEL CALL MODE's amendment: under --global-frame
    // the routes' lanes are WHOLE PANEL PATHS (a territory's member
    // rows span hundreds of chromosome-scale lanes), so the whole-
    // lane policy would hold the panel in RAM per territory. The
    // span-scoped policy serves the same bytes from the touched
    // SPAN per lane — the derived bound is the crops the call
    // itself requests; the committed local frame keeps the
    // whole-lane policy verbatim.)
    let lane_cache = LaneCache {
        lane_lengths: sources.lanes.iter().map(|&(_, len)| len).collect(),
        sources: std::sync::Arc::clone(sources),
        whole_lane: !options.global_frame,
        seqs: (0..inputs.graph.lanes.len())
            .map(|_| std::sync::Mutex::new(None))
            .collect(),
        loads: AtomicU64::new(0),
        bytes: AtomicU64::new(0),
        nanos: AtomicU64::new(0),
        lane_trace: Default::default(),
    };
    inputs_laps.lap("lane-cache-skeleton");
    let fetch_seq = |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
        let path_idx = *path_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the syng name map"))?;
        let path_len = panel.name_map.path_to_length[path_idx];
        let hi = hi.min(path_len);
        let source = *source_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the route sources"))?;
        if lo >= hi {
            return Ok(Vec::new());
        }
        fetch_calls.fetch_add(1, Ordering::Relaxed);
        let seq = lane_cache.fetch(source, lo, hi)?;
        fetch_bytes.fetch_add(seq.len() as u64, Ordering::Relaxed);
        Ok(seq)
    };

    // The census receipt (slice A's input, verbatim; not needed in the
    // frame-audit-only mode — the skeleton audit touches no reads).
    let mut census_records: Vec<CensusRecordLine> = Vec::new();
    if !options.frame_audit_only {
        {
            let (reader, _) = niffler::get_reader(Box::new(File::open(&options.census)?))
                .map_err(io::Error::other)?;
            for line in BufReader::new(reader).lines() {
                census_records.push(serde_json::from_str(&line?)?);
            }
        }
        ensure(
            census_records
                .windows(2)
                .all(|w| w[0].record + 1 == w[1].record),
            "census record ids are not dense",
        )?;
    }
    inputs_laps.lap("census");

    // The partition graph maps; the AXIS ROWS (the component path's
    // own member rows) ranked by start map to the component's window
    // ids in order (slice A's derivation, generalized: one window per
    // axis ROW, not per partition — a repeat-locality partition can
    // carry SEVERAL of the component's windows, first measured at
    // chrIV where 167 windows tile over 158 axis partitions; the
    // window->partition map may repeat a partition and every window
    // keeps its OWN axis row, which the truth-fold derivation below
    // matches exactly).
    let mut maps: BTreeMap<u32, PartitionMap> = BTreeMap::new();
    let mut axis_rows: Vec<(u32, u64, u64)> = Vec::new();
    for entry in std::fs::read_dir(&options.partition_graphs)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.ends_with(".gfa.map.json") {
            continue;
        }
        let text = std::fs::read_to_string(&path)?;
        let map: PartitionMap = serde_json::from_str(&text)?;
        for axis in map
            .members
            .iter()
            .filter(|m| m.path_name == options.component)
        {
            axis_rows.push((map.partition, axis.start, axis.end));
        }
        maps.insert(map.partition, map);
    }
    axis_rows.sort_by_key(|&(partition, start, end)| (start, end, partition));
    let n_windows = if options.frame_audit_only {
        axis_rows.len()
    } else {
        census_records
            .iter()
            .flat_map(|line| line.occurrences.iter())
            .flat_map(|occ| occ.partitions.iter().copied())
            .max()
            .map(|w| w as usize + 1)
            .unwrap_or(0)
    };
    ensure(
        axis_rows.len() == n_windows,
        "axis row count does not match the census window span",
    )?;
    let window_partition: Vec<u32> =
        axis_rows.iter().map(|&(partition, _, _)| partition).collect();
    // The window's OWN axis row (the truth-first fold is matched on
    // it, window-aware).
    let window_axis: Vec<(u64, u64)> =
        axis_rows.iter().map(|&(_, start, end)| (start, end)).collect();
    // THE REGION-AS-INTERFACE RESOLUTION (the owner's ruling: the
    // CLI's interface is COORDINATES - "detect any locus, any
    // position, any allele"): the region resolved against the
    // component's own axis rows to the covered windows - the same
    // --loci convention, computed instead of typed. (Resolved BEFORE
    // the global-frame translation: the region is typed in the
    // component's own LOCAL frame and the window ids are
    // frame-independent.)
    let loci: Vec<u32> = match loci_request {
        LociRequest::Listed(loci) => loci,
        LociRequest::Region(region) => {
            let loci = resolve_region_loci(&options.component, &window_axis, &region)?;
            eprintln!(
                "[score] region {region:?} resolved to loci {loci:?} (frame {})",
                options.component
            );
            loci
        }
    };
    // THE GLOBAL-PANEL FRAME TRANSLATION (the campaign's lever 1:
    // --panel names the ONE global panel; the locality maps are the
    // committed territory-local substrate). The record header IS the
    // frame statement: every member row "path:start-end" and the
    // component lane translate to the global frame by the interval
    // offset (conversion-not-materialization, the re-anchor
    // principle); the product emission translates back via the
    // inverse table. Default (unset): the maps verbatim, the
    // component verbatim — the committed per-territory path
    // byte-for-byte.
    let (territory_frame, panel_component): (
        Option<territory_frame::TerritoryFrame>,
        String,
    ) = if options.global_frame {
        let frame = territory_frame::TerritoryFrame::from_local_names(
            maps
                .values()
                .flat_map(|map| map.members.iter().map(|row| row.path_name.as_str())),
            options.record_origin,
        )?;
        let panel_component = frame.global_name(&options.component)?;
        for map in maps.values_mut() {
            for row in &mut map.members {
                let translated =
                    territory_frame::translate_row(&row.path_name, row.start, row.end, options.record_origin)?;
                row.path_name = translated.path_name;
                row.start = translated.start;
                row.end = translated.end;
            }
        }
        for row in axis_rows.iter_mut() {
            let translated =
                territory_frame::translate_row(&options.component, row.1, row.2, options.record_origin)?;
            row.1 = translated.start;
            row.2 = translated.end;
        }
        eprintln!(
            "[score] global-panel frame: component {:?} -> {panel_component:?} \
             over {} records (record origin +{})",
            options.component,
            frame.records_len(),
            options.record_origin,
        );
        (Some(frame), panel_component)
    } else {
        (None, options.component.clone())
    };
    // The window's OWN axis row, in the panel's frame (the machinery
    // consumes the translated extents under the global frame; the
    // region resolution above consumed the local ones).
    let window_axis: Vec<(u64, u64)> =
        axis_rows.iter().map(|&(_, start, end)| (start, end)).collect();
    // The partition membership per window: a partition hosting more
    // than one of the component's windows cannot attribute an SK1
    // member row to a window by placement alone (measured at chrIV's
    // partition 110, whose two SK1 rows are OTHER windows' orthologs)
    // — the truth-second fold is derived only at single-window
    // partitions, where the committed path-name rule is exact.
    let mut windows_of_partition: HashMap<u32, usize> = HashMap::new();
    for &partition in &window_partition {
        *windows_of_partition.entry(partition).or_insert(0) += 1;
    }
    for &locus in &loci {
        ensure((locus as usize) < n_windows, "locus outside the window span")?;
    }
    inputs_laps.lap("maps-frame-loci");
    // (Stage 2, rule 8's EDGE REFINEMENT: the per-window axis-node
    // attribution, built ONCE per component from the panel's own
    // stored walks over every window's axis row. A fold's image
    // extends to its own row edge over edge material that NO OTHER
    // window's axis row claims — the two measured boundary cases:
    // chrMT L2's winner-row head is 131/158 steps shared with the
    // NEIGHBOR window's walk (neighbor territory, excluded), chrV
    // L24's SK1 ortholog head is 12/12 steps shared with no window's
    // walk (a true edge variant, in-territory by elimination).
    // rep = the first window whose axis row carries the node;
    // count = how many windows' axis rows carry it. A node is
    // attributable to ANOTHER window iff rep != this locus or the
    // node is carried by more than one window.)
    let axis_node_attribution: Option<(HashMap<u32, u32>, HashMap<u32, u32>)> =
        if territory_normalized {
            let path_idx = *path_of_name
                .get(&panel_component)
                .ok_or_else(|| invalid("component path absent from the syng"))?;
            let path_len = panel.name_map.path_to_length[path_idx];
            let mut rep: HashMap<u32, u32> = HashMap::new();
            let mut count: HashMap<u32, u32> = HashMap::new();
            for (window, &(start, end)) in window_axis.iter().enumerate() {
                let end = end.min(path_len);
                ensure(start < end, "degenerate axis row")?;
                for (node, _bp) in panel.walk_path_range(path_idx, start, end)? {
                    let key = node.unsigned_abs();
                    rep.entry(key).or_insert(window as u32);
                    *count.entry(key).or_insert(0) += 1;
                }
            }
            Some((rep, count))
        } else {
            None
        };
    inputs_laps.lap("axis-attribution");
    let inputs_seconds = inputs_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] phase inputs: panel+routes+census+partition maps \
         ({} census records, {} windows) [{inputs_seconds:.1}s]",
        census_records.len(),
        n_windows,
    );
    inputs_laps.emit(inputs_seconds);
    rss.probe("inputs")?;

    // ------------------------------------- slice D gate 1: the frame audit
    // (the frame-blindness diagnosis, pure emission: per axis-partition
    // row of the given loci, the canonical-scheme pin skeleton beside
    // the stored path walk, with the per-position diff — the class
    // census of the frame-blinded anchors. Every canonical step is
    // sequence-verified in-process against the fetched panel sequence.)
    if options.frame_audit_only {
        let audit_started = Instant::now();
        let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq, k, false);
        for &locus in &loci {
            let partition = window_partition[locus as usize];
            let map = maps
                .get(&partition)
                .ok_or_else(|| invalid("axis partition map missing"))?;
            store.ensure_partition(&map.members, partition)?;
        }
        let mut audit = BufWriter::new(File::create(&options.out)?);
        let mut rows_total = 0u64;
        let mut rows_affected = 0u64;
        let mut steps_forward = 0u64;
        let mut steps_reverse = 0u64;
        let mut added_total = 0u64;
        let mut replaced_total = 0u64;
        let mut dropped_total = 0u64;
        for row in &store.folds {
            let diff = contained_skeleton_diff(row, k);
            rows_total += 1;
            steps_forward += diff.kept_forward;
            steps_reverse += diff.kept_reverse;
            added_total += diff.added;
            replaced_total += diff.replaced;
            dropped_total += diff.dropped;
            if diff.added > 0 || diff.replaced > 0 || diff.dropped > 0 {
                rows_affected += 1;
            }
            serde_json::to_writer(
                &mut audit,
                &json!({
                    "partitions": row.members.iter().map(|&(p, _)| p).collect::<Vec<_>>(),
                    "path_name": panel.name_map.path_to_name[row.path],
                    "start": row.start,
                    "end": row.end,
                    "length": row.end - row.start,
                    "steps": diff.canonical.iter().zip(diff.frames.iter())
                        .map(|(&(bp, node), &(_, frame))| json!([bp, node, frame]))
                        .collect::<Vec<_>>(),
                    "stored": diff.stored.iter()
                        .map(|&(bp, node)| json!([bp, node]))
                        .collect::<Vec<_>>(),
                    "added": diff.added,
                    "replaced": diff.replaced,
                    "dropped": diff.dropped,
                    "kept_forward": diff.kept_forward,
                    "kept_reverse": diff.kept_reverse,
                }),
            )?;
            writeln!(audit)?;
        }
        audit.flush()?;
        let partitions: BTreeSet<u32> = loci
            .iter()
            .map(|&l| window_partition[l as usize])
            .collect();
        eprintln!(
            "[score] frame audit: {} partitions, {} rows, {} affected \
             (added {added_total}, replaced {replaced_total}, dropped {dropped_total}; \
             kept forward {steps_forward}, reverse {steps_reverse}) [{:.1}s]",
            partitions.len(),
            rows_total,
            rows_affected,
            audit_started.elapsed().as_secs_f64(),
        );
        rss.probe("frame_audit")?;
        return Ok(());
    }
    eprintln!(
        "[score] inputs: {} census records, k {k}, {} windows, loci {loci:?}",
        census_records.len(),
        n_windows,
    );

    // ------------------------------------------- the scoring rate derivation
    let quality_started = Instant::now();
    let mut quality_bytes: BTreeSet<u8> = BTreeSet::new();
    let mut quality_counts: BTreeMap<u8, u64> = BTreeMap::new();
    let mut total_reads = 0u64;
    stream_fastq(&options.reads, |_seq, qual| {
        total_reads += 1;
        for &q in qual {
            quality_bytes.insert(q);
            *quality_counts.entry(q).or_insert(0) += 1;
        }
        Ok(())
    })?;
    // The sample's own quality statement: the mean per-base error
    // probability over every quality byte in the FASTQ (the yeast
    // validation's uniform-quality reads reduce EXACTLY to their single
    // Phred; real reads - the HG002 pilot - carry their measured
    // per-base qualities).
    ensure(
        !quality_bytes.is_empty(),
        "the sample carries no read qualities",
    )?;
    let mut base_count = 0u64;
    let mut error_sum = 0.0f64;
    for (&q, &count) in quality_counts.iter() {
        base_count += count;
        error_sum += count as f64 * 10f64.powf(-((q - b'!') as f64) / 10.0);
    }
    ensure(base_count > 0, "the sample carries no read qualities")?;
    let epsilon0 = error_sum / base_count as f64;
    ensure(epsilon0 > 0.0 && epsilon0 < 1.0, "degenerate read qualities")?;
    let phred = (-10.0 * epsilon0.log10()).round() as u8;
    let epsilon = 10f64.powf(-(phred as f64) / 10.0);
    let scoring = Scoring {
        a: (1.0 - epsilon).ln(),
        b: (epsilon / 3.0).ln(),
        epsilon,
        phred,
    };
    let quality_seconds = quality_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] scoring rate: uniform Phred {phred} over {total_reads} reads \
         (epsilon {epsilon}, A {}, B {}) [{quality_seconds:.1}s]",
        scoring.a,
        scoring.b,
    );
    rss.probe("quality")?;

    // ------------------------- the ELSEWHERE branch derivation (slice C)
    // The donor's genome size, measured from the panel (the
    // instrument's own genome universe — truth-free): the panel's
    // paths are grouped per (strain, haplotype) — "S288C#0#chrI",
    // "ATV#3#block43_contig1" — so the measured unit is the mean
    // total path length per HAPLOTYPE path set, doubled for the
    // diploid generation model (the balanced diploid sample's
    // documented construction: two haploid homologs).
    let mut strain_totals: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    let elsewhere_started = Instant::now();
    // THE GLOBAL-FRAME RECORD CONTEXT (the first L1 gate's named
    // divergence class, fixed by the supervisor ruling): the committed
    // local panel IS the territory's member record set — its name map
    // iterates exactly the ~251kb `path:start-end` slices, so the
    // committed E is measured at the TERRITORY'S RECORD context. The
    // global panel's name map iterates whole chromosomes (mean haploid
    // ~3.01 Gbp over the same strain families) — a DIFFERENT model
    // quantity: the first L1 gate receipted E -14.418156 (local, 251150.0
    // bp mean haploid) vs -23.810253 (global name map, 3011203926.4 bp)
    // — a uniform -9.392 per-component additive shift, QUAL-inflating
    // at 12+ loci. Under --global-frame the branch is therefore
    // computed at the SAME RECORD CONTEXT as the committed path: the
    // territory's own record statements (the frame table's records
    // enumeration — the frame-invariant extents `end - start`, grouped
    // by the same haplotype rule), reproducing the committed model
    // quantity exactly. The whole-chromosome context remains the
    // global-frame-NATIVE profile — arguably the more truthful
    // universe (the 251kb slices were the materialization-era
    // artifact) — named in the design doc as the OWNER's product
    // decision, never switched silently.
    if let Some(frame) = territory_frame.as_ref() {
        for (path, start, end) in frame.records() {
            let haplotype: String = path
                .splitn(3, '#')
                .take(2)
                .collect::<Vec<_>>()
                .join("#");
            let entry = strain_totals.entry(haplotype).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += end - start;
        }
    } else {
        for (name, &len) in panel
            .name_map
            .path_to_name
            .iter()
            .zip(panel.name_map.path_to_length.iter())
        {
            let haplotype: String = name
                .splitn(3, '#')
                .take(2)
                .collect::<Vec<_>>()
                .join("#");
            let entry = strain_totals.entry(haplotype).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += len;
        }
    }
    ensure(!strain_totals.is_empty(), "the panel has no paths")?;
    let panel_mean_haploid_bp: f64 = strain_totals
        .values()
        .map(|&(_, bp)| bp as f64)
        .sum::<f64>()
        / strain_totals.len() as f64;
    let genome_diploid_bp = 2.0 * panel_mean_haploid_bp;
    let elsewhere = elsewhere_log_prob(read_len(), scoring.a, genome_diploid_bp);
    let strain_list: Vec<serde_json::Value> = strain_totals
        .iter()
        .map(|(strain, &(paths, bp))| json!({"strain": strain, "paths": paths, "bp": bp}))
        .collect();
    eprintln!(
        "[score] elsewhere branch: E = {elsewhere:.6} \
         (panel strains {}, mean haploid {panel_mean_haploid_bp:.1} bp, \
         diploid {genome_diploid_bp:.1} bp{})",
        strain_totals.len(),
        if territory_frame.is_some() {
            "; the territory's RECORD context (the committed local panel's \
             own record set; the global panel's whole-chromosome context is \
             the named native-profile alternative)"
        } else {
            ""
        },
    );
    eprintln!(
        "[score] phase elsewhere: E derivation over {} haplotype sets \
         [{:.3}s]",
        strain_totals.len(),
        elsewhere_started.elapsed().as_secs_f64(),
    );
    let elsewhere_seconds = elsewhere_started.elapsed().as_secs_f64();

    // ------------------------------------- the read-record re-derivation
    let derive_started = Instant::now();
    ensure(options.derive_cache.exists(), "derive cache absent")?;
    let (key_tokens, key_reads, reads, read_multiplicity, raw_read_records): (
        Vec<Vec<u64>>,
        Vec<Vec<usize>>,
        Vec<Vec<u8>>,
        Vec<u32>,
        Vec<Vec<(u32, Vec<(i32, u32)>)>>,
    ) = {
        let (tokens, lists, seqs, mults, records) = read_derive_cache(&options.derive_cache)?;
        ensure(!seqs.is_empty(), "empty derive cache")?;
        let declared = seqs[0].len();
        ensure(
            declared > 0 && seqs.iter().all(|seq| seq.len() == declared),
            "cached reads do not share one length",
        )?;
        read_length = declared;
        (tokens, lists, seqs, mults, records)
    };
    let derive_seconds = derive_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] derive cache loaded: {} reads, {} keys ({derive_seconds:.1}s)",
        reads.len(),
        key_tokens.len(),
    );
    struct ReadRecord {
        key: u32,
        walk: Vec<(i32, u32)>,
    }
    let read_records: Vec<Vec<ReadRecord>> = raw_read_records
        .into_iter()
        .map(|records| {
            records
                .into_iter()
                .map(|(key, walk)| ReadRecord { key, walk })
                .collect()
        })
        .collect();
    rss.probe("derivation")?;

    // The identity gate (slice D): under IMPG_REALIGN_STORED_WALK_SKELETON
    // the pin skeleton stays the stored path walk (the slice-B/C
    // before-record; assessment-side diagnostic for the before/after
    // pairing). The repair — the canonical-scheme skeleton — is the
    // default.
    let stored_frame = std::env::var("IMPG_REALIGN_STORED_WALK_SKELETON").is_ok();
    let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq, k, stored_frame);
    let contig = panel_component
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| invalid("component lacks a contig suffix"))?
        .to_string();
    // (THE PANEL-SIDE PER-LOCUS PRECOMPUTE SERVE, SS10.6.23 slice 1 +
    // SS10.6.25 slice 2 — the labeled-toggle doctrine:
    // IMPG_PANEL_PRECOMPUTE unset serves when the asset beside the
    // partition maps is present and EVERY content key matches the
    // live inputs; =off restores the live path VERBATIM (the
    // IMPG_CLASS_SCAN_MEMO=off precedent); =on REQUIRES the serve (an
    // absent/mismatched asset is a loud error, the gate form —
    // including the slice-2 bind-scaffolding section: a format-/1
    // asset predates it and is refused under =on); any other value is
    // a loud error. A key mismatch under the default is a LOUD
    // LIVE-PATH FALLBACK (the routes-inventory precedent: never serve
    // stale structure, never fail a call because a panel-side asset
    // is stale). The serve changes NO emitted byte — the answer gates
    // are its proof — so no model label is warranted; it states
    // itself here and in the phase lines. The discovery sits ABOVE
    // the bind phase by slice 2's design: the bind scaffolding
    // serves from this same asset, before the first record binds.)
    let precompute_dir = options.partition_graphs.join("panel-precompute-loci");
    let precompute_mode = match std::env::var("IMPG_PANEL_PRECOMPUTE") {
        Err(_) => PanelPrecomputeMode::Auto,
        Ok(value) if value == "off" => PanelPrecomputeMode::Off,
        Ok(value) if value == "on" => PanelPrecomputeMode::On,
        Ok(other) => {
            return Err(invalid(&format!(
                "unknown IMPG_PANEL_PRECOMPUTE {other:?} (off | on)"
            )))
        }
    };
    let mut served: Option<PrecomputeServe> = None;
    if !matches!(precompute_mode, PanelPrecomputeMode::Off)
        && precompute_dir.join("manifest.json").is_file()
    {
        let maps_digest = partition_maps_digest(&options.partition_graphs)?;
        let routes_digest = routes_graph_digest(&options.routes.join("graph.json"))?;
        let attempt = (|| -> io::Result<()> {
            let keys = PanelPrecomputeKeys {
                panel: identity.to_owned(),
                k,
                component: panel_component.clone(),
                frame: if options.global_frame {
                    "global".to_string()
                } else {
                    "local".to_string()
                },
                record_origin: options.record_origin,
                stored_frame,
                maps_digest,
                routes_digest,
            };
            let asset = load_panel_precompute(&precompute_dir, &keys)?;
            ensure(
                asset.windows.len() == window_partition.len(),
                &format!(
                    "the panel-precompute asset carries {} windows against the call's {} \
                     (the window grid drifted)",
                    asset.windows.len(),
                    window_partition.len(),
                ),
            )?;
            // (SS10.6.25 slice 2 — THE BIND SCAFFOLDING SECTION: the
            // /2 format carries it; a /1 asset serves the row context
            // + fold layouts only and the bind derives live, loudly.
            // The =on gate form requires the section — the serve is
            // never partial without naming it.)
            let bind = match asset.bind {
                Some(paths) => {
                    Some(BindScaffoldServe::from_asset(paths, &path_of_name, panel)?)
                }
                None => {
                    ensure(
                        !matches!(precompute_mode, PanelPrecomputeMode::On),
                        "the panel-precompute asset predates the bind-scaffolding section \
                         (format /1) — rebuild it (IMPG_PANEL_PRECOMPUTE=on requires the \
                         serve)",
                    )?;
                    None
                }
            };
            // The served store, built from the asset rows in the
            // asset's row order (the builder's own creation order —
            // window order, member order), with the sequences left
            // in the sources (fetched live per representative row
            // at serve time, FNV-verified).
            for (row_index, row) in asset.rows.iter().enumerate() {
                let path_idx = *path_of_name
                    .get(&row.path_name)
                    .ok_or_else(|| invalid("the asset names a path absent from the syng"))?;
                ensure(
                    row.walk.len() == row.frames.len(),
                    "the asset's walk and frame tags disagree",
                )?;
                store.folds.push(RowFold {
                    path: path_idx,
                    start: row.start,
                    end: row.end,
                    seq: Vec::new(),
                    walk: row.walk.clone(),
                    frames: row.frames.clone(),
                    stored_walk: row.stored_walk.clone(),
                    members: row.members.clone(),
                });
                ensure(
                    store
                        .dedup
                        .insert((path_idx, row.start, row.end), row_index)
                        .is_none(),
                    "the asset carries a duplicate row identity",
                )?;
            }
            for (window, asset_window) in asset.windows.iter().enumerate() {
                ensure(
                    asset_window.partition == window_partition[window]
                        && asset_window.axis_start == window_axis[window].0
                        && asset_window.axis_end == window_axis[window].1,
                    &format!(
                        "the panel-precompute asset's window {window} disagrees with the \
                         live window grid (content keys drifted)",
                    ),
                )?;
                let map = maps
                    .get(&window_partition[window])
                    .ok_or_else(|| invalid("axis partition map missing"))?;
                ensure(
                    asset_window.member_fold.len() == map.members.len(),
                    &format!(
                        "the panel-precompute asset's member-fold map for window {window} \
                         does not match the partition maps (content keys drifted)",
                    ),
                )?;
            }
            served = Some(PrecomputeServe {
                windows: asset.windows,
                seq_fnv: asset.rows.iter().map(|row| row.seq_fnv).collect(),
                seqs: Mutex::new(HashMap::new()),
                bind,
            });
            Ok(())
        })();
        match attempt {
            Ok(()) => {
                let bind_state = match served.as_ref().and_then(|s| s.bind.as_ref()) {
                    Some(bind) => format!(
                        "; bind scaffolding: {} paths, {} ranges, {} steps",
                        bind.by_path.len(),
                        bind.by_path.values().map(|p| p.ranges.len()).sum::<usize>(),
                        bind.by_path.values().map(|p| p.steps.len()).sum::<usize>(),
                    ),
                    None => "; the bind-scaffolding section is absent — the bind derives \
                           live"
                        .to_string(),
                };
                eprintln!(
                    "[score] panel-precompute: SERVING the row context + fold layouts \
                     from {} ({} rows, {} windows, {} folds; panel identity, maps digest \
                     {:016x}, routes digest {:016x}){bind_state}",
                    precompute_dir.display(),
                    served.as_ref().map(|s| s.seq_fnv.len()).unwrap_or(0),
                    served.as_ref().map(|s| s.windows.len()).unwrap_or(0),
                    served
                        .as_ref()
                        .map(|s| s.windows.iter().map(|w| w.fold_rep.len()).sum::<usize>())
                        .unwrap_or(0),
                    maps_digest,
                    routes_digest,
                );
            }
            Err(e) if matches!(precompute_mode, PanelPrecomputeMode::On) => return Err(e),
            Err(e) => {
                eprintln!(
                    "[score] panel-precompute: NOT SERVED — {e}; the LIVE path runs \
                     (rebuild the asset or remove it; never a silent stale serve)",
                );
                store.folds.clear();
                store.dedup.clear();
                served = None;
            }
        }
    }

    // The routing's own canonical-scheme territory steps around one
    // panel position (ported from build_territory_index_rows, slice
    // A's binding lesson: the census occurrence starts live in THIS
    // convention, not the stored path walk's single frame) — extracted
    // VERBATIM into the one-spelling `canonical_steps_over_range`
    // (SS10.6.25 slice 2), so the live bind, the precompute builder,
    // and the serve's uncovered-position fallback all run THE SAME
    // code. THE BIND SCAFFOLDING SERVE: the scaffold (from the panel
    // asset's bind-scaffolding section, when it is present and keyed)
    // answers the canonical step at (path, position) from the
    // builder's recorded derivation when the position lies inside a
    // served band — the locality lemma's own exactness — and every
    // uncovered position (the census's own set is sample-side, A4's
    // own statement) derives LIVE, verbatim, counted loudly.
    let bind_scaffold: Option<&BindScaffoldServe> =
        served.as_ref().and_then(|serve| serve.bind.as_ref());

    // ------------------------- the pilot record set and the binding
    // The locus's evidence set: records with at least one occurrence
    // TOUCHING a scored window (the routing's own territory-touch rule).
    // (Phase instrumentation: the binding sub-timers split the phase
    // into its natural cost drivers — the per-record canonical-scheme
    // step extraction, the per-candidate key-shape re-derivation, and
    // the per-occurrence AGC fetch verification.)
    let binding_started = Instant::now();
    let binding_fetch_calls = fetch_calls.load(Ordering::Relaxed);
    let binding_fetch_bytes = fetch_bytes.load(Ordering::Relaxed);
    let binding_lane_loads = lane_cache.loads.load(Ordering::Relaxed);
    let binding_lane_bytes = lane_cache.bytes.load(Ordering::Relaxed);
    struct BindTimers {
        shape_seconds: f64,
        verify_seconds: f64,
        shape_calls: u64,
        verify_calls: u64,
        range_fetches: u64,
    }
    let mut bind_timers = BindTimers {
        shape_seconds: 0.0,
        verify_seconds: 0.0,
        shape_calls: 0,
        verify_calls: 0,
        range_fetches: 0,
    };
    let mut bind_steps_seconds = 0.0f64;
    let mut bind_steps_calls = 0u64;
    let mut bind_candidates = 0u64;
    let mut bind_occurrences = 0u64;
    let mut bind_shaped_keys: BTreeSet<u32> = BTreeSet::new();
    let mut pilot_records: BTreeSet<usize> = BTreeSet::new();
    // (THE FAR-BIN PRUNE, the owner's routing principle at the call's
    // census view: the non-adjacent multi-window records leave the
    // census view ENTIRELY — bind, variants, and every locus's
    // scoring — receipted here by count and mass. The prune is
    // per-record and the record ids stay the census's own dense ids,
    // so every locus the pruned records never touched is byte-
    // identical to the committed receipts by construction.)
    let mut far_bin_pruned_records = 0u64;
    let mut far_bin_pruned_mult = 0u64;
    let mut census_view_mult = 0u64;
    for (record, line) in census_records.iter().enumerate() {
        if options.far_bin_prune && census_record_far_bin_scattered(line) {
            far_bin_pruned_records += 1;
            far_bin_pruned_mult += line.multiplicity;
            continue;
        }
        census_view_mult += line.multiplicity;
        for occ in &line.occurrences {
            // (Locality mode: a record can vote at a pilot locus only
            // by overlapping a locality fold row, and every fold row
            // lies inside a window partition's own territory rows —
            // so the census's own window touch remains a valid
            // SUPERSET pre-filter over the full window set; the
            // exhaustive gate runs bind every touching record exactly
            // as the committed exhaustive runs do. The per-locus
            // selection below applies the exact locality overlap.)
            let touches = if locality_domain {
                !occ.partitions.is_empty()
            } else {
                occ.partitions.iter().any(|w| loci.contains(w))
            };
            if touches {
                pilot_records.insert(record);
                break;
            }
        }
    }
    eprintln!("[score] pilot records: {}", pilot_records.len());
    if options.far_bin_prune {
        let total_mult = census_view_mult + far_bin_pruned_mult;
        eprintln!(
            "[score] far-bin prune: {far_bin_pruned_records} records pruned, \
             census mult {far_bin_pruned_mult} of {total_mult} \
             ({:.1}%) — the non-adjacent-bin class",
            100.0 * far_bin_pruned_mult as f64 / total_mult.max(1) as f64,
        );
    }

    // Bind every pilot record to its derived key (slice A's binding,
    // verbatim, restricted to the pilot set).
    let decode_node = |token: u64| -> io::Result<i32> {
        let zigzag = token
            .checked_sub(2)
            .ok_or_else(|| invalid("record node token"))?
            / 2;
        Ok(((zigzag >> 1) as i64 ^ -(zigzag as i64 & 1)) as i32)
    };
    let mut key_by_first: HashMap<(i32, u64), Vec<u32>> = HashMap::new();
    let mut key_by_last: HashMap<(i32, u64), Vec<u32>> = HashMap::new();
    for (key, tokens) in key_tokens.iter().enumerate() {
        let first_node = decode_node(tokens[0])?;
        let last_node = decode_node(*tokens.last().unwrap())?;
        let anchor_count = (tokens.len() / 2 + 1) as u64;
        key_by_first
            .entry((first_node, anchor_count))
            .or_default()
            .push(key as u32);
        key_by_last
            .entry((last_node, anchor_count))
            .or_default()
            .push(key as u32);
    }
    fn key_shape(
        key: u32,
        key_tokens: &[Vec<u64>],
        key_reads: &[Vec<usize>],
        read_records: &[Vec<ReadRecord>],
        k: u64,
    ) -> io::Result<KeyShape> {
        let tokens = &key_tokens[key as usize];
        let canonical = decode_record_tokens(tokens)?;
        ensure(!canonical.is_empty(), "empty record walk")?;
        let span = walk_span(&canonical, k);
        let read = *key_reads[key as usize]
            .first()
            .ok_or_else(|| invalid("record key has no reads"))?;
        let own = read_records[read]
            .iter()
            .find(|record| record.key == key)
            .map(|record| record.walk.clone())
            .ok_or_else(|| invalid("read record mismatch"))?;
        let encoded = crate::sample_mem_bwt::encode_walk(
            &own
                .iter()
                .map(|&(node, pos)| (node, pos as u64))
                .collect::<Vec<_>>(),
        )?;
        let mirrored = encoded != tokens.as_slice();
        let w_lo = own.first().map(|&(_, p)| p as u64).unwrap_or(0);
        let own_positions = own_positions_of(&canonical, mirrored, w_lo);
        let m = canonical.len();
        ensure(own.len() == m, "own walk length differs from the canonical walk")?;
        for (t, &(own_node, own_pos)) in own.iter().enumerate() {
            let j = if mirrored { m - 1 - t } else { t };
            let expect_node = if mirrored { -canonical[j].0 } else { canonical[j].0 };
            ensure(
                own_node == expect_node && own_pos as u64 == own_positions[j],
                "the own walk disagrees with the mirror formula",
            )?;
        }
        Ok(KeyShape {
            canonical,
            span,
            read,
            mirrored,
            own_positions,
        })
    }
    // Sequence verification of one key at one occurrence: the hoisted
    // (phase 1, lever 1) union-range verify — see `verify_key_at` at
    // the top level. The fetch is performed once per (occurrence,
    // candidate) below, covering all three origin shifts.
    let mut census_key: Vec<Option<u32>> = vec![None; census_records.len()];
    let mut record_shifts: Vec<Vec<i64>> = vec![Vec::new(); census_records.len()];
    // (Phase 3 — the parallelism pilot, SEAM 1: the per-record binding
    // loop. The body below is a PURE function of the record: it reads
    // only immutable shared state (the census lines, the key indexes,
    // the derive cache, the panel, the reads) plus the thread-safe
    // fetch path (the atomic fetch counters and the mutex-per-lane
    // LaneCache, whose loads are serialized per lane by construction
    // and idempotent), and writes only its own result slot. Results
    // are collected in pilot_records order and merged serially, so
    // the bound keys and shifts are thread-interleaving-independent
    // by construction; the serial path is the committed phase-1 code
    // path with identical semantics, and the byte-identity gate
    // against the serial receipts is the race detector.)
    struct BindOne {
        key: Option<u32>,
        shifts: Vec<i64>,
        occurrences: u64,
        candidates: u64,
        shaped_keys: Vec<u32>,
        timers: BindTimers,
        steps_seconds: f64,
    }
    let bind_one = |record: usize| -> io::Result<BindOne> {
        let line = &census_records[record];
        let occ = line
            .occurrences
            .first()
            .ok_or_else(|| invalid("census record has no occurrences"))?;
        let occurrences = line.occurrences.len() as u64;
        let steps_started = Instant::now();
        let steps: Vec<i32> = match bind_scaffold.and_then(|scaffold| {
            scaffold.step_at(occ.path, occ.start, panel.name_map.path_to_length[occ.path])
        }) {
            Some(Some((_, node))) => vec![node],
            Some(None) => Vec::new(),
            None => canonical_steps_over_range(
                panel,
                &fetch_seq,
                k,
                occ.path,
                occ.start.saturating_sub(150),
                occ.start + 150,
            )?
            .into_iter()
            .filter(|&(bp, _)| bp == occ.start)
            .map(|(_, node)| node)
            .collect(),
        };
        let steps_seconds = steps_started.elapsed().as_secs_f64();
        ensure(!steps.is_empty(), "census occurrence start has no panel step")?;
        let mut candidates: BTreeSet<u32> = BTreeSet::new();
        for &step in &steps {
            let index: &HashMap<(i32, u64), Vec<u32>> = if occ.orientation == 0 {
                &key_by_first
            } else {
                &key_by_last
            };
            let node = if occ.orientation == 0 { step } else { -step };
            if let Some(keys) = index.get(&(node, line.anchors)) {
                for &key in keys {
                    candidates.insert(key);
                }
            }
        }
        let candidates_count = candidates.len() as u64;
        let shifts = [0i64, 1, -1];
        let mut timers = BindTimers {
            shape_seconds: 0.0,
            verify_seconds: 0.0,
            shape_calls: 0,
            verify_calls: 0,
            range_fetches: 0,
        };
        let verify_all = |key: u32, t: &mut BindTimers| -> io::Result<Option<Vec<i64>>> {
            let shape_started = Instant::now();
            let shape = key_shape(key, &key_tokens, &key_reads, &read_records, k)?;
            t.shape_seconds += shape_started.elapsed().as_secs_f64();
            t.shape_calls += 1;
            let mut per_occurrence = Vec::with_capacity(line.occurrences.len());
            for occ in &line.occurrences {
                let mut found: Vec<i64> = Vec::new();
                // (Phase 1, lever 1 — hoist the re-fetch: ONE range
                // fetch per (occurrence, candidate). Phase 0 measured
                // the three-shift loop as a 3.01x pure re-fetch
                // redundancy — the three per-shift ranges differ only
                // in the low bound, all inside
                // [max(0, start-2), start + span + 1). The union fetch
                // replaces the three per-shift fetches; each shift
                // then verifies against the exact slice its own
                // fetch returned before the rebuild.)
                let verify_started = Instant::now();
                let union_lo = (occ.start as i64 - 2).max(0) as u64;
                let union_hi = occ.start + shape.span + 1;
                let seq = fetch_seq(
                    &panel.name_map.path_to_name[occ.path],
                    union_lo,
                    union_hi,
                )?;
                for &shift in &shifts {
                    let verified =
                        verify_key_at(&shape, occ, &reads, &seq, union_lo, k, shift);
                    t.verify_calls += 1;
                    if verified {
                        found.push(shift);
                    }
                }
                t.verify_seconds += verify_started.elapsed().as_secs_f64();
                t.range_fetches += 1;
                match found.len() {
                    0 => return Ok(None),
                    1 => per_occurrence.push(found[0]),
                    _ => return Err(invalid(&format!(
                        "key verifies at two origin shifts (census record {record})"
                    ))),
                }
            }
            Ok(Some(per_occurrence))
        };
        let mut surviving: Vec<(u32, Vec<i64>)> = Vec::new();
        // Every CANDIDATE key is shaped before its verify (the serial
        // phase-1 code inserted each candidate into bind_shaped_keys
        // before the verify, surviving or not — the SET is the
        // semantic, collected here and unioned serially in the merge).
        let shaped_keys: Vec<u32> = candidates.iter().copied().collect();
        for key in candidates {
            if let Some(per_occurrence) = verify_all(key, &mut timers)? {
                surviving.push((key, per_occurrence));
            }
        }
        // Env-gated bind diagnostic (the held-out masking experiment's
        // walk-lock characterization): on a failed bind, dump the
        // record's first occurrence and every candidate key's shape
        // plus the panel slice the verify saw; for the first candidate
        // also scan EVERY occurrence and name the first one that fails
        // at all three origin shifts. Default (env unset): behavior
        // identical to the committed path.
        if surviving.len() != 1
            && std::env::var("IMPG_REALIGN_BIND_DEBUG")
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        {
            // The fail-closed receipt's own precision: name the
            // surviving keys (the ambiguity itself) before the
            // candidate-shape dump.
            eprintln!(
                "[bind-debug] record {record}: SURVIVING keys (the ambiguity): {:?}",
                surviving.iter().map(|(key, shifts)| (*key, shifts.clone())).collect::<Vec<_>>()
            );
            if let Some(&key) = shaped_keys.first() {
                if let Ok(shape) = key_shape(key, &key_tokens, &key_reads, &read_records, k) {
                    for (index, occ) in line.occurrences.iter().enumerate() {
                        let union_lo = (occ.start as i64 - 2).max(0) as u64;
                        let union_hi = occ.start + shape.span + 1;
                        let Ok(seq) = fetch_seq(
                            &panel.name_map.path_to_name[occ.path],
                            union_lo,
                            union_hi,
                        ) else {
                            eprintln!(
                                "[bind-debug] record {record}: occ {index}: fetch failed"
                            );
                            break;
                        };
                        let mut found = 0;
                        for &shift in &shifts {
                            if verify_key_at(&shape, occ, &reads, &seq, union_lo, k, shift) {
                                found += 1;
                            }
                        }
                        if found == 0 {
                            eprintln!(
                                "[bind-debug] record {record}: FIRST FAILING occ {index}: \
                                 path {} ({}) start {} orient {} span {} seq.len {}",
                                occ.path,
                                &panel.name_map.path_to_name[occ.path],
                                occ.start,
                                occ.orientation,
                                shape.span,
                                seq.len(),
                            );
                            break;
                        }
                    }
                }
            }
            for key in &shaped_keys {
                let Ok(shape) = key_shape(*key, &key_tokens, &key_reads, &read_records, k)
                else {
                    eprintln!("[bind-debug] record {record}: candidate {key}: shape error");
                    continue;
                };
                let seq = fetch_seq(
                    &panel.name_map.path_to_name[occ.path],
                    occ.start.saturating_sub(2),
                    occ.start + shape.span + 1,
                )
                .unwrap_or_default();
                let read = &reads[shape.read];
                eprintln!(
                    "[bind-debug] record {record}: occ path {} start {} orient {} \
                     intervals[0] {:?}; candidate {key}: anchors {} span {} mirrored \
                     {} read {} own_pos0 {:?} canonical0 {:?}; panel[0..16] {:?} \
                     read[own0..own0+16] {:?}",
                    occ.path,
                    occ.start,
                    occ.orientation,
                    occ.intervals.first().map(|x| &x[..]),
                    shape.canonical.len(),
                    shape.span,
                    shape.mirrored,
                    shape.read,
                    shape.own_positions.first(),
                    shape.canonical.first(),
                    String::from_utf8_lossy(&seq[..seq.len().min(16)]),
                    String::from_utf8_lossy(
                        &read[shape.own_positions[0] as usize
                            ..(shape.own_positions[0] as usize + 16).min(read.len())],
                    ),
                );
            }
        }
        ensure(
            surviving.len() == 1,
            &format!(
                "census record has no unique verifying key (census record {record}; \
                 {} surviving of {} candidates)",
                surviving.len(),
                candidates_count
            ),
        )?;
        Ok(BindOne {
            key: Some(surviving[0].0),
            shifts: surviving.swap_remove(0).1,
            occurrences,
            candidates: candidates_count,
            shaped_keys,
            timers,
            steps_seconds,
        })
    };
    let pilot_list: Vec<usize> = pilot_records.iter().copied().collect();
    let bound: Vec<BindOne> = if seam_parallel {
        seam_map(&pilot_list, seam_width, &|&record| bind_one(record))?
    } else {
        pilot_list
            .iter()
            .map(|&record| bind_one(record))
            .collect::<io::Result<Vec<_>>>()?
    };
    for (&record, one) in pilot_records.iter().zip(bound) {
        census_key[record] = one.key;
        record_shifts[record] = one.shifts;
        bind_occurrences += one.occurrences;
        bind_candidates += one.candidates;
        bind_steps_seconds += one.steps_seconds;
        bind_steps_calls += 1;
        bind_shaped_keys.extend(one.shaped_keys);
        bind_timers.shape_seconds += one.timers.shape_seconds;
        bind_timers.verify_seconds += one.timers.verify_seconds;
        bind_timers.shape_calls += one.timers.shape_calls;
        bind_timers.verify_calls += one.timers.verify_calls;
        bind_timers.range_fetches += one.timers.range_fetches;
    }
    eprintln!(
        "[score] pilot records bound: {} ({:.1}s)",
        pilot_records.len(),
        started.elapsed().as_secs_f64(),
    );
    eprintln!(
        "[score] phase binding: {} records / {} occurrences / {} candidate keys \
         ({} distinct shaped) / {} shape derivations / {} verify probes over \
         {} range fetches; steps {:.1}s ({} extractions), shapes {:.1}s, verify \
         {:.1}s, fetches {} ({:.1} MB), lane loads {} ({:.1} MB), phase {:.1}s; bind \
         scaffold: {} served / {} live-derived (uncovered)",
        pilot_records.len(),
        bind_occurrences,
        bind_candidates,
        bind_shaped_keys.len(),
        bind_timers.shape_calls,
        bind_timers.verify_calls,
        bind_timers.range_fetches,
        bind_steps_seconds,
        bind_steps_calls,
        bind_timers.shape_seconds,
        bind_timers.verify_seconds,
        fetch_calls.load(Ordering::Relaxed) - binding_fetch_calls,
        (fetch_bytes.load(Ordering::Relaxed) - binding_fetch_bytes) as f64 / (1024.0 * 1024.0),
        lane_cache.loads.load(Ordering::Relaxed) - binding_lane_loads,
        (lane_cache.bytes.load(Ordering::Relaxed) - binding_lane_bytes) as f64
            / (1024.0 * 1024.0),
        binding_started.elapsed().as_secs_f64(),
        bind_scaffold
            .map(|s| s.served.load(Ordering::Relaxed))
            .unwrap_or(0),
        bind_scaffold
            .map(|s| s.fallbacks.load(Ordering::Relaxed))
            .unwrap_or(pilot_records.len() as u64),
    );
    let binding_seconds = binding_started.elapsed().as_secs_f64();
    rss.probe("binding")?;

    // ----------------------- per-record read variants (slice A, verbatim)
    struct ReadVariant {
        mirror: bool,
        read: usize,
        w_lo: u64,
        count: u64,
    }
    let mut record_variants: HashMap<usize, Vec<ReadVariant>> = HashMap::new();
    let variants_started = Instant::now();
    for &record in &pilot_records {
        let line = &census_records[record];
        let key = census_key[record].unwrap();
        let mut variants: BTreeMap<(bool, u64, u64), ReadVariant> = BTreeMap::new();
        let mut instance_cursor: HashMap<usize, usize> = HashMap::new();
        for &read in &key_reads[key as usize] {
            let cursor = instance_cursor.entry(read).or_insert(0);
            let mut matching = read_records[read]
                .iter()
                .filter(|record| record.key == key)
                .skip(*cursor);
            let entry = matching
                .next()
                .ok_or_else(|| invalid("key_reads entry without a matching record"))?;
            *cursor += 1;
            let own = &entry.walk;
            let encoded = crate::sample_mem_bwt::encode_walk(
                &own
                    .iter()
                    .map(|&(node, pos)| (node, pos as u64))
                    .collect::<Vec<_>>(),
            )?;
            let mirrored = encoded != key_tokens[key as usize].as_slice();
            let w_lo = own.first().map(|&(_, p)| p as u64).unwrap_or(0);
            let vkey = (mirrored, w_lo, fnv1a64(&reads[read]));
            let weight = read_multiplicity[read] as u64;
            match variants.get_mut(&vkey) {
                Some(variant) => variant.count += weight,
                None => {
                    variants.insert(
                        vkey,
                        ReadVariant {
                            mirror: mirrored,
                            read,
                            w_lo,
                            count: weight,
                        },
                    );
                }
            }
        }
        let variants: Vec<ReadVariant> = variants.into_values().collect();
        let total: u64 = variants.iter().map(|v| v.count).sum();
        ensure(
            total == line.multiplicity,
            "variant instance counts do not sum to the census multiplicity",
        )?;
        record_variants.insert(record, variants);
    }
    eprintln!(
        "[score] phase variants: {} records, per-record read variants \
         [{:.1}s]",
        record_variants.len(),
        variants_started.elapsed().as_secs_f64(),
    );
    let variants_seconds = variants_started.elapsed().as_secs_f64();
    rss.probe("variants")?;

    // ------------------------------------------------------------------ receipts
    let mut report = BufWriter::new(File::create(&options.out)?);
    let mut exactness = BufWriter::new(File::create(&options.exactness_out)?);
    let mut anatomy = match &options.anatomy_out {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    let ingredients_path = format!("{}.ingredients.jsonl", options.out.display());
    let mut ingredients = BufWriter::new(File::create(&ingredients_path)?);
    let records_path = format!("{}.records.jsonl", options.out.display());
    let mut records_file = BufWriter::new(File::create(&records_path)?);
    // THE PRODUCT EMISSION (truth-free per-locus diplotype calls) and
    // the TEST-MODE artifact (the truth-referenced sequence QV,
    // derived from the product path: <product>.truth-qv.jsonl — a
    // SEPARATE stream, never mixed into the product output).
    let product_emission = options.product_calls.is_some();
    let mut product_calls = match &options.product_calls {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    let truth_qv_mode = !options.truth_qv.is_empty();
    ensure(
        !truth_qv_mode || product_emission,
        "--truth-qv requires --product-calls",
    )?;
    // THE DIPLOTYPE SEQUENCE EMISSION (the product output at every
    // called locus): <product>.sequences.jsonl beside the product
    // calls, truth-free. The product calls stream itself is
    // byte-identical whether or not this is set.
    ensure(
        !options.product_sequences || product_emission,
        "--product-sequences requires --product-calls",
    )?;
    let mut product_sequences_out = match (&options.product_calls, options.product_sequences) {
        (Some(path), true) => Some(BufWriter::new(File::create(
            format!("{}.sequences.jsonl", path.display()),
        )?)),
        _ => None,
    };
    ensure(
        !truth_qv_mode || options.truth_qv.len() == 2,
        "the diplotype truth-QV test mode requires exactly two truth haplotypes",
    )?;
    let mut truth_qv_out = match &options.product_calls {
        Some(path) if truth_qv_mode => Some(BufWriter::new(File::create(
            format!("{}.truth-qv.jsonl", path.display()),
        )?)),
        _ => None,
    };
    let truth_qv_aligner = match (truth_qv_mode, &truth_qv_out) {
        (true, Some(_)) => Some(std::sync::Mutex::new(QvAligner::spawn(&options.qv_biwfa)?)),
        _ => None,
    };
    // THE OFFICIAL SEQUENCE-QV CONVENTION (the owner's ruling, the
    // likegt max-qv form) — env-gated beside the committed biWFA
    // computation: IMPG_TRUTH_QV_ALLWAVE=1 adds the official
    // allwave-derived fields to the truth-QV test artifact (the
    // committed biWFA fields stay byte-identical and comparable);
    // IMPG_TRUTH_QV_ALLWAVE_BIN overrides the allwave binary
    // (default "allwave" from PATH). The default path (env unset)
    // is byte-identical to the committed receipts.
    let truth_qv_allwave = match (truth_qv_mode, std::env::var("IMPG_TRUTH_QV_ALLWAVE")) {
        (true, Ok(_)) => {
            let binary = std::env::var_os("IMPG_TRUTH_QV_ALLWAVE_BIN")
                .unwrap_or_else(|| "allwave".into());
            eprintln!(
                "[score] truth-QV official convention enabled (the likegt max-qv form): \
                 allwave binary {binary:?}, sparsification tree:5:0:0, QV cap 50"
            );
            Some(AllwaveQv {
                binary,
                // likegt's own max-qv thread default (-t 4, the
                // command's documented default).
                threads: 4,
            })
        }
        _ => None,
    };

    // The identity gate (slice D): under IMPG_REALIGN_STORED_WALK_SKELETON
    // the pin skeleton stays the stored path walk (the slice-B/C
    // before-record; assessment-side diagnostic for the before/after
    // pairing). The repair — the canonical-scheme skeleton — is the
    // default. (SS10.6.25 slice 2 moved the stored-frame convention,
    // the row store, and the serve discovery ABOVE the bind phase —
    // the bind scaffolding serves from the same discovered asset.)
    let mut skeleton = match &options.skeleton_out {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    // The axis partitions of every window the pilot records touch
    // (the locality classification needs their rows).
    // (Phase instrumentation: this is the CONTEXT ASSEMBLY phase —
    // per-row sequence fetch, stored-walk extraction, and the
    // canonical-scheme skeleton extraction with in-process sequence
    // verification. Under the precompute serve the phase is the
    // INDEX WALK ONLY: the partition maps' member rows resolve to the
    // asset's rows (the by-partition and by-partition-path indices,
    // the live insertion order verbatim), with the row content
    // already served and the representative rows' sequences fetched
    // per scored locus at serve time.)
    let rows_started = Instant::now();
    let rows_fetch_calls = fetch_calls.load(Ordering::Relaxed);
    let rows_fetch_bytes = fetch_bytes.load(Ordering::Relaxed);
    let rows_lane_loads = lane_cache.loads.load(Ordering::Relaxed);
    let rows_lane_bytes = lane_cache.bytes.load(Ordering::Relaxed);
    if served.is_some() {
        for &partition in &window_partition {
            if store.by_partition.contains_key(&partition) {
                continue;
            }
            let map = maps
                .get(&partition)
                .ok_or_else(|| invalid("axis partition map missing"))?;
            let mut indices = Vec::with_capacity(map.members.len());
            for member in &map.members {
                let path_idx = *path_of_name
                    .get(&member.path_name)
                    .ok_or_else(|| invalid("member path absent from the syng"))?;
                let row_index = *store
                    .dedup
                    .get(&(
                        path_idx,
                        member.start,
                        member.end.min(panel.name_map.path_to_length[path_idx]),
                    ))
                    .ok_or_else(|| {
                        invalid("the precompute asset lacks a member row (content keys drifted)")
                    })?;
                store
                    .by_partition_path
                    .entry((partition, path_idx))
                    .or_default()
                    .push(row_index);
                indices.push(row_index);
            }
            store.by_partition.insert(partition, indices);
        }
    } else {
        for &partition in &window_partition {
            let map = maps
                .get(&partition)
                .ok_or_else(|| invalid("axis partition map missing"))?;
            store.ensure_partition(&map.members, partition)?;
        }
    }
    eprintln!(
        "[score] phase rows: {} distinct rows over {} partitions ({} windows), context \
         assembly (canonical-scheme skeletons; fetches {} / {:.1} MB, \
         lane loads {} / {:.1} MB){} [{:.1}s]",
        store.folds.len(),
        windows_of_partition.len(),
        window_partition.len(),
        fetch_calls.load(Ordering::Relaxed) - rows_fetch_calls,
        (fetch_bytes.load(Ordering::Relaxed) - rows_fetch_bytes) as f64 / (1024.0 * 1024.0),
        lane_cache.loads.load(Ordering::Relaxed) - rows_lane_loads,
        (lane_cache.bytes.load(Ordering::Relaxed) - rows_lane_bytes) as f64
            / (1024.0 * 1024.0),
        if served.is_some() {
            " — SERVED from the panel-precompute asset (the index walk only; the \
             row content is the asset's own)"
        } else {
            ""
        },
        rows_started.elapsed().as_secs_f64(),
    );
    let rows_seconds = rows_started.elapsed().as_secs_f64();
    rss.probe("rows")?;

    let mut total_factorized_checked = 0u64;
    let mut total_factorized_equal = 0u64;
    let mut loci_seconds_total = 0.0f64;
    let skeleton_enabled = skeleton.is_some();
    let anatomy_enabled = anatomy.is_some();

    // (Phase 3 — the parallelism pilot, SEAM 2: the per-locus loop.
    // The body below is a PURE function of the locus: it reads only
    // immutable shared state (the partition maps, the assembled row
    // store, the census lines, the bound keys/shifts/variants, the
    // reads, the derive cache, the panel) plus the thread-safe fetch
    // path (atomic counters; the mutex-per-lane LaneCache) and the
    // read-only RSS guard — and it writes only its own per-locus
    // buffers. Every sidecar line is buffered per locus and written
    // in LOCUS ORDER after the loop, so parallel execution cannot
    // interleave lines and the receipts are thread-interleaving-
    // independent by construction. The serial path is the committed
    // phase-1 code path; the byte-identity gate against the serial
    // receipts is the race detector.)
    // The receipt's model string (hoisted: the product record carries
    // the same model identity; the receipt emission is unchanged).
    let model_name = if stored_frame {
        "anchor-realign-v2-marginal"
    } else if locality_domain {
        "anchor-realign-v2-marginal-frame-territory-locality"
    } else if territory_normalized {
        "anchor-realign-v2-marginal-frame-territory"
    } else {
        "anchor-realign-v2-marginal-frame"
    };
    // (The alignment-connected fold relaxation's marker — the
    // receipt states its own model honestly; env-unset keeps the
    // committed strings byte-identical.)
    let model_name: std::borrow::Cow<str> = if alignment_connected {
        std::borrow::Cow::Owned(format!("{model_name}-alignfold"))
    } else {
        std::borrow::Cow::Borrowed(model_name)
    };
    // (The class-scan order's own label — the fast mode states
    // itself in every emitted model string (the per-locus receipt,
    // the product calls, the sequences), so a sparse-skip output can
    // never be mistaken for an exact-mode one. Env-unset keeps the
    // committed strings byte-identical.)
    let model_name: std::borrow::Cow<str> =
        if matches!(class_scan_order, ClassScanOrder::SparseSkip) {
            std::borrow::Cow::Owned(format!("{model_name}-classscan-sparse"))
        } else {
            model_name
        };
    // (THE SANCTIONED BOUNDED CLASS SCAN's own label, the sparse-skip
    // precedent's shape (SS10.6.21): the bounded mode states itself
    // in every emitted model string — the per-locus receipt, the
    // product calls, the sequences — so a bounded output can never
    // be mistaken for an exhaustive one. THE PRODUCTION DEFAULT (the
    // owner: "take the better, faster one"); the --truth-qv
    // assessment mode keeps its EXHAUSTIVE form (the truth-rank
    // threshold keeps nearly everything at the receipted shapes), so
    // the assessment artifacts keep the committed model strings
    // byte-identical; IMPG_CLASS_SCAN_ORDER=exact restores the
    // exhaustive form with the committed unlabeled strings — the
    // toggle doctrine.)
    let class_scan_bounded =
        matches!(class_scan_order, ClassScanOrder::Bounded) && !truth_qv_mode;
    let model_name: std::borrow::Cow<str> = if class_scan_bounded {
        std::borrow::Cow::Owned(format!("{model_name}-classscan-bounded"))
    } else {
        model_name
    };
    // (The truth haplotypes as panel path names: under the global
    // frame a name carrying the territory-local interval form
    // translates with the maps — the committed CLI input shape is
    // preserved verbatim, and a name already in the panel's own
    // frame passes through.)
    let truth_qv_haplotypes: Vec<String> = options
        .truth_qv
        .iter()
        .map(|name| {
            if options.global_frame {
                territory_frame::translate_row(name, 0, 1, options.record_origin)
                    .map(|row| row.path_name)
                    .unwrap_or_else(|_| name.clone())
            } else {
                name.clone()
            }
        })
        .collect();
    let truth_qv_haplotypes = &truth_qv_haplotypes;
    struct LocusOutput {
        skeleton: Vec<u8>,
        exactness: Vec<u8>,
        records: Vec<u8>,
        // (Candidate C: the ingredients line is streamed to a
        // per-locus SPILL file at the receipts phase and copied into
        // the ordered ingredients file below — the ll_matrix JSON
        // never materializes in RAM. The final file's bytes are the
        // committed form's exactly.)
        ingredients_spill: PathBuf,
        report: Vec<u8>,
        anatomy: Vec<u8>,
        product: Vec<u8>,
        product_sequences: Vec<u8>,
        truth_qv: Vec<u8>,
        factorized_checked: u64,
        factorized_equal: u64,
        locus_seconds: f64,
    }
    // (The serve reference handed to the per-locus seam: the serve is
    // read-only shared state (the window groupings, the rows' content
    // keys, the Mutex-guarded sequence cache), Sync by construction.)
    let served_precompute = served.as_ref();
    let compute_locus = |locus: u32| -> io::Result<LocusOutput> {
        let locus_started = Instant::now();
        let mut skeleton_buf: Vec<u8> = Vec::new();
        let mut exactness_buf: Vec<u8> = Vec::new();
        let mut records_buf: Vec<u8> = Vec::new();
        let mut report_buf: Vec<u8> = Vec::new();
        let mut anatomy_buf: Vec<u8> = Vec::new();
        let mut product_buf: Vec<u8> = Vec::new();
        let mut product_sequences_buf: Vec<u8> = Vec::new();
        let mut truth_qv_buf: Vec<u8> = Vec::new();
        let window = locus as usize;
        let partition = window_partition[window];
        let (axis_start, axis_end) = window_axis[window];
        let map = &maps[&partition];

        // ------------------------------------------------------- the folds
        // Fold the partition's member rows by identical sequence AND
        // walk (coalesced copies are ONE hypothesis). Under the
        // panel-precompute serve (SS10.6.23 slice 1) the grouping is
        // SERVED from the asset — the masking commute lemma's regroup
        // over the live member list, the fold content derived through
        // the same one-spelling helpers, the representative rows'
        // sequences fetched live and FNV-verified; the live path
        // (and IMPG_PANEL_PRECOMPUTE=off) runs the extracted
        // build_window_folds_live VERBATIM — one spelling of the
        // derivation, two sources of the grouping.
        let fold_started = Instant::now();
        let (folds, excluded_rows, flank_excluded_rows): (Vec<Fold>, u64, u64) =
            if let Some(serve) = served_precompute {
                serve.serve_window_folds(
                    &store.folds,
                    window,
                    &map.members,
                    &panel_component,
                    axis_start,
                    axis_end,
                    options.flank_overlap,
                    &options.exclude_paths,
                    k,
                    &fetch_seq,
                    &panel.name_map.path_to_name,
                )?
            } else {
                let (folds, _rep_rows, _member_fold, excluded_rows, flank_excluded_rows) =
                    build_window_folds_live(
                        &store.folds,
                        &store.dedup,
                        &map.members,
                        &panel_component,
                        axis_start,
                        axis_end,
                        options.flank_overlap,
                        &options.exclude_paths,
                        &path_of_name,
                        &|path_idx| panel.name_map.path_to_length[path_idx],
                        k,
                    )?;
                (folds, excluded_rows, flank_excluded_rows)
            };
        let n_folds = folds.len();
        // ------------------------------------------- the territory images
        // (Stage 2, rule 8: the locus's normalized comparison domain.
        // The territory is the LOCUS'S OWN COORDINATE — the window's
        // axis row (the component's own member row at [axis_start,
        // axis_end)); per fold the image is the interval spanned by
        // the fold's walk steps shared with the axis row's walk —
        // the alignment-induced correspondence. Under the default
        // convention every fold's image is its FULL EXTENT and the
        // rule degenerates exactly to the committed comparison.)
        let territory_axis_fold: Option<usize>;
        let territory: Vec<(u64, u64)> = if territory_normalized {
            let hits: Vec<usize> = folds
                .iter()
                .enumerate()
                .filter(|(_, fold)| {
                    fold.members.iter().any(|m| {
                        m.path_name == panel_component
                            && m.start == axis_start
                            && m.end == axis_end
                    })
                })
                .map(|(index, _)| index)
                .collect();
            ensure(
                hits.len() == 1,
                "the window's axis row is not exactly one fold at this locus",
            )?;
            let axis_nodes: BTreeSet<u32> = folds[hits[0]]
                .walk
                .iter()
                .map(|&(_, node)| node.unsigned_abs())
                .collect();
            territory_axis_fold = Some(hits[0]);
            // (The EDGE REFINEMENT's attribution witness: node ->
            // (first window's axis row carrying it, how many windows'
            // axis rows carry it) — a fold's edge pocket extends to
            // its own row edge iff NO step in the pocket is claimed
            // by ANOTHER window's axis row.)
            let (rep, count) = axis_node_attribution
                .as_ref()
                .expect("attribution built under the normalized convention");
            let attributable_elsewhere = |node: u32| -> bool {
                match (rep.get(&node), count.get(&node)) {
                    (Some(&first), Some(&windows)) => first != locus || windows > 1,
                    _ => false,
                }
            };
            folds
                .iter()
                .enumerate()
                .map(|(fold_index, fold)| {
                    // The window's own axis fold IS the territory by
                    // definition: its image is its full extent.
                    if fold_index == hits[0] {
                        return (0u64, fold.len);
                    }
                    let shared: Vec<u64> = fold
                        .walk
                        .iter()
                        .filter(|&&(_, node)| axis_nodes.contains(&node.unsigned_abs()))
                        .map(|&(bp, _)| bp)
                        .collect();
                    let (lo, hi) = match (shared.first(), shared.last()) {
                        (Some(&lo), Some(&hi)) => (lo, hi + k),
                        // No shared anchor: the whole row is non-own
                        // material — own BY ELIMINATION iff no other
                        // window's axis row claims any of its steps
                        // (the empty image otherwise).
                        _ => {
                            if fold
                                .walk
                                .iter()
                                .all(|&(_, node)| !attributable_elsewhere(node.unsigned_abs()))
                            {
                                return (0u64, fold.len);
                            }
                            return (0u64, 0u64);
                        }
                    };
                    // The head pocket (steps before the first shared
                    // anchor): extend to the row edge iff no step is
                    // another window's territory material.
                    let head_claimed = fold.walk.iter().any(|&(bp, node)| {
                        bp < lo && attributable_elsewhere(node.unsigned_abs())
                    });
                    let tail_claimed = fold.walk.iter().any(|&(bp, node)| {
                        bp >= hi && attributable_elsewhere(node.unsigned_abs())
                    });
                    (
                        if head_claimed { lo } else { 0 },
                        if tail_claimed { hi } else { fold.len },
                    )
                })
                .collect()
        } else {
            territory_axis_fold = None;
            folds.iter().map(|fold| (0u64, fold.len)).collect()
        };
        // The pin-skeleton sidecar (slice D's frame-repair audit): per
        // fold, the canonical-scheme contained steps with frame tags
        // beside the stored walk's, with the position diff — and the
        // per-locus totals for the receipt. Under the identity gate
        // (stored_frame) the canonical rows were not extracted for
        // scoring, so the sidecar is not emitted there.
        let mut skeleton_rows = 0u64;
        let mut skeleton_affected = 0u64;
        let mut skeleton_added = 0u64;
        let mut skeleton_replaced = 0u64;
        let mut skeleton_dropped = 0u64;
        let mut skeleton_forward = 0u64;
        let mut skeleton_reverse = 0u64;
        if skeleton_enabled {
            for (fold_index, fold) in folds.iter().enumerate() {
                let first = &fold.members[0];
                let path_idx = path_of_name[&first.path_name];
                let clipped_end =
                    first.end.min(panel.name_map.path_to_length[path_idx]);
                let row_index = *store
                    .dedup
                    .get(&(path_idx, first.start, clipped_end))
                    .ok_or_else(|| invalid("fold row missing from the store"))?;
                let row = &store.folds[row_index];
                let diff = contained_skeleton_diff(row, k);
                skeleton_rows += 1;
                skeleton_added += diff.added;
                skeleton_replaced += diff.replaced;
                skeleton_dropped += diff.dropped;
                skeleton_forward += diff.kept_forward;
                skeleton_reverse += diff.kept_reverse;
                if diff.added > 0 || diff.replaced > 0 || diff.dropped > 0 {
                    skeleton_affected += 1;
                }
                serde_json::to_writer(
                    &mut skeleton_buf,
                    &json!({
                        "locus": locus,
                        "partition": partition,
                        "fold": fold_index,
                        "members": fold.members.iter().map(|m| json!({
                            "path_name": m.path_name, "start": m.start, "end": m.end,
                        })).collect::<Vec<_>>(),
                        "length": fold.len,
                        "steps": diff.canonical.iter().zip(diff.frames.iter())
                            .map(|(&(bp, node), &(_, frame))| json!([bp, node, frame]))
                            .collect::<Vec<_>>(),
                        "stored": diff.stored.iter()
                            .map(|&(bp, node)| json!([bp, node]))
                            .collect::<Vec<_>>(),
                        "added": diff.added,
                        "replaced": diff.replaced,
                        "dropped": diff.dropped,
                        "kept_forward": diff.kept_forward,
                        "kept_reverse": diff.kept_reverse,
                    }),
                )?;
                skeleton_buf.push(b'\n');
            }
        }
        eprintln!(
            "[score] locus {locus}: partition {partition}, {} members -> {n_folds} folds{} \
             [{:.1}s]",
            map.members.len(),
            if served_precompute.is_some() {
                " (SERVED from the panel-precompute asset)"
            } else {
                ""
            },
            fold_started.elapsed().as_secs_f64(),
        );

        // ---------------------------------- the locus's records and units
        let units_started = Instant::now();
        // (Locality mode: the territory-touch convention on the NEW
        // locality extents — an occurrence votes at the locus iff its
        // interval overlaps a locality fold row on the occurrence's
        // OWN path, the census's own coordinate-overlap rule over the
        // locality's rows. The committed partition mode keeps the
        // census's window ids verbatim.)
        let locality_touch = |position: usize, occ: &CensusOccurrence,
                               span: u64,
                               shifts: &Vec<i64>| -> bool {
            if !locality_domain {
                return occ.partitions.contains(&locus);
            }
            let origin =
                ((occ.start as i64 + shifts[position]).max(0)) as u64;
            store
                .rows_of_partition_path(locus, occ.path)
                .iter()
                .any(|&row_index| {
                    let row = &store.folds[row_index];
                    row.start < origin + span && origin < row.end
                })
        };
        let mut locus_records: Vec<usize> = Vec::new();
        for &record in &pilot_records {
            let line = &census_records[record];
            let span = walk_span(
                &key_shape(
                    census_key[record].unwrap(),
                    &key_tokens,
                    &key_reads,
                    &read_records,
                    k,
                )?
                .canonical,
                k,
            );
            let shifts = &record_shifts[record];
            if line
                .occurrences
                .iter()
                .enumerate()
                .any(|(position, occ)| locality_touch(position, occ, span, shifts))
            {
                locus_records.push(record);
            }
        }
        // Per record: the touching occurrences' orientations, covered
        // keys, and the locality classification.
        struct TouchingData {
            orientations: BTreeSet<u8>,
            covered_nodes: BTreeSet<u64>,
            covered_edges: BTreeSet<u64>,
            in_axis: u64,
            overhang: u64,
            extension: u64,
        }
        let universe_nodes: BTreeSet<u64> =
            folds.iter().flat_map(|f| f.nodes.iter().map(|n| n.0)).collect();
        let universe_edges: BTreeSet<u64> =
            folds.iter().flat_map(|f| f.edges.iter().map(|e| e.0)).collect();
        let mut touching: HashMap<usize, TouchingData> = HashMap::new();
        let mut locality_examples: Vec<serde_json::Value> = Vec::new();
        for &record in &locus_records {
            let line = &census_records[record];
            let span = walk_span(
                &key_shape(
                    census_key[record].unwrap(),
                    &key_tokens,
                    &key_reads,
                    &read_records,
                    k,
                )?
                .canonical,
                k,
            );
            let shifts = &record_shifts[record];
            let mut data = TouchingData {
                orientations: BTreeSet::new(),
                covered_nodes: BTreeSet::new(),
                covered_edges: BTreeSet::new(),
                in_axis: 0,
                overhang: 0,
                extension: 0,
            };
            for (position, occ) in line.occurrences.iter().enumerate() {
                if !locality_touch(position, occ, span, shifts) {
                    continue;
                }
                data.orientations.insert(occ.orientation);
                for interval in &occ.intervals {
                    for &node in interval {
                        if universe_nodes.contains(&(node as u64)) {
                            data.covered_nodes.insert(node as u64);
                        }
                    }
                    for pair in interval.windows(2) {
                        let edge = pack_edge(pair[0], pair[1]);
                        if universe_edges.contains(&edge) {
                            data.covered_edges.insert(edge);
                        }
                    }
                }
                // The locality classification (slice A's own-row rule):
                // containment in a member row of the locus's domain
                // rows, on the occurrence's own path. (Locality mode:
                // the domain rows are the locality folds; the committed
                // partition mode keeps every touched window's axis
                // partition rows.)
                let origin = (occ.start as i64 + shifts[position]).max(0) as u64;
                let mut contained = false;
                let mut overlapped = false;
                if locality_domain {
                    for &row_index in store.rows_of_partition_path(locus, occ.path) {
                        let row = &store.folds[row_index];
                        if row.start <= origin && origin + span <= row.end {
                            contained = true;
                        }
                        if row.start < origin + span && origin < row.end {
                            overlapped = true;
                        }
                    }
                } else {
                    for &w in &occ.partitions {
                        if (w as usize) >= n_windows {
                            continue;
                        }
                        for &row_index in
                            store.rows_of_partition_path(window_partition[w as usize], occ.path)
                        {
                            let row = &store.folds[row_index];
                            if row.start <= origin && origin + span <= row.end {
                                contained = true;
                            }
                            if row.start < origin + span && origin < row.end {
                                overlapped = true;
                            }
                        }
                    }
                }
                if contained {
                    data.in_axis += 1;
                } else if overlapped {
                    data.overhang += 1;
                    if locality_examples.len() < 8 {
                        locality_examples.push(json!({
                            "record": record, "occ": position,
                            "path": occ.path, "origin": origin, "span": span,
                            "class": "overhang",
                        }));
                    }
                } else {
                    data.extension += 1;
                    if locality_examples.len() < 8 {
                        locality_examples.push(json!({
                            "record": record, "occ": position,
                            "path": occ.path, "origin": origin, "span": span,
                            "class": "extension",
                        }));
                    }
                }
            }
            touching.insert(record, data);
        }
        // The evidence units: (record, variant).
        let mut units: Vec<Unit> = Vec::new();
        for &record in &locus_records {
            for (variant, v) in record_variants[&record].iter().enumerate() {
                units.push(Unit {
                    record,
                    variant,
                    mirror: v.mirror,
                    read: v.read,
                    w_lo: v.w_lo,
                    count: v.count,
                });
            }
        }
        let n_units = units.len();
        eprintln!(
            "[score] locus {locus}: {} records, {n_units} units \
             (touching/locality classification) [{:.1}s]",
            locus_records.len(),
            units_started.elapsed().as_secs_f64(),
        );

        // ------------------------------- the scoring matrix (parallel folds)
        let scoring_started = Instant::now();
        struct FoldScore {
            // (Candidate C: the row in its per-row smaller encoding —
            // the dense triple below is the pair loop's accumulator,
            // consumed by ScoreRow::from_pair_loop at the row's end.)
            row: ScoreRow,
            factorized_checked: u64,
            factorized_equal: u64,
            edge_discarded: u64,
            no_pin_records: u64,
            abstained_hull_bp: u64,
        }
        // (THE PER-RECORD HOIST — the marathon profile's named
        // fold-independent redundancy class: the canonical token
        // decode, the own anchor positions, and the merged
        // anchor-window hull are PURE functions of the record, but
        // the committed pair loop recomputed them at every (fold,
        // unit) pair — up to n_folds times per unit, up to ~2.4e9
        // recomputes at the receipted marathon worst loci, the
        // measured no-pin path's whole 0.2-0.3us/pair floor. Computed
        // ONCE per record here; the pair loop reads the table.
        // Byte-identity: the hoisted values are the same
        // pure-function outputs the pair loop computed; the scoring
        // math, the per-pair counter accumulation
        // (abstained_hull_bp adds the same hull_covered - backbone_bp
        // per pair), and every iteration order are unchanged.)
        struct RecordData {
            canonical: Vec<(i32, u64)>,
            own_flat: Vec<u64>,
            hull_covered: u64,
        }
        let mut hoist_max_record = 0usize;
        for unit in &units {
            hoist_max_record = hoist_max_record.max(unit.record);
        }
        let mut record_data: Vec<Option<RecordData>> =
            (0..=hoist_max_record).map(|_| None).collect();
        for unit in &units {
            let record = unit.record;
            if record_data[record].is_some() {
                continue;
            }
            let canonical = decode_record_tokens(
                &key_tokens[census_key[record].unwrap() as usize],
            )
            .expect("bound record");
            let own_flat = own_positions_of(&canonical, false, 0);
            let hull_windows: Vec<(u64, u64)> =
                own_flat.iter().map(|&p| (p, p + k)).collect();
            let hull_covered: u64 =
                merge_windows(hull_windows).iter().map(|&(lo, hi)| hi - lo).sum();
            record_data[record] = Some(RecordData {
                canonical,
                own_flat,
                hull_covered,
            });
        }
        // (The call-core census accumulators — per LOCUS, shared
        // atomically across the fold par_iter; only touched when the
        // env instrument is on, so the default path never pays an
        // atomic or an Instant.)
        let census_anchor_ns = AtomicU64::new(0);
        let census_placement_ns = AtomicU64::new(0);
        let census_assertion_ns = AtomicU64::new(0);
        let census_marginal_ns = AtomicU64::new(0);
        let census_anchor_calls = AtomicU64::new(0);
        let census_placement_calls = AtomicU64::new(0);
        let census_assertion_calls = AtomicU64::new(0);
        // (THE ANCHOR-CLASS HOIST's record slot: locus_records-position
        // per census record id, the O(1) unit→pins index.)
        let mut slot_of_record = vec![u32::MAX; census_records.len()];
        for (slot, &record) in locus_records.iter().enumerate() {
            slot_of_record[record] = slot as u32;
        }
        let scored: Vec<FoldScore> = (0..n_folds)
            .into_par_iter()
            .map(|fold_index| {
                let fold = &folds[fold_index];
                let mut ll = Vec::with_capacity(n_units);
                let mut pinned: Vec<bool> = Vec::with_capacity(n_units);
                let mut mismatches = Vec::with_capacity(n_units);
                let mut factorized_checked = 0u64;
                let mut factorized_equal = 0u64;
                let mut edge_discarded = 0u64;
                let mut no_pin_records = 0u64;
                let mut abstained_hull_bp = 0u64;
                let positions_of = |node: u32| fold.positions_of(node);
                // (THE ANCHOR-CLASS HOIST — the §10.6.14 autopsy's
                // re-anchorable COORDINATE-LOOKUP class, implemented as
                // the pure-function collapse it is: the pin build —
                // anchor_correspondence's monotone DP over the fold's
                // own node positions, the merged pinned-window
                // backbone, the windows-inside test — is a PURE
                // function of (record, fold, orientation): the
                // canonical walk, the touching orientations, and the
                // fold's node positions are all record- or fold-level
                // state and NO unit (variant) input enters it, so
                // every variant of one record shares one result per
                // fold. The committed form rebuilt it at every
                // (unit, fold, orientation) — 54.2M calls at L1,
                // 16.8 core-s of pure re-computation at 0.31 µs/call,
                // against ~150x fewer (record, fold, orientation)
                // builds. The hoist computes the pins ONCE per
                // (record, fold) before the unit loop; every per-unit
                // term — the abstained-hull accumulation (an integer
                // add, exact in any order), the no-pin elsewhere
                // branch, the per-(unit, pin) placement — stays in the
                // committed unit order with the committed values, so
                // the pair loop's outputs are bit-identical by
                // construction. The scoring-census instrument's
                // anchor timer/call counters now wrap the (record,
                // fold, orientation) builds — the class's true shape.)
                struct FoldPins {
                    pins: Vec<Option<Pin>>,
                    any_pin: bool,
                    // Per orientation IN THE COMMITTED ITERATION ORDER
                    // (the touching orientations' BTreeSet order) that
                    // pinned at all: its merged backbone bp — the
                    // per-unit abstained-hull term's own input.
                    abstain_backbone_bp: Vec<u64>,
                }
                let census_anchor_started =
                    if locus_scoring_census { Some(Instant::now()) } else { None };
                let mut pins_by_record: Vec<FoldPins> =
                    Vec::with_capacity(locus_records.len());
                for &record in &locus_records {
                    let data = record_data[record]
                        .as_ref()
                        .expect("bound record (the per-record hoist)");
                    let canonical = &data.canonical;
                    // (Stage 2: this fold's territory image — the
                    // placement-validity domain and the prior's
                    // support. Under the default convention it is the
                    // fold's full extent.)
                    let (img_lo, img_hi) = territory[fold_index];
                    let mut pins: Vec<Option<Pin>> = vec![None, None];
                    let mut any_pin = false;
                    let mut abstain_backbone_bp: Vec<u64> = Vec::new();
                    for &orientation in &touching[&record].orientations {
                        let (pinned, _ambiguous) =
                            anchor_correspondence(canonical, orientation, &positions_of);
                        if locus_scoring_census {
                            census_anchor_calls.fetch_add(1, Ordering::Relaxed);
                        }
                        if pinned.iter().all(|p| p.is_none()) {
                            continue;
                        }
                        // The merged pinned-window backbone (hull
                        // coordinates, mirror-invariant).
                        let own_flat = &data.own_flat;
                        let windows: Vec<(u64, u64)> = pinned
                            .iter()
                            .enumerate()
                            .filter_map(|(j, &p)| p.map(|r| (own_flat[j], r as u64)))
                            .map(|(p, _)| (p, p + k))
                            .collect();
                        let backbone_bp: u64 =
                            merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
                        let windows_inside = pinned.iter().flatten().all(|&r| {
                            r >= img_lo && r + k <= img_hi
                        });
                        // The abstained hull bp (the named remainder:
                        // read hull bases inside unpinned anchor
                        // windows) — the per-record hoisted value (the
                        // hull is fold-independent).
                        abstain_backbone_bp.push(backbone_bp);
                        if !windows_inside {
                            continue;
                        }
                        any_pin = true;
                        let pinned_i: Vec<Option<i64>> =
                            pinned.iter().map(|p| p.map(|v| v as i64)).collect();
                        pins[orientation as usize] = Some(Pin {
                            orientation,
                            pinned: pinned_i,
                            backbone_bp,
                        });
                    }
                    pins_by_record.push(FoldPins {
                        pins,
                        any_pin,
                        abstain_backbone_bp,
                    });
                }
                if let Some(started) = census_anchor_started {
                    census_anchor_ns
                        .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                for unit in &units {
                    let record = unit.record;
                    let data = record_data[record]
                        .as_ref()
                        .expect("bound record (the per-record hoist)");
                    let canonical = &data.canonical;
                    // (Stage 2: this fold's territory image — the
                    // placement-validity domain and the prior's
                    // support. Under the default convention it is the
                    // fold's full extent.)
                    let (img_lo, img_hi) = territory[fold_index];
                    // The record's memoized pins for THIS fold (the
                    // anchor-class hoist above): one build per
                    // (record, fold), shared by every variant.
                    // (The unit's record is one of the locus's own
                    // records — the committed code's touching[&record]
                    // lookup already assumed it; an out-of-set record
                    // fails loudly at the pins index.)
                    let slot = slot_of_record[record] as usize;
                    let built = &pins_by_record[slot];
                    for &backbone_bp in &built.abstain_backbone_bp {
                        abstained_hull_bp += data.hull_covered - backbone_bp;
                    }
                    if !built.any_pin {
                        // No anchored placement: the ELSEWHERE branch
                        // (the generated-elsewhere marginalization's
                        // floor; a read a candidate cannot place is
                        // explained by the genome outside the locale at
                        // the derived background rate, not by the
                        // all-mismatch minimum).
                        no_pin_records += 1;
                        ll.push(elsewhere);
                        pinned.push(false);
                        mismatches.push(0);
                        continue;
                    }
                    let mut scores: Vec<f64> = Vec::new();
                    let mut unit_mismatches = 0u32;
                    let mut pinned_here = false;
                    for pin in built.pins.iter().flatten() {
                        let census_placement_started =
                            if locus_scoring_census { Some(Instant::now()) } else { None };
                        let placement = placement_score(
                            unit, canonical, pin, &fold.seq, fold.len, &reads, k,
                            img_lo as i64, img_hi as i64,
                        );
                        if let Some(started) = census_placement_started {
                            census_placement_ns
                                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                            census_placement_calls.fetch_add(1, Ordering::Relaxed);
                        }
                        let Some((m, c)) = placement else {
                            edge_discarded += 1;
                            continue;
                        };
                        // (THE OWNER'S RULING 1: the direct_score
                        // exactness assertion is a receipt/gate/test
                        // instrument — IMPG_SCORE_EXACTNESS_ASSERT
                        // restores the committed always-on
                        // re-verification verbatim; the production
                        // default TRUSTS the factorized (m, c) here.
                        // No product/output bit ever depended on
                        // (dm, dc).)
                        if exactness_assert {
                            let census_assertion_started =
                                if locus_scoring_census { Some(Instant::now()) } else { None };
                            let (dm, dc) = direct_score(
                                unit, canonical, pin, &fold.seq, fold.len, &reads, k,
                            )
                            .expect("direct score");
                            if let Some(started) = census_assertion_started {
                                census_assertion_ns.fetch_add(
                                    started.elapsed().as_nanos() as u64,
                                    Ordering::Relaxed,
                                );
                                census_assertion_calls.fetch_add(1, Ordering::Relaxed);
                            }
                            factorized_checked += 1;
                            if (m, c) == (dm, dc) {
                                factorized_equal += 1;
                            }
                        }
                        unit_mismatches = unit_mismatches.max(c);
                        pinned_here = true;
                        scores.push(scoring.placement(m, c));
                    }
                    pinned.push(pinned_here);
                    // THE MARGINALIZATION: LL(read | fold) =
                    // logsumexp( the local-spell branch (the uniform
                    // placement prior + the anchored placements),
                    // E(read) — the candidate-independent
                    // generated-elsewhere branch ). A local branch
                    // that scores no better than E is absorbed by it
                    // (the read does not vote beyond E); the
                    // logsumexp2 is monotone in the local branch, so
                    // per-unit orderings the local branch preserves
                    // (a full-match truth placement vs a rival's
                    // worse placement) survive the marginalization.
                    let census_marginal_started =
                        if locus_scoring_census { Some(Instant::now()) } else { None };
                    ll.push(logsumexp2(
                        unit_log_likelihood(&scores, img_hi - img_lo, read_length),
                        elsewhere,
                    ));
                    mismatches.push(unit_mismatches);
                    if let Some(started) = census_marginal_started {
                        census_marginal_ns
                            .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    }
                }
                let census_row_started =
                    if locus_scoring_census { Some(Instant::now()) } else { None };
                let row =
                    ScoreRow::from_pair_loop(ll, pinned, mismatches, elsewhere);
                if let Some(started) = census_row_started {
                    census_marginal_ns
                        .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }
                FoldScore {
                    row,
                    factorized_checked,
                    factorized_equal,
                    edge_discarded,
                    no_pin_records,
                    abstained_hull_bp,
                }
            })
            .collect();
        let scoring_seconds = scoring_started.elapsed().as_secs_f64();
        // (Candidate C: the matrix is the per-row smaller encoding —
        // the committed dense clone-collect (a second 13B/pair copy,
        // transiently TWO copies of the whole matrix) is gone; the
        // pinned-majority rows keep the committed 13B/unit size and
        // the no-pin-majority rows collapse to 16B/pinned entry.)
        let factorized_checked: u64 = scored.iter().map(|s| s.factorized_checked).sum();
        let factorized_equal: u64 = scored.iter().map(|s| s.factorized_equal).sum();
        let edge_discarded: u64 = scored.iter().map(|s| s.edge_discarded).sum();
        let no_pin: u64 = scored.iter().map(|s| s.no_pin_records).sum();
        let abstained_hull_bp: u64 = scored.iter().map(|s| s.abstained_hull_bp).sum();
        let matrix = ScoreMatrix {
            n_units,
            elsewhere,
            rows: scored.into_iter().map(|s| s.row).collect(),
            unit_fold: None,
        };
        eprintln!(
            "[score] locus {locus}: matrix {} folds x {} units = {} pairs, \
             {} live MB (row-compressed; the dense form is {} MB)",
            matrix.n_folds(),
            n_units,
            matrix.n_folds() as u64 * n_units as u64,
            matrix.live_bytes() / (1024 * 1024),
            matrix.n_folds() as u64 * n_units as u64 * 13 / (1024 * 1024),
        );
        ensure(
            factorized_checked == factorized_equal,
            &format!(
                "FACTORIZATION GATE FAILED at locus {locus}: \
                 {factorized_equal} of {factorized_checked} equal"
            ),
        )?;
        eprintln!(
            "[score] locus {locus}: factorized {factorized_equal}/{factorized_checked} exact, \
             edge-discarded placements {edge_discarded}, no-pin (unit,fold) {no_pin}, \
             abstained hull bp {abstained_hull_bp} [{scoring_seconds:.1}s]",
        );
        rss.probe(&format!("locus_{locus}_scored"))?;

        // --------------------------- the observed census mass (moved up
        // from below the truth folds: pure code motion — the merge pass
        // reads it; the values are identical to the committed site's,
        // and every later use is unchanged.)
        let observed_nodes: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for &record in &locus_records {
                let multiplicity = census_records[record].multiplicity as f64;
                for &node in &touching[&record].covered_nodes {
                    *map.entry(node).or_default() += multiplicity;
                }
            }
            map
        };
        let observed_edges: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for &record in &locus_records {
                let multiplicity = census_records[record].multiplicity as f64;
                for &edge in &touching[&record].covered_edges {
                    *map.entry(edge).or_default() += multiplicity;
                }
            }
            map
        };
        let observed_node_mass: f64 = observed_nodes.values().sum();
        let observed_edge_mass: f64 = observed_edges.values().sum();

        // --------------- the alignment-connected fold relaxation pass
        // (env IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS; the derived rule
        // documented at the merge function. Env-unset: the exact-rule
        // structures pass through untouched — the committed receipts
        // byte-identical. Env-set: the exact folds merge by alignment
        // connection subject to the read-level indistinguishability
        // safety property, and every fold-indexed structure regroups
        // onto the relaxed space (the representative's rows — asserted
        // bit-identical to every member's).)
        let merge_axis_fold = territory_axis_fold.or_else(|| {
            folds.iter().position(|fold| {
                fold.members.iter().any(|m| {
                    m.path_name == panel_component
                        && m.start == axis_start
                        && m.end == axis_end
                })
            })
        });
        // (The phase census's merge timer — 0.0 and untouched with the
        // census env unset or the merge env unset.)
        let mut merge_seconds = 0.0f64;
        let (
            folds,
            mut matrix,
            territory,
            territory_axis_fold,
            merge_groups,
            merge_counters,
        ) = if alignment_connected {
            let merge_started = locus_phase_census.then(Instant::now);
            let merged = alignment_connected_merge(
                folds,
                matrix,
                territory,
                merge_axis_fold,
                territory_axis_fold,
                &observed_nodes,
                &observed_edges,
            );
            merge_seconds = merge_started
                .map(|started| started.elapsed().as_secs_f64())
                .unwrap_or(0.0);
            merged
        } else {
            (
                folds,
                matrix,
                territory,
                territory_axis_fold,
                Vec::new(),
                MergeCounters {
                    equal_column_pairs: 0,
                    refused_outside_bubble: 0,
                    refused_unshared: 0,
                    refused_pocket_mass: 0,
                    merged_pairs: 0,
                },
            )
        };
        let n_folds = folds.len();
        if alignment_connected && !merge_groups.is_empty() {
            eprintln!(
                "[score] locus {locus}: alignment-connected folds: {} merge groups \
                 ({} admitted pairs of {} equal-column candidates; refusals: \
                 {} outside-bubble, {} unshared, {} pocket-mass)",
                merge_groups.len(),
                merge_counters.merged_pairs,
                merge_counters.equal_column_pairs,
                merge_counters.refused_outside_bubble,
                merge_counters.refused_unshared,
                merge_counters.refused_pocket_mass,
            );
        }

        // ------------------------------------------- the exhaustive classes
        // (THE OWNER'S RULING 3 at the scan itself: EXACT = the
        // committed unit-order accumulation verbatim (the toggle
        // doctrine — IMPG_CLASS_SCAN_ORDER=exact), bit-reproducible;
        // SPARSE-SKIP = the legacy fast mode — the union of the two
        // rows' pinned units plus one constant-mass term (superseded
        // by the bounded default; unchanged). The compare receipt
        // (IMPG_CLASS_SCAN_COMPARE) sums BOTH ways and uses the EXACT
        // value downstream, so it changes no output byte.
        // THE SANCTIONED BOUNDED CLASS SCAN (SS10.6.21, the owner's
        // "take the better, faster one") runs INSTEAD of the
        // exhaustive/sparse forms at the production default: STAGE
        // ONE here — the deterministic branch-and-bound search for
        // best over the per-fold excess-mass bound, the strict
        // running-best form (a pair whose bound falls below the
        // running best can neither win nor tie the emitted argmax);
        // STAGE TWO below, once the truth pair is resolved — the
        // kept set {bound >= best + the derived exp-underflow
        // point} (every skipped pair's QUAL terms are exactly +0.0,
        // the committed bit-identity argument) extended by the truth
        // window where the truth pair sits outside it (the receipt's
        // truth fields need every pair that can be strictly above
        // or tied with the truth pair). Every computed value is the
        // exhaustive scan's own — the committed unit-order
        // accumulation, bit for bit — so the ANSWERS (the MAP pair,
        // the called set, the QUAL terms, the truth fields, the
        // best-outside rival) are the exhaustive run's by the
        // bound-plus-slack's own arithmetic and the final-value
        // tie-set rebuild (SS10.6.22). The skipped pairs' instrument
        // entries become the labeled bounded receipt (the null
        // skip certificate + the per-pair bounds and their derived
        // slack); the compare receipt has no bounded form (its
        // both-ways sum IS the bounded scan's own bookkeeping)
        // and is ignored here.)
        let class_started = Instant::now();
        let pair_count = n_folds * (n_folds + 1) / 2;
        let unit_counts: Vec<f64> = units.iter().map(|u| u.count as f64).collect();
        // THE UNIT-FOLD REDUCTION (env IMPG_UNIT_FOLD — see
        // `UnitFoldPlan`): the plan is attached to the matrix ONCE, on
        // the FINAL (post-merge) matrix, so every exact-sum site below
        // — the bounded scan's two stages, its best-outside fallback
        // and the exhaustive arm — walks the same surviving groups and
        // adds the same constant scalar; there is no site that can keep
        // the old walk by omission.
        if unit_fold_enabled {
            match matrix.fold_units(&unit_counts) {
                Some(plan) => eprintln!(
                    "[unit-fold] locus {locus}: U {} -> {} columns \
                     ({} constant, {} merged), pair terms {} -> {} \
                     [constant units {}, non-constant units {}, \
                     build {:.3}s]",
                    plan.units,
                    plan.distinct_columns,
                    plan.constant_columns,
                    plan.merged_columns(),
                    plan.units,
                    plan.pair_terms,
                    plan.constant_units,
                    plan.units - plan.constant_units,
                    plan.build_seconds,
                ),
                None => eprintln!(
                    "[unit-fold] locus {locus}: declined (one fold \
                     — every column is trivially constant, so the whole \
                     sum would be one product)"
                ),
            }
        }
        let sparse_skip = matches!(class_scan_order, ClassScanOrder::SparseSkip);
        // THE BIT-EQUAL-ROW MEMOIZATION (the fourth expressibility
        // reduction — see `ClassScanMemo`): built once per locus on
        // the FINAL (post-merge) matrix, before every exact-sum site
        // — the exhaustive scan, the bounded scan's two stages and
        // its best-outside fallback all share one memo table per
        // locus. IMPG_CLASS_SCAN_MEMO=off restores the direct
        // per-pair form verbatim (the toggle doctrine's measurement
        // restore path); the memo changes no emitted byte either
        // way, so no output label is warranted.
        let mut scan_memo: Option<ClassScanMemo> = if scan_memo_enabled {
            Some(ClassScanMemo::build(&matrix))
        } else {
            None
        };
        // (The compare receipt's per-locus atomics: class pairs
        // summed both ways, pairs whose f64 bits differ, the max ulp
        // delta, and the constant (fold, unit) entries the fast
        // form skips. Declared at the scan's own scope — the
        // sparse/compare receipt line below reads them after the
        // branch; under the bounded mode they stay zero and the
        // line does not print.)
        let compare_pairs = AtomicU64::new(0);
        let compare_diff = AtomicU64::new(0);
        let compare_max_ulp = AtomicU64::new(0);
        let compare_skipped = AtomicU64::new(0);
        // THE PINNED-UNION OPTION FOR THE BOUNDED PATH (the production
        // default's own scan). The tail locus is 96.9% constant, so its
        // pair sums need ~2,700 corrections rather than 45,423 — 11x less
        // EXACT work (measured: T00048 l8, 3.6e8 vs 4.2e9 unit-pairs,
        // 17x on that locus alone). The union walk reassociates the float
        // sum, so with IMPG_CLASS_SCAN_UNION it IS the value (the decision
        // gate is calls.jsonl), and with IMPG_CLASS_SCAN_COMPARE both forms
        // are computed and the delta is receipted while the exact value
        // is kept (no answer change).
        let bounded_union = std::env::var("IMPG_CLASS_SCAN_UNION").is_ok();
        // The per-locus bounded state, stage one -> stage two.
        let mut bounded_state: Option<BoundedScanState> = None;
        let (mut class_lls, mut best_ll, mut called): (Vec<f64>, f64, Vec<(usize, usize)>) =
            if class_scan_bounded {
                let (pair_bounds, pair_slack, computed, best, memo_values) =
                    bounded_stage_one(&matrix, &unit_counts, scan_memo.as_ref());
                let stage_one_seconds = class_started.elapsed().as_secs_f64();
                // The PROVISIONAL called set from stage one's computed
                // prefix, in flat order (the exhaustive form's own
                // order) — it carries the winner and the scan-site
                // invariants (best finite, called non-empty) until
                // stage two's kept window completes the exact values;
                // the emitted set is REBUILT from the final class_lls
                // after stage two (the T03226 locus-16 fix — never a
                // function of the bound's rounding).
                let called: Vec<(usize, usize)> = computed
                    .iter()
                    .enumerate()
                    .filter_map(|(flat, value)| {
                        let sum = (*value)?;
                        (sum == best).then(|| unrank_pair(flat))
                    })
                    .collect();
                let class_lls: Vec<f64> = computed
                    .iter()
                    .map(|value| value.unwrap_or(f64::NEG_INFINITY))
                    .collect();
                bounded_state = Some(BoundedScanState {
                    pair_bounds,
                    pair_slack,
                    computed,
                    best,
                    stage_one_seconds,
                    scan_memo: scan_memo.take(),
                    memo_values,
                });
                (class_lls, best, called)
            } else {
                let scan_fast = if sparse_skip || class_scan_compare {
                    Some(ClassScanFastIndex::build(&matrix, &unit_counts))
                } else {
                    None
                };
                // THE MEMOIZED EXHAUSTIVE SCAN (the exact arm — the
                // assessment form's own): phase A computes each
                // DISTINCT canonical pair sum exactly once, in
                // parallel, by the committed unit-order accumulation
                // on the representative rows; phase B fills every
                // pair's own flat slot from the memo — the values are
                // the direct scan's own bits (the construction proof
                // at `ClassScanMemo`), the flat order untouched. The
                // legacy sparse-skip/compare arms keep the direct
                // per-pair form (they are measurement forms, off at
                // the defaults).
                let class_lls: Vec<f64> =
                    if let (Some(memo), None) = (scan_memo.as_ref(), scan_fast.as_ref()) {
                        let distinct: Vec<usize> = (0..pair_count)
                            .into_par_iter()
                            .filter(|&flat| {
                                let (first, second) = unrank_pair(flat);
                                memo.canonical_flat(first, second) == flat
                            })
                            .collect();
                        memo.note_computed(distinct.len() as u64);
                        memo.note_hits(pair_count as u64 - distinct.len() as u64);
                        let sums: Vec<(usize, f64)> = distinct
                            .into_par_iter()
                            .map(|flat| {
                                let (first, second) = unrank_pair(flat);
                                (flat, matrix.class_pair_sum(first, second, &unit_counts))
                            })
                            .collect();
                        let mut values = vec![f64::NAN; pair_count];
                        for (flat, sum) in sums {
                            values[flat] = sum;
                        }
                        (0..pair_count)
                            .into_par_iter()
                            .map(|flat| {
                                let (first, second) = unrank_pair(flat);
                                values[memo.canonical_flat(first, second)]
                            })
                            .collect()
                    } else {
                        (0..pair_count)
                            .into_par_iter()
                            .map(|flat| {
                                let (first, second) = unrank_pair(flat);
                                match scan_fast.as_ref() {
                                    None => matrix.class_pair_sum(first, second, &unit_counts),
                                    Some(fast) => {
                                        let fast_sum =
                                            fast.class_pair_sum(first, second, &unit_counts);
                                        if !class_scan_compare {
                                            return fast_sum;
                                        }
                                        let exact =
                                            matrix.class_pair_sum(first, second, &unit_counts);
                                        compare_pairs.fetch_add(1, Ordering::Relaxed);
                                        compare_skipped.fetch_add(
                                            fast.constant_units_skipped(first, second),
                                            Ordering::Relaxed,
                                        );
                                        if exact.to_bits() != fast_sum.to_bits() {
                                            compare_diff.fetch_add(1, Ordering::Relaxed);
                                            let ulp = ulp_distance(exact, fast_sum);
                                            let mut current =
                                                compare_max_ulp.load(Ordering::Relaxed);
                                            while ulp > current {
                                                match compare_max_ulp.compare_exchange_weak(
                                                    current,
                                                    ulp,
                                                    Ordering::Relaxed,
                                                    Ordering::Relaxed,
                                                ) {
                                                    Ok(_) => break,
                                                    Err(seen) => current = seen,
                                                }
                                            }
                                        }
                                        if sparse_skip {
                                            fast_sum
                                        } else {
                                            exact
                                        }
                                    }
                                }
                            })
                            .collect()
                    };
                let (best_ll, called) = called_set_from_class_lls(&class_lls);
                (class_lls, best_ll, called)
            };
        let mut class_seconds = class_started.elapsed().as_secs_f64();
        if locus_phase_census && !class_scan_bounded {
            eprintln!(
                "[score] locus {locus}: class scan {pair_count} pairs over {n_units} units \
                 [{class_seconds:.1}s] (phase census)"
            );
        }
        if !class_scan_bounded && (sparse_skip || class_scan_compare) {
            eprintln!(
                "[score] locus {locus}: class scan {} pairs in {class_seconds:.1}s \
                 ({}{}; {} of {} pairs' sums differ bitwise, max {} ulp; \
                 {} constant unit-pairs skipped)",
                pair_count,
                if sparse_skip { "sparse-skip" } else { "exact" },
                if class_scan_compare { " + compare" } else { "" },
                compare_diff.load(Ordering::Relaxed),
                compare_pairs.load(Ordering::Relaxed),
                compare_max_ulp.load(Ordering::Relaxed),
                compare_skipped.load(Ordering::Relaxed),
            );
        }
        if class_scan_bounded && (bounded_union || class_scan_compare) {
            eprintln!(
                "[score] locus {locus}: class scan bounded{}; {} of {} kept-pair sums \
                 differ bitwise from the exact form, max {} ulp; {} constant unit-pairs skipped",
                if bounded_union { " + pinned-union" } else { " + compare" },
                UNION_CENSUS_DIFF.load(Ordering::Relaxed),
                UNION_CENSUS_PAIRS.load(Ordering::Relaxed),
                UNION_CENSUS_MAXULP.load(Ordering::Relaxed),
                UNION_CENSUS_SKIPPED.load(Ordering::Relaxed),
            );
        }
        ensure(best_ll.is_finite(), "no class has anchored evidence at this locus")?;
        let mut winner_pair = called[0];

        // -------------------------------------------------- the truth folds
        // (Window-aware, the chrIV generalization: the truth FIRST
        // fold is the fold carrying THIS WINDOW'S OWN axis row —
        // path + interval — because a repeat-locality partition can
        // host several of the component's windows, each in its own
        // fold, and the path-name rule would be ambiguous there. The
        // truth SECOND fold keeps the committed path-name rule, but
        // only at partitions hosting a SINGLE window: at a
        // multi-window partition an SK1 member row cannot be
        // attributed to this window by placement alone (measured:
        // chrIV partition 110 carries two SK1 rows that are OTHER
        // windows' orthologs), so no truth pair is fabricated — the
        // locus is scored with the truth pair inexpressible.)
        let truth_first = {
            let hits: Vec<usize> = folds
                .iter()
                .enumerate()
                .filter(|(_, fold)| {
                    fold.members.iter().any(|m| {
                        m.path_name == panel_component
                            && m.start == axis_start
                            && m.end == axis_end
                    })
                })
                .map(|(index, _)| index)
                .collect();
            ensure(hits.len() <= 1, "multiple folds carry the window's axis row")?;
            hits.first().copied()
        };
        let (truth_second, truth_second_correspondence) = if alignment_connected {
            // (The relaxation's attribution clause: the truth-second
            // fold is the fold carrying the haplotype's rows with the
            // MAXIMUM shared-node mass with the axis fold — the
            // alignment connection itself, replacing the exact rule's
            // single-window/exactly-one-fold refusal; ties and
            // disconnected candidates still refuse — no fabrication.)
            connected_truth_second(
                &folds,
                truth_first,
                &format!("SK1#0#{contig}"),
            )
        } else if windows_of_partition[&partition] > 1 {
            (None, None)
        } else {
            // (THE FLEET GENERALIZATION, measured at chrXIV/chrXVI/
            // chrXVII-class repeat-locality partitions: a SINGLE-window
            // partition can still carry SEVERAL folds with SK1 member
            // rows (chrXIV partition 16 carries SK1 [99039,104912) and
            // [558497,563620); chrXVI partitions 117/257 similar) — the
            // path-name rule cannot attribute WHICH SK1 fold is this
            // window's ortholog, so NO truth pair is fabricated, the
            // same no-fabrication convention as the multi-window
            // partitions. The committed single-fold case is unchanged;
            // chrMT/chrI/chrIV never hit the ambiguous case (their
            // single-window partitions carry at most one SK1 fold —
            // the identity gates prove the change is a no-op there.)
            let hits: Vec<usize> = folds
                .iter()
                .enumerate()
                .filter(|(_, fold)| {
                    fold.members
                        .iter()
                        .any(|m| m.path_name == format!("SK1#0#{contig}"))
                })
                .map(|(index, _)| index)
                .collect();
            (
                if hits.len() == 1 {
                    hits.first().copied()
                } else {
                    None
                },
                None,
            )
        };
        let truth_pair = match (truth_first, truth_second) {
            (Some(a), Some(b)) => Some((a.min(b), a.max(b))),
            _ => None,
        };
        // THE TEST-MODE TRUTH FOLDS (the --truth-qv input made
        // explicit): the first truth haplotype must be the component
        // lane itself — the window's own axis row (the committed
        // convention's truth-first rule, verbatim); the second is the
        // input haplotype's fold under the committed path-name rule
        // (single-window partitions only; exactly one fold). The
        // committed instrument hard-coded the validation sample's
        // SK1 second haplotype; the test mode requires it as input.
        // (The alignment-connected relaxation's test-mode
        // correspondence receipt, emitted in the truth-QV artifact.)
        let mut truth_second_correspondence_test: Option<serde_json::Value> = None;
        let truth_pair_test = if truth_qv_mode {
            ensure(
                truth_qv_haplotypes[0] == panel_component,
                "the first --truth-qv haplotype must be the component lane \
                 (the window's own axis row is the truth-first fold)",
            )?;
            let second_name = truth_qv_haplotypes[1].as_str();
            let (second, correspondence) = if alignment_connected {
                // (The relaxation's attribution clause, test mode: the
                // same alignment-connected correspondence as the
                // committed truth pair's — the four segdup windows'
                // expressibility clause.)
                connected_truth_second(&folds, truth_first, second_name)
            } else if windows_of_partition[&partition] > 1 {
                (None, None)
            } else {
                let hits: Vec<usize> = folds
                    .iter()
                    .enumerate()
                    .filter(|(_, fold)| {
                        fold.members
                            .iter()
                            .any(|m| m.path_name == second_name)
                    })
                    .map(|(index, _)| index)
                    .collect();
                ((hits.len() == 1).then(|| hits[0]), None)
            };
            truth_second_correspondence_test = correspondence;
            match (truth_first, second) {
                (Some(a), Some(b)) => Some((a.min(b), a.max(b))),
                _ => None,
            }
        } else {
            None
        };
        // The test-mode pair's LL rank (the committed receipts'
        // truth_rank/rank1 fields, for the same convention): 1 + the
        // count of classes strictly above the pair's class LL.
        let truth_rank_test = truth_pair_test.map(|(a, b)| {
            let flat = b * (b + 1) / 2 + a;
            // (the union path: the truth pair may sit outside the
            // certification set, so its own value is recomputed in the
            // committed order — the rank is emitted, so it must be exact)
            let value = if std::env::var("IMPG_CLASS_SCAN_UNION").is_ok()
                && !std::env::var("IMPG_CLASS_SCAN_COMPARE").is_ok()
            {
                matrix.class_pair_sum(a, b, &unit_counts)
            } else {
                class_lls[flat]
            };
            1 + class_lls.iter().filter(|&&ll| ll > value).count()
        });
        let truth_flat = truth_pair.map(|(a, b)| b * (b + 1) / 2 + a);
        // ------------------------- the bounded scan's STAGE TWO (the kept set)
        // (The truth pair is now resolved, so the derived keep
        // threshold is final: the QUAL window {bound + the derived
        // rounding-error slack >= best + the exp-underflow point}
        // extended by the truth window where the truth pair sits
        // outside it. The truth pair itself is computed FIRST — it
        // is always kept (its own bound + slack >= its exact ll) —
        // so the receipt's truth fields read the exact value.
        // Every pair the stage-two pass computes is summed by the
        // committed unit-order accumulation, so every value is the
        // exhaustive scan's own bit for bit; the kept/skipped sets
        // are deterministic functions of the matrix alone.)
        if let Some(state) = bounded_state.as_mut() {
            let stage_two_started = Instant::now();
            // THE MEMOIZED EXACT SUM (the fourth reduction at the
            // bounded path's own sites): the pair's canonical flat's
            // value — computed on the representative rows if the
            // slot is still empty, read from the memo otherwise —
            // bit-identical to the direct per-pair accumulation by
            // construction (`ClassScanMemo`'s proof).
            let memoized_sum = |state: &mut BoundedScanState,
                                first: usize,
                                second: usize|
             -> f64 {
                match state.scan_memo.as_ref() {
                    None => matrix.class_pair_sum(first, second, &unit_counts),
                    Some(memo) => {
                        let canon = memo.canonical_flat(first, second);
                        match state.memo_values[canon] {
                            Some(value) => {
                                memo.note_hits(1);
                                value
                            }
                            None => {
                                let (r_first, r_second) = unrank_pair(canon);
                                let value =
                                    matrix.class_pair_sum(r_first, r_second, &unit_counts);
                                memo.note_computed(1);
                                state.memo_values[canon] = Some(value);
                                value
                            }
                        }
                    }
                }
            };
            let truth_ll: Option<f64> = match truth_flat {
                Some(flat) => match state.computed[flat] {
                    Some(value) => Some(value),
                    None => {
                        let (first, second) = unrank_pair(flat);
                        let value = memoized_sum(state, first, second);
                        state.computed[flat] = Some(value);
                        class_lls[flat] = value;
                        Some(value)
                    }
                },
                None => None,
            };
            let keep_threshold = bounded_keep_threshold(state.best, truth_ll);
            let kept: Vec<usize> = (0..pair_count)
                .filter(|&flat| {
                    state.computed[flat].is_none()
                        && state.pair_bounds[flat] + state.pair_slack[flat] >= keep_threshold
                })
                .collect();
            // the kept set's distinct canonical sums first (each
            // computed exactly once, in parallel), then every kept
            // pair reads its own flat slot — the values shared,
            // never reordered
            let sums: Vec<(usize, f64)> = match state.scan_memo.as_ref() {
                None => kept
                    .par_iter()
                    .map(|&flat| {
                        let (first, second) = unrank_pair(flat);
                        (flat, matrix.class_pair_sum(first, second, &unit_counts))
                    })
                    .collect(),
                Some(memo) => {
                    let keyed: Vec<(usize, usize)> = kept
                        .iter()
                        .map(|&flat| {
                            let (first, second) = unrank_pair(flat);
                            (flat, memo.canonical_flat(first, second))
                        })
                        .collect();
                    let mut missing: Vec<usize> =
                        keyed.iter().map(|&(_, canon)| canon).collect();
                    missing.sort_unstable();
                    missing.dedup();
                    let owned: Vec<usize> = missing
                        .iter()
                        .copied()
                        .filter(|&canon| state.memo_values[canon].is_none())
                        .collect();
                    memo.note_computed(owned.len() as u64);
                    memo.note_hits(keyed.len() as u64 - owned.len() as u64);
                    let fresh: Vec<(usize, f64)> = owned
                        .par_iter()
                        .map(|&canon| {
                            let (first, second) = unrank_pair(canon);
                            (canon, matrix.class_pair_sum(first, second, &unit_counts))
                        })
                        .collect();
                    for (canon, sum) in fresh {
                        state.memo_values[canon] = Some(sum);
                    }
                    keyed
                        .into_iter()
                        .map(|(flat, canon)| {
                            (
                                flat,
                                state.memo_values[canon]
                                    .expect("the kept set's canonical sum is filled"),
                            )
                        })
                        .collect()
                }
            };
            for (flat, sum) in sums {
                state.computed[flat] = Some(sum);
                class_lls[flat] = sum;
            }
            // THE CERTIFICATION (the pinned-union path only): the
            // decisive values are recomputed in the COMMITTED order, so
            // every emitted byte that depends on a class LL is the exact
            // form's. The set is: the winner's float-window (a pair whose
            // reassociated value lands within it could flip the called set
            // or the rank) plus the TOP 32 by fast value (the alternative
            // and the tie breakers come from the top of the spectrum).
            // Measured at the tail locus: ONE pair falls in the window.
            if std::env::var("IMPG_CLASS_SCAN_UNION").is_ok() && !std::env::var("IMPG_CLASS_SCAN_COMPARE").is_ok() {
                let n_units = matrix.n_units as f64;
                let best_fast = class_lls
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let window = 2.0 * n_units * f64::EPSILON * best_fast.abs().max(1.0);
                let mut order: Vec<usize> = (0..class_lls.len()).collect();
                order.sort_unstable_by(|&a, &b| {
                    class_lls[b].partial_cmp(&class_lls[a]).expect("finite class LLs")
                });
                let mut rechecks = 0usize;
                for &flat in order.iter().take(32) {
                    let (first, second) = unrank_pair(flat);
                    class_lls[flat] = matrix.class_pair_sum(first, second, &unit_counts);
                    rechecks += 1;
                }
                for &flat in order.iter() {
                    if class_lls[flat] < best_fast - window {
                        break;
                    }
                    if rechecks > 0 && order.iter().take(32).any(|&f| f == flat) {
                        continue;
                    }
                    let (first, second) = unrank_pair(flat);
                    class_lls[flat] = matrix.class_pair_sum(first, second, &unit_counts);
                    rechecks += 1;
                }
                UNION_CENSUS_RECHECKS.fetch_add(rechecks as u64, Ordering::Relaxed);
            }
            // THE TIE-SET REBUILD (the T03226 locus-16 fix,
            // SS10.6.22): the emitted called set, best, and the
            // winner they name are REBUILT from the FINAL class_lls
            // — a pure function of the final exact values, never of
            // the bound's rounding. Stage two's kept window has
            // already recovered every tie bit-identical (a tie's
            // bound + slack sits at best, far above the keep
            // threshold), so this equals the exhaustive called set;
            // the best-outside fallback that runs later can only
            // patch pairs whose exact sums sit below the keep
            // threshold — far below best — so it can never change
            // the rebuilt set.
            let (rebuilt_best, rebuilt_called) = called_set_from_class_lls(&class_lls);
            best_ll = rebuilt_best;
            called = rebuilt_called;
            winner_pair = called[0];
            class_seconds += stage_two_started.elapsed().as_secs_f64();
            let known = state
                .computed
                .iter()
                .filter(|value| value.is_some())
                .count();
            let skipped = pair_count - known;
            eprintln!(
                "[score] locus {locus}: class scan bounded: {pair_count} pairs in \
                 {class_seconds:.1}s (stage one {:.1}s; kept {known}, skipped \
                 {skipped} ({:.2}% pruned); the derived keep threshold and the \
                 per-pair bounds are the receipt's skip certificate)",
                state.stage_one_seconds,
                100.0 * skipped as f64 / pair_count as f64,
            );
        }
        let truth_log_likelihood = truth_flat.map(|flat| class_lls[flat]);
        let truth_rank = truth_log_likelihood.map(|value| {
            1 + class_lls
                .iter()
                .filter(|&&ll| ll > value)
                .count()
        });
        let truth_tied_classes = truth_log_likelihood.map(|value| {
            class_lls
                .iter()
                .enumerate()
                .filter(|&(flat, &ll)| ll == value && Some(flat) != truth_flat)
                .count()
        });
        let log_gap = match (best_ll, truth_log_likelihood) {
            (best, Some(truth)) => Some(best - truth),
            _ => None,
        };

        // ------------------------------------------- observed mass (no shares)
        // (Moved above the fold-relaxation pass — the merge reads it;
        // the committed values are reused verbatim below.)
        let observed_node_fn = |key: u64| observed_nodes.get(&key).copied().unwrap_or(0.0);
        let observed_edge_fn = |key: u64| observed_edges.get(&key).copied().unwrap_or(0.0);

        // ------------------------------------------------- the QUAL cluster
        let qual_started = Instant::now();
        let winner_nodes = {
            let (first, second) = winner_pair;
            merged_multiset(&folds[first].nodes, &folds[second].nodes)
        };
        let winner_edges = {
            let (first, second) = winner_pair;
            merged_multiset(&folds[first].edges, &folds[second].edges)
        };
        let called_score = Some(best_ll);
        let spectrum_started = locus_phase_census.then(Instant::now);
        let spectrum: Vec<(f64, f64, f64, usize, usize, usize)> =
            match (winner_pair, called_score) {
                ((winner_first, winner_second), Some(score)) => {
                    let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
                        .into_par_iter()
                        .filter_map(|flat| {
                            let (first, second) = unrank_pair(flat);
                            if first == winner_first && second == winner_second {
                                return None;
                            }
                            let ll = class_lls[flat];
                            let nodes =
                                merged_multiset(&folds[first].nodes, &folds[second].nodes);
                            let edges =
                                merged_multiset(&folds[first].edges, &folds[second].edges);
                            let distance = differing_observed_mass(&nodes, &winner_nodes, &observed_node_fn)
                                + differing_observed_mass(&edges, &winner_edges, &observed_edge_fn);
                            let signature =
                                signature_cosine_distance(&nodes, &edges, &winner_nodes, &winner_edges);
                            Some((distance, (ll - score).exp(), signature, first, second, flat))
                        })
                        .collect();
                    spectrum
                        .sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
                    spectrum
                }
                _ => Vec::new(),
            };
        // (Stage 2, THE WALL FIX — the degenerate-locus tie
        // certificate. At a fully degenerate locus every fold's every
        // unit scores the candidate-independent E branch (every
        // territory image is empty or shorter than a read, so no
        // candidate places anything), every one of the pair_count
        // classes ties at the pure-E score and called.len() ==
        // pair_count — an exact structural equality, not a size
        // threshold. The naive certificate materializes |called| - 1
        // bands x |spectrum| observed-mass distances — measured at the
        // repeat-domain monster loci (the partition-110 windows
        // chrXII L86/chrXV L14/chrX L52/chrVII L56): 35,510 x 35,510 =
        // 1.26e9 multiset merges, qual_seconds 4,802s in the serial
        // run, and a ~6GB qual_tied_bands DOM that breached the
        // 64GiB guard in every 4-wide that reached it. The cluster
        // machinery consumes from each band ONLY its <=cut predicate
        // set, so the degenerate path computes exactly that: each
        // (tied class, spectrum entry) distance is decided by SOUND
        // two-sided O(1) bounds — D = Dn + De >= |Sm_n(t) - Sm_n(e)| +
        // |Sm_e(t) - Sm_e(e)| and D <= Sm_n(t) + Sm_n(e) + Sm_e(t) +
        // Sm_e(e), where Sm is the observed mass (the bounds are
        // exact-deciding: the upper bound <= cut proves the predicate
        // TRUE, the lower bound > cut proves it FALSE) — with the
        // exact differing_observed_mass as the fallback for the
        // undecided middle band. The emitted certificate is the
        // compact cut-index list; the answers (k, excluded,
        // alternative, cluster_size, shape, p) are identical to the
        // naive path BY CONSTRUCTION and unit-gated against a
        // verbatim oracle below.
        //
        // THE PER-PAIR CONSTANT'S FIX (the fold-merge open term,
        // closed): the NON-degenerate path took the SAME naive band
        // form (|called| - 1 bands x |spectrum| entries, every entry
        // a fresh merged-multiset pair plus two observed-mass walks)
        // — the ~6.4 ms-per-class-pair constant the PRE-MERGE
        // fold-heavy shapes receipt (the ledger's T03226 l15: 228.4s
        // serial pre-merge vs 0.5s committed, the tie class the
        // committed merge collapses 455x). The cluster machinery's
        // consumption contract is unchanged (only each band's <=cut
        // predicate set), so BOTH paths now take the cut-index-list
        // certificate — per-fold mass signatures computed once, O(1)
        // two-sided bound decisions per (tied class, spectrum entry),
        // the exact walk only in the undecided middle band — and the
        // answers are the naive band matrix's own at every shape
        // (the non-degenerate oracle test below proves it). The
        // emitted receipt member changes at NON-degenerate tie loci:
        // qual_tied_bands (the full band vectors) is replaced by
        // qual_tied_cut_indices (the compact predicate lists) — the
        // degenerate precedent's own form, sanctioned with the fix.)
        let degenerate = called.len() == pair_count;
        let mut cut_index_lists: Vec<(usize, CutSet)> = Vec::new();
        // THE O(f^4) IS THE ENUMERATION, NOT THE ANSWER. With
        // IMPG_TIED_FOLD_FAST the certificate's per-class predicate sets are
        // never materialized: `tie_cluster_fast` answers k (connected
        // components of the d<=cut relation) and the excluded union from the
        // same predicates, so the lists have no consumer left. They are
        // skipped here, and the report member `qual_tied_cut_indices` —
        // which was the certificate's own enumeration — is not produced by
        // that path. Everything the answer depends on (k, excluded,
        // alternative, p, cluster_size) is computed identically, which is
        // what the band-oracle tests compare.
        let tied_fold_fast = std::env::var("IMPG_TIED_FOLD_FAST").is_ok();
        let certificate_started = locus_phase_census.then(Instant::now);
        if !tied_fold_fast {
            let mut sorted: Vec<f64> = spectrum
                .iter()
                .map(|&(distance, ..)| distance)
                .collect();
            sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
            // (The same knee the cluster form itself takes — the
            // predicate sets are decided at exactly the cut the
            // cluster consumes them at.)
            let cut = spectrum_knee(&sorted).cut;
            cut_index_lists = tied_class_cut_lists(
                &folds,
                &called,
                &spectrum,
                &observed_node_fn,
                &observed_edge_fn,
                cut,
            );
        }
        // (THE CERTIFICATE'S OWN RSS PROBE, the receipts-phase memory
        // audit's blind window: `cut_index_lists` is the ONE structure
        // here whose cardinality is (|called| - 1) x |spectrum|, and at
        // the structurally degenerate shape (called.len() == pair_count,
        // the equality this code already relies on) that is
        // pair_count x pair_count. Whether it is the reported deep-locus
        // peak is a measurement, not an argument: the delta against the
        // `_scored` bookend prices it, and `_tie_emitted` prices the
        // Value-tree emission downstream.)
        rss.probe(&format!("locus_{locus}_certificate"))?;
        let tie_certificate = TieCertificate::CutIndices(&cut_index_lists);
        // (The phase census's certificate timer — one path for both
        // shapes now; the before-fix receipts' qual_bands_seconds/
        // qual_degenerate_seconds are this same certificate wall,
        // branch-exclusive.)
        let qual_certificate_seconds = certificate_started
            .map(|started| started.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        let qual_spectrum_seconds = spectrum_started
            .map(|started| started.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        let qual_cluster_started = locus_phase_census.then(Instant::now);
        let cluster = called_score
            .filter(|_| !spectrum.is_empty() || !called.is_empty())
            .map(|_| {
                if std::env::var("IMPG_TIED_FOLD_FAST").is_ok() {
                    tie_cluster_fast(
                        &folds,
                        &called,
                        &spectrum,
                        &observed_node_fn,
                        &observed_edge_fn,
                    )
                } else {
                    let pairs: Vec<(f64, f64)> = spectrum
                        .iter()
                        .map(|&(distance, score, _, _, _, _)| (distance, score))
                        .collect();
                    cluster_form_qual(1.0, &pairs, &tie_certificate)
                }
            });
        let qual_cluster_seconds = qual_cluster_started
            .map(|started| started.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        let cluster_knee = cluster.as_ref().and_then(|state| state.knee);
        let cluster_shape = cluster.as_ref().map(|state| state.shape);
        let cluster_size = cluster.as_ref().map(|state| state.cluster_size);
        let cluster_k = cluster.as_ref().map(|state| state.k);
        let cluster_alternative = cluster.as_ref().and_then(|state| state.alternative);
        let cluster_excluded = cluster
            .as_ref()
            .map(|state| state.excluded.iter().filter(|e| **e).count());
        let confidence = cluster.as_ref().and_then(|state| state.p);
        let qual = confidence.and_then(qual_from_p);
        let qual_unbounded = confidence.is_some_and(|value| value >= 1.0);
        let posterior_logsumexp = called_score.map(|score| {
            let total: f64 = class_lls.iter().map(|&ll| (ll - score).exp()).sum();
            score + total.ln()
        });
        let top_cluster_posterior = match (&cluster, posterior_logsumexp) {
            (Some(state), Some(logsumexp)) => {
                let mut mass = 0.0f64;
                for &(distance, score, _, _, _, _) in &spectrum {
                    if distance <= state.cut {
                        mass += score;
                    }
                }
                Some((1.0 + mass) / (logsumexp - called_score.unwrap_or(0.0)).exp())
            }
            _ => None,
        };
        let truth_in_called_set = truth_flat.map(|flat| {
            let (first, second) = unrank_pair(flat);
            called.contains(&(first, second))
        });
        // The best divergent rival outside the excluded band.
        // (THE BOUNDED SCAN'S BEST-OUTSIDE FALLBACK — the receipted
        // branch-and-bound form, SS10.6.21 (d), slack-corrected the
        // SS10.6.22 way: the skipped pairs' exact sums are unknown,
        // but each is bounded above by its own bound + the derived
        // rounding-error slack, so the argmax over the SKIPPED
        // non-excluded pairs is walked in bound-plus-slack
        // DESCENDING order (ties by flat index — deterministic)
        // until the next dominated bound falls to the running max
        // or below: every pair beyond has exact <= bound + slack <=
        // running max, so the walk's max is the exhaustive run's
        // own value. A second deterministic sweep then computes the
        // exact sums of every skipped non-excluded pair whose
        // bound + slack is at or above the found max — the tied
        // set's own members (a pair below the max's dominated bound
        // cannot tie it) — and patches them into class_lls, so the
        // emitted tied classes and the records' named pairs are the
        // exhaustive run's own. Where nothing is skipped this is the
        // committed fold verbatim.)
        let best_outside_started = locus_phase_census.then(Instant::now);
        let best_outside_log_likelihood = cluster.as_ref().map(|state| {
            let kept_max = spectrum
                .iter()
                .zip(&state.excluded)
                .filter(|(_, excluded)| !**excluded)
                .map(|(entry, _)| class_lls[entry.5])
                .fold(f64::NEG_INFINITY, f64::max);
            let skipped_nonexcluded = bounded_state
                .as_mut()
                .map(|bounded| {
                    let mut candidates: Vec<usize> = spectrum
                        .iter()
                        .zip(&state.excluded)
                        .filter(|(_, excluded)| !**excluded)
                        .map(|(entry, _)| entry.5)
                        .filter(|&flat| bounded.computed[flat].is_none())
                        .collect();
                    candidates.sort_by(|&a, &b| {
                        (bounded.pair_bounds[b] + bounded.pair_slack[b])
                            .total_cmp(&(bounded.pair_bounds[a] + bounded.pair_slack[a]))
                            .then(a.cmp(&b))
                    });
                    candidates
                })
                .unwrap_or_default();
            let mut running_max = kept_max;
            let mut next = 0usize;
            while next < skipped_nonexcluded.len() {
                let flat = skipped_nonexcluded[next];
                // the dominated bound: the pair's own bound + its
                // derived rounding-error slack (the same
                // slack-corrected dominance premise as stage one's
                // stop — the walk continues while it is strictly
                // above the running max)
                let bound = bounded_state
                    .as_ref()
                    .map(|bounded| bounded.pair_bounds[flat] + bounded.pair_slack[flat])
                    .unwrap_or(f64::NEG_INFINITY);
                if bound.partial_cmp(&running_max) != Some(std::cmp::Ordering::Greater) {
                    break;
                }
                let (first, second) = unrank_pair(flat);
                // the memoized exact sum (the fourth reduction at the
                // best-outside walk's own site): the canonical flat's
                // value, computed on the representative rows on a
                // miss — bit-identical to the direct form by
                // construction
                let sum = if let Some(bounded) = bounded_state.as_mut() {
                    match bounded.scan_memo.as_ref() {
                        None => matrix.class_pair_sum(first, second, &unit_counts),
                        Some(memo) => {
                            let canon = memo.canonical_flat(first, second);
                            match bounded.memo_values[canon] {
                                Some(value) => {
                                    memo.note_hits(1);
                                    value
                                }
                                None => {
                                    let (r_first, r_second) = unrank_pair(canon);
                                    let value =
                                        matrix.class_pair_sum(r_first, r_second, &unit_counts);
                                    memo.note_computed(1);
                                    bounded.memo_values[canon] = Some(value);
                                    value
                                }
                            }
                        }
                    }
                } else {
                    matrix.class_pair_sum(first, second, &unit_counts)
                };
                if let Some(bounded) = bounded_state.as_mut() {
                    bounded.computed[flat] = Some(sum);
                }
                class_lls[flat] = sum;
                if sum > running_max {
                    running_max = sum;
                }
                next += 1;
            }
            // the tied-set sweep: every skipped non-excluded pair at
            // or above the found max, exact — deterministic.
            if let Some(bounded) = bounded_state.as_mut() {
                let tied: Vec<usize> = skipped_nonexcluded
                    .iter()
                    .copied()
                    .filter(|&flat| {
                        bounded.computed[flat].is_none()
                            && bounded.pair_bounds[flat] + bounded.pair_slack[flat] >= running_max
                    })
                    .collect();
                for flat in tied {
                    let (first, second) = unrank_pair(flat);
                    // the memoized exact sum — the same canonical-slot
                    // discipline as the walk above
                    let sum = match bounded.scan_memo.as_ref() {
                        None => matrix.class_pair_sum(first, second, &unit_counts),
                        Some(memo) => {
                            let canon = memo.canonical_flat(first, second);
                            match bounded.memo_values[canon] {
                                Some(value) => {
                                    memo.note_hits(1);
                                    value
                                }
                                None => {
                                    let (r_first, r_second) = unrank_pair(canon);
                                    let value =
                                        matrix.class_pair_sum(r_first, r_second, &unit_counts);
                                    memo.note_computed(1);
                                    bounded.memo_values[canon] = Some(value);
                                    value
                                }
                            }
                        }
                    };
                    bounded.computed[flat] = Some(sum);
                    class_lls[flat] = sum;
                }
            }
            running_max
        })
        .filter(|value| value.is_finite());
        let best_outside_seconds = best_outside_started
            .map(|started| started.elapsed().as_secs_f64())
            .unwrap_or(0.0);
        // (The phase census's shape counters — computed here at write
        // time, O(folds + lists): the entry counts that price the
        // certificate phases and the per-fold multiset lengths that
        // price every multiset walk. Zero and untouched with the
        // census env unset.)
        let phase_census_tied_list_entries: u64 = if locus_phase_census {
            cut_index_lists.iter().map(|(_, list)| list.len() as u64).sum()
        } else {
            0
        };
        let phase_census_tied_list_count: usize = if locus_phase_census {
            cut_index_lists.len()
        } else {
            0
        };
        let phase_census_fold_multiset_entries: u64 = if locus_phase_census {
            folds
                .iter()
                .map(|fold| (fold.nodes.len() + fold.edges.len()) as u64)
                .sum()
        } else {
            0
        };
        let phase_census_fold_multiset_max: usize = if locus_phase_census {
            folds
                .iter()
                .map(|fold| fold.nodes.len() + fold.edges.len())
                .max()
                .unwrap_or(0)
        } else {
            0
        };
        let mut qual_alternative_classes: Vec<serde_json::Value> = Vec::new();
        let alternative_log_gap =
            best_outside_log_likelihood.map(|outside| best_ll - outside);
        if let (Some(state), Some(best_outside)) = (&cluster, best_outside_log_likelihood) {
            for (index, entry) in spectrum.iter().enumerate() {
                let &(distance, score, signature_distance, first, second, flat) = entry;
                if !state.excluded[index] && class_lls[flat] == best_outside {
                    qual_alternative_classes.push(json!({
                        "spectrum_index": index,
                        "distance": distance,
                        "signature_distance": signature_distance,
                        "relative_likelihood": score,
                        "log_likelihood": class_lls[flat],
                        "fold_indices": [first, second],
                    }));
                }
            }
        }
        let nearest_rival_distance = spectrum
            .iter()
            .find(|&&(distance, ..)| distance > 0.0)
            .map(|&(distance, ..)| distance);
        let spectrum_zero_distance = called_score
            .map(|_| spectrum.iter().filter(|&&(distance, ..)| distance == 0.0).count());
        let qual_seconds = qual_started.elapsed().as_secs_f64();

        // ---------------------------------------- the exactness sample
        // Strata: the evenly-spaced valid (unit, fold) pairs; the
        // evenly-spaced pairs WITH mismatch votes (the variant-pocket
        // class); the truth folds; the identical-through-graph fold
        // (the named BTE#3/#4 block28_contig1 pair, when it is a member
        // of this locus's partition). The in-process factorization
        // assertion already covers the WHOLE pilot domain; the sample
        // carries vote-level detail for the independent checker.
        let exactness_started = Instant::now();
        // (Candidate C: the pinned-pair enumeration walks the rows'
        // own pinned entries — the same (fold, unit) pairs in the same
        // order the dense filter produced, with the mismatches read
        // from the entry itself.)
        //
        // THE TWO VECTORS ARE NOT MATERIALIZED (the receipts-phase
        // memory audit's Finding 3): the pinned-pair list is used for
        // exactly two things - its LENGTH (to derive each stratum's
        // stride) and every stride-th entry - so it is walked twice
        // and never held. At a deep locus it was ~49.6M pairs x 16B =
        // 794MB (up to ~1.6GB mid-growth, plus the same again for the
        // mismatch vector), retained past the whole QUAL phase into
        // the receipt emission. The walk order is identical, so the
        // sampled set is identical, by construction.
        let mut valid_len = 0usize;
        let mut mismatch_len = 0usize;
        for fold in 0..n_folds {
            for (_, mismatches) in matrix.pinned_units(fold) {
                valid_len += 1;
                if mismatches > 0 {
                    mismatch_len += 1;
                }
            }
        }
        let stride_of = |len: usize, want: u64| -> usize {
            if len == 0 || want == 0 {
                usize::MAX
            } else {
                ((len + want as usize - 1) / want as usize).max(1)
            }
        };
        let mut sampled: BTreeSet<(usize, usize)> = BTreeSet::new();
        let mut stratum_counts: BTreeMap<String, u64> = BTreeMap::new();
        {
            let stride = stride_of(valid_len, options.exactness_sample);
            let mut taken = 0u64;
            let mut index = 0usize;
            for fold in 0..n_folds {
                for (unit, _) in matrix.pinned_units(fold) {
                    if index % stride == 0 {
                        sampled.insert((fold, unit));
                        taken += 1;
                    }
                    index += 1;
                }
            }
            stratum_counts.insert("evenly_spaced_valid".to_string(), taken);
        }
        {
            let stride = stride_of(mismatch_len, options.exactness_sample);
            let mut taken = 0u64;
            let mut index = 0usize;
            for fold in 0..n_folds {
                for (unit, mismatches) in matrix.pinned_units(fold) {
                    if mismatches > 0 {
                        if index % stride == 0 && sampled.insert((fold, unit)) {
                            taken += 1;
                        }
                        index += 1;
                    }
                }
            }
            stratum_counts.insert("variant_pocket_mismatch".to_string(), taken);
        }
        let mut named_folds: Vec<(String, usize)> = Vec::new();
        if let Some(first) = truth_first {
            named_folds.push(("truth_first".to_string(), first));
        }
        if let Some(second) = truth_second {
            named_folds.push(("truth_second".to_string(), second));
        }
        if let Some(index) = folds.iter().position(|fold| {
            fold.members
                .iter()
                .any(|m| m.path_name == "BTE#3#block28_contig1")
        }) {
            named_folds.push(("identical_pair".to_string(), index));
        }
        for (name, fold_index) in &named_folds {
            let pairs: Vec<(usize, usize)> = matrix
                .pinned_units(*fold_index)
                .map(|(unit, _)| (*fold_index, unit))
                .collect();
            let stride = stride_of(pairs.len(), 25);
            let mut taken = 0u64;
            for (index, pair) in pairs.iter().enumerate() {
                if index % stride == 0 && sampled.insert(*pair) {
                    taken += 1;
                }
            }
            stratum_counts.insert(format!("named_{name}"), taken);
        }
        for &(fold_index, unit_index) in &sampled {
            let unit = &units[unit_index];
            let record = unit.record;
            let canonical =
                decode_record_tokens(&key_tokens[census_key[record].unwrap() as usize])?;
            let fold = &folds[fold_index];
            let (img_lo, img_hi) = territory[fold_index];
            let positions_of = |node: u32| fold.positions_of(node);
            let mut placements: Vec<serde_json::Value> = Vec::new();
            let mut local_scores: Vec<f64> = Vec::new();
            for &orientation in &touching[&record].orientations {
                let (pinned, _ambiguous) =
                    anchor_correspondence(&canonical, orientation, &positions_of);
                if pinned.iter().all(|p| p.is_none()) {
                    continue;
                }
                let windows_inside =
                    pinned.iter().flatten().all(|&r| r >= img_lo && r + k <= img_hi);
                if !windows_inside {
                    continue;
                }
                let own_flat = own_positions_of(&canonical, false, 0);
                let windows: Vec<(u64, u64)> = pinned
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &p)| p.map(|r| (own_flat[j], r as u64)))
                    .map(|(p, _)| (p, p + k))
                    .collect();
                let backbone_bp: u64 =
                    merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
                let pinned_i: Vec<Option<i64>> =
                    pinned.iter().map(|p| p.map(|v| v as i64)).collect();
                let pin = Pin {
                    orientation,
                    pinned: pinned_i,
                    backbone_bp,
                };
                let Some((m, c)) = placement_score(
                    unit, &canonical, &pin, &fold.seq, fold.len, &reads, k,
                    img_lo as i64, img_hi as i64,
                ) else {
                    continue;
                };
                // The vote-level detail (the skipped bases, their
                // projected coordinates and their votes).
                let forward = physical_forward(orientation, unit.mirror);
                let own_positions = own_positions_of(&canonical, unit.mirror, unit.w_lo);
                let own_hull = own_positions_of(&canonical, unit.mirror, 0);
                let skipped = skipped_positions(&own_positions, k, read_length as u64);
                let mut votes: Vec<serde_json::Value> = Vec::new();
                for &(lo, hi) in &skipped {
                    for i in lo..hi {
                        let d = i as i64 - unit.w_lo as i64;
                        let (index, own_pos) =
                            serving_anchor(&pin.pinned, &own_hull, d).expect("serving anchor");
                        let r = pin.pinned[index].unwrap();
                        let coord = if forward {
                            r + (d - own_pos)
                        } else {
                            r + (own_pos + k as i64 - 1 - d)
                        };
                        let base = fold.seq[coord as usize];
                        let read_base = if forward {
                            reads[unit.read][i as usize]
                        } else {
                            complement(reads[unit.read][i as usize])
                        };
                        votes.push(json!([i, coord, if read_base == base { 1 } else { 0 }]));
                    }
                }
                placements.push(json!({
                    "orientation": orientation,
                    "pinned": pin.pinned,
                    "backbone_bp": backbone_bp,
                    "m": m,
                    "c": c,
                    "score": scoring.placement(m, c),
                    "votes": votes,
                }));
                local_scores.push(scoring.placement(m, c));
            }
            let read = &reads[unit.read];
            serde_json::to_writer(
                &mut exactness_buf,
                &json!({
                    "locus": locus,
                    "unit": unit_index,
                    "record": record,
                    "variant": unit.variant,
                    "mirror": if unit.mirror { 1 } else { 0 },
                    "w_lo": unit.w_lo,
                    "count": unit.count,
                    "read_fnv": fnv1a64(read),
                    "read": String::from_utf8(read.clone()).unwrap(),
                    "fold": fold_index,
                    "ll": matrix.ll(fold_index, unit_index),
                    "ll_local": unit_log_likelihood(&local_scores, img_hi - img_lo, read_length),
                    "elsewhere": elsewhere,
                    "placements": placements,
                }),
            )?;
            exactness_buf.push(b'\n');
        }
        let exactness_seconds = exactness_started.elapsed().as_secs_f64();
        eprintln!(
            "[score] locus {locus}: exactness sample {} pairs of {} valid \
             ({} with mismatch votes) [{exactness_seconds:.1}s]",
            sampled.len(),
            valid_len,
            mismatch_len,
        );

        // -------------------------------------- the pocket-read anatomy
        // (slice C gate 1: THE ANATOMY THAT DECIDES THE FORM of the
        // generated-here/generated-elsewhere marginalization. For
        // every unit, classified by whether the winner and truth
        // pairs place it, the placement anatomy on the winner and
        // truth folds: the anchored/skipped/abstained decomposition,
        // the collinearity verdict, the slice-B piecewise (m, c)
        // BESIDE the one-continuous-path whole-read score, and the
        // flank votes with their projected coordinates — plus the
        // read's full verified-occurrence census and, on the
        // sample's own donor paths (S288C/SK1, the balanced 15x/15x
        // construction), the FULL-READ match at each donor
        // occurrence: the generated-elsewhere evidence. Pure
        // emission; the scoring is unchanged.)
        let anatomy_seconds = if anatomy_enabled {
            let anatomy_started = Instant::now();
            let (winner_first, winner_second) = winner_pair;
            let mut anatomy_folds: Vec<usize> = vec![winner_first, winner_second];
            if let Some((a, b)) = truth_pair {
                for fold in [a, b] {
                    if !anatomy_folds.contains(&fold) {
                        anatomy_folds.push(fold);
                    }
                }
            }
            let truth_places_unit = |u: usize| {
                truth_pair
                    .map(|(a, b)| matrix.pinned(a, u) || matrix.pinned(b, u))
                    .unwrap_or(false)
            };
            let winner_places_unit =
                |u: usize| matrix.pinned(winner_first, u) || matrix.pinned(winner_second, u);
            let is_donor_path = |name: &str| {
                name.starts_with("S288C#0#") || name.starts_with("SK1#0#")
            };
            let mut anatomy_units = 0u64;
            for (u, unit) in units.iter().enumerate() {
                let record = unit.record;
                let canonical =
                    decode_record_tokens(&key_tokens[census_key[record].unwrap() as usize])?;
                let line = &census_records[record];
                let shifts = &record_shifts[record];
                let span = walk_span(&canonical, k);
                let class = match (winner_places_unit(u), truth_places_unit(u)) {
                    (true, false) => "winner_only",
                    (true, true) => "both",
                    (false, true) => "truth_only",
                    (false, false) => "neither",
                };
                let pinning_folds: Vec<usize> = (0..n_folds)
                    .filter(|&f| matrix.pinned(f, u))
                    .collect();
                // The donor-path full-read verification: at each
                // verified occurrence of this record on a donor
                // path, the whole read (all 150 bases, strand-aware,
                // origin-shift-corrected) against the donor path's
                // own sequence.
                let mut donor_checks: Vec<serde_json::Value> = Vec::new();
                for (position, occ) in line.occurrences.iter().enumerate() {
                    let name = &panel.name_map.path_to_name[occ.path];
                    if !is_donor_path(name) {
                        continue;
                    }
                    let origin = (occ.start as i64 + shifts[position]).max(0) as u64;
                    let path_len = panel.name_map.path_to_length[occ.path];
                    let forward = physical_forward(occ.orientation, unit.mirror);
                    let lo = origin.saturating_sub(300);
                    let hi = (origin + span + 300).min(path_len);
                    let seq = fetch_seq(name, lo, hi)?;
                    let read = &reads[unit.read];
                    let mut full_m = 0u32;
                    let mut full_c = 0u32;
                    let mut full_oob = 0u32;
                    for i in 0..read.len() as i64 {
                        let d = i - unit.w_lo as i64;
                        let coord = if forward {
                            origin as i64 + d
                        } else {
                            origin as i64 + span as i64 - 1 - d
                        };
                        if coord < 0
                            || coord >= path_len as i64
                            || (coord as u64) < lo
                            || coord as u64 >= hi
                        {
                            full_oob += 1;
                            continue;
                        }
                        let base = seq[(coord - lo as i64) as usize];
                        let read_base = if forward {
                            read[i as usize]
                        } else {
                            complement(read[i as usize])
                        };
                        if read_base == base {
                            full_m += 1;
                        } else {
                            full_c += 1;
                        }
                    }
                    donor_checks.push(json!({
                        "path": name, "origin": origin,
                        "orientation": occ.orientation,
                        "m": full_m, "c": full_c, "oob": full_oob,
                    }));
                }
                // The placement anatomy on the winner and truth
                // folds (the pocket rows and the rows the truth
                // spells).
                let emit_votes = class == "winner_only" || u % 16 == 0;
                let mut placements: Vec<serde_json::Value> = Vec::new();
                for &fold_index in &anatomy_folds {
                    let fold = &folds[fold_index];
                    let (img_lo, img_hi) = territory[fold_index];
                    let positions_of = |node: u32| fold.positions_of(node);
                    for &orientation in &touching[&record].orientations {
                        let (pinned_vec, _ambiguous) =
                            anchor_correspondence(&canonical, orientation, &positions_of);
                        // The per-anchor candidate census (the diagnosis of
                        // all-None correspondences: which anchors the
                        // fold's walk carries, with which signs).
                        let anchor_candidates: Vec<serde_json::Value> = canonical
                            .iter()
                            .map(|&(node, _)| {
                                let required = required_sign(node, orientation);
                                let positions = positions_of(required.unsigned_abs());
                                let with_sign = positions
                                    .iter()
                                    .filter(|&&(_, sign)| sign == required.signum())
                                    .count();
                                json!([
                                    required,
                                    with_sign,
                                    positions.len(),
                                ])
                            })
                            .collect();
                        if pinned_vec.iter().all(|p| p.is_none()) {
                            placements.push(json!({
                                "fold": fold_index,
                                "orientation": orientation,
                                "valid": false,
                                "windows_inside": false,
                                "pinned": pinned_vec,
                                "anchor_candidates": anchor_candidates,
                                "backbone_bp": 0,
                                "abstained_bp": 0,
                                "skipped_bp": 0,
                                "collinear": false,
                                "offset": null,
                                "m": 0, "c": 0,
                                "single": null,
                                "votes": [],
                            }));
                            continue;
                        }
                        let windows_inside =
                            pinned_vec.iter().flatten().all(|&r| r >= img_lo && r + k <= img_hi);
                        let own_flat = own_positions_of(&canonical, false, 0);
                        let backbone_windows: Vec<(u64, u64)> = pinned_vec
                            .iter()
                            .enumerate()
                            .filter_map(|(j, &r)| r.map(|r2| (own_flat[j], r2)))
                            .map(|(p, _)| (p, p + k))
                            .collect();
                        let backbone_bp: u64 =
                            merge_windows(backbone_windows).iter().map(|&(lo, hi)| hi - lo).sum();
                        let hull_windows: Vec<(u64, u64)> =
                            own_flat.iter().map(|&p| (p, p + k)).collect();
                        let hull_covered: u64 = merge_windows(hull_windows)
                            .iter()
                            .map(|&(lo, hi)| hi - lo)
                            .sum();
                        let pinned_i: Vec<Option<i64>> =
                            pinned_vec.iter().map(|p| p.map(|v| v as i64)).collect();
                        let pin = Pin {
                            orientation,
                            pinned: pinned_i.clone(),
                            backbone_bp,
                        };
                        let forward = physical_forward(orientation, unit.mirror);
                        let own_hull = own_positions_of(&canonical, unit.mirror, 0);
                        let (collinear, offset_const) =
                            placement_offset(&pin.pinned, &own_hull, forward);
                        let own_positions =
                            own_positions_of(&canonical, unit.mirror, unit.w_lo);
                        let skipped = skipped_positions(&own_positions, k, read_length as u64);
                        let skipped_bp: u64 =
                            skipped.iter().map(|&(lo, hi)| hi - lo).sum();
                        let piecewise = placement_score(
                            unit, &canonical, &pin, &fold.seq, fold.len, &reads, k,
                            img_lo as i64, img_hi as i64,
                        );
                        let (m, c, valid) = match piecewise {
                            Some((m, c)) => (m, c, true),
                            None => (0, 0, false),
                        };
                        let single = if collinear {
                            Some(single_offset_score(
                                unit,
                                offset_const.unwrap(),
                                forward,
                                &fold.seq,
                                fold.len,
                                &reads,
                                k,
                            ))
                        } else {
                            None
                        };
                        let mut votes: Vec<serde_json::Value> = Vec::new();
                        if emit_votes {
                            if let Some((m, c)) = piecewise {
                                let _ = (m, c);
                                for &(lo, hi) in &skipped {
                                    for i in lo..hi {
                                        let d = i as i64 - unit.w_lo as i64;
                                        let (index, own_pos) = match serving_anchor(
                                            &pin.pinned, &own_hull, d,
                                        ) {
                                            Some(found) => found,
                                            None => continue,
                                        };
                                        let r = match pin.pinned[index] {
                                            Some(r) => r,
                                            None => continue,
                                        };
                                        let coord = if forward {
                                            r + (d - own_pos)
                                        } else {
                                            r + (own_pos + k as i64 - 1 - d)
                                        };
                                        if coord < 0 || coord >= fold.len as i64 {
                                            votes.push(json!([i, -1, 0]));
                                            continue;
                                        }
                                        let base = fold.seq[coord as usize];
                                        let read_base = if forward {
                                            reads[unit.read][i as usize]
                                        } else {
                                            complement(reads[unit.read][i as usize])
                                        };
                                        votes.push(json!([
                                            i,
                                            coord,
                                            if read_base == base { 1 } else { 0 },
                                        ]));
                                    }
                                }
                            }
                        }
                        placements.push(json!({
                            "fold": fold_index,
                            "orientation": orientation,
                            "valid": valid,
                            "windows_inside": windows_inside,
                            "pinned": pin.pinned,
                            "anchor_candidates": anchor_candidates,
                            "backbone_bp": backbone_bp,
                            "abstained_bp": hull_covered - backbone_bp,
                            "skipped_bp": skipped_bp,
                            "collinear": collinear,
                            "offset": offset_const,
                            "m": m, "c": c,
                            "single": single,
                            "votes": votes,
                        }));
                    }
                }
                let read = &reads[unit.read];
                serde_json::to_writer(
                    &mut anatomy_buf,
                    &json!({
                        "type": "unit",
                        "locus": locus,
                        "unit": u,
                        "record": record,
                        "variant": unit.variant,
                        "mirror": if unit.mirror { 1 } else { 0 },
                        "w_lo": unit.w_lo,
                        "count": unit.count,
                        "read_fnv": fnv1a64(read),
                        "read": String::from_utf8(read.clone()).unwrap(),
                        "class": class,
                        "pinning_folds": pinning_folds,
                        "donor_checks": donor_checks,
                        "placements": placements,
                    }),
                )?;
                anatomy_buf.push(b'\n');
                anatomy_units += 1;
            }
            eprintln!(
                "[score] locus {locus}: anatomy {anatomy_units} units \
                 [{:.1}s]",
                anatomy_started.elapsed().as_secs_f64(),
            );
            anatomy_started.elapsed().as_secs_f64()
        } else {
            0.0
        };

        // ------------------------------------------------- the named classes
        // (Phase instrumentation: the receipts/IO phase — the named
        // classes, the ingredients (with the full ll_matrix), and the
        // main receipt write.)
        let receipts_started = Instant::now();
        let mut named: Vec<(usize, usize)> = Vec::new();
        named.push(winner_pair);
        if let Some(flat) = truth_flat {
            let pair = unrank_pair(flat);
            if !named.contains(&pair) {
                named.push(pair);
            }
        }
        if let (Some(state), Some(_)) = (&cluster, best_outside_log_likelihood) {
            for (index, entry) in spectrum.iter().enumerate() {
                if !state.excluded[index] && class_lls[entry.5] == best_outside_log_likelihood.unwrap() {
                    let pair = (entry.3, entry.4);
                    if !named.contains(&pair) {
                        named.push(pair);
                    }
                    break;
                }
            }
        }
        let mut classes_evidence: Vec<serde_json::Value> = Vec::new();
        for &(first, second) in &named {
            let per_unit: Vec<serde_json::Value> = units
                .iter()
                .enumerate()
                .map(|(u, unit)| {
                    json!({
                        "record": unit.record,
                        "variant": unit.variant,
                        "count": unit.count,
                        "fold_ll": [matrix.ll(first, u), matrix.ll(second, u)],
                        "pair_term": mix_logsumexp(matrix.ll(first, u), matrix.ll(second, u)),
                    })
                })
                .collect();
            classes_evidence.push(json!({
                "fold_indices": [first, second],
                "log_likelihood": class_lls[second * (second + 1) / 2 + first],
                "per_unit": per_unit,
            }));
        }
        serde_json::to_writer(&mut records_buf, &json!({
            "locus": locus,
            "classes": classes_evidence,
        }))?;
        records_buf.push(b'\n');

        // ----------------------------------------------------- the ingredients
        // (Phase 1, lever 3 — the streaming serializer: the ll_matrix is
        // serialized DIRECTLY from the typed matrix, not materialized
        // as a serde_json value tree first. Phase 0 measured the L16
        // tree at +1.15GB live RSS over the scoring baseline (34.5M
        // entries); the streamed bytes are identical — the json! macro
        // emits object keys in serde_json's default (non-preserved)
        // map order, i.e. alphabetical, and every field here
        // serializes through the same writer in that same order. The
        // identity gate proves the sidecar byte-identical.)
        // (THE INGREDIENT ARRAYS STREAM PER ELEMENT - the receipts-phase
        // memory audit's Finding 2, the lever-3 form extended from
        // ll_matrix to its two siblings. Building them as
        // `Vec<serde_json::Value>` first cost several GB at a deep locus
        // (per fold: members, the whole sequence, a walk/nodes/edges
        // array pair per step; per unit: six keys including a heap copy
        // of the read sequence) and held it through the whole receipt
        // emission. Each element is now built, serialized and dropped in
        // turn; the bytes are identical, because a Vec's serialization
        // IS its elements' serializations joined by ',' inside "[" "]",
        // and json!'s maps are alphabetical in both forms.)
        // (Candidate C: the ingredients line — including the full
        // ll_matrix — streams to a per-locus spill file, the dense
        // JSON text never held in RAM; the ordered write below copies
        // the spill into the ingredients file, so the emitted bytes
        // are the committed form's exactly.)
        let ingredients_spill_path =
            ingredients_spill_path(&options.out, locus);
        let mut ingredients_spill =
            BufWriter::new(File::create(&ingredients_spill_path)?);
        ingredients_spill.write_all(b"{\"folds\":[")?;
        for (position, fold) in folds.iter().enumerate() {
            if position > 0 {
                ingredients_spill.write_all(b",")?;
            }
            serde_json::to_writer(
                &mut ingredients_spill,
                &json!({
                    "members": fold.members.iter().map(|m| json!({
                        "path_name": m.path_name, "start": m.start, "end": m.end,
                    })).collect::<Vec<_>>(),
                    "length": fold.len,
                    "sequence": String::from_utf8(fold.seq.clone()).unwrap(),
                    "walk": fold.walk.iter().map(|&(bp, node)| json!([bp, node]))
                        .collect::<Vec<_>>(),
                    "nodes": fold.nodes,
                    "edges": fold.edges,
                }),
            )?;
        }
        ingredients_spill.write_all(b"]")?;
        // The per-unit per-fold log-likelihood matrix (folds x units,
        // row-major; entries at E = no valid anchored placement).
        ingredients_spill.write_all(b",\"ll_matrix\":")?;
        matrix.write_ll_matrix_json(&mut ingredients_spill)?;
        ingredients_spill.write_all(b",\"locus\":")?;
        serde_json::to_writer(&mut ingredients_spill, &locus)?;
        ingredients_spill.write_all(b",\"partition\":")?;
        serde_json::to_writer(&mut ingredients_spill, &partition)?;
        ingredients_spill.write_all(b",\"units\":[")?;
        for (position, unit) in units.iter().enumerate() {
            if position > 0 {
                ingredients_spill.write_all(b",")?;
            }
            serde_json::to_writer(
                &mut ingredients_spill,
                &json!({
                    "record": unit.record,
                    "variant": unit.variant,
                    "mirror": if unit.mirror { 1 } else { 0 },
                    "w_lo": unit.w_lo,
                    "count": unit.count,
                    "read": String::from_utf8(reads[unit.read].clone()).unwrap(),
                }),
            )?;
        }
        ingredients_spill.write_all(b"]")?;
        ingredients_spill.write_all(b"}")?;
        ingredients_spill.write_all(b"\n")?;
        ingredients_spill.flush()?;
        drop(ingredients_spill);
        rss.probe(&format!("locus_{locus}_receipts"))?;

        // ---------------------------------------------------- the main receipt
        let called_classes_json: Vec<serde_json::Value> = called
            .iter()
            .map(|&(first, second)| {
                json!({
                    "fold_indices": [first, second],
                    "physical_pairs": class_physical_pairs(
                        first == second,
                        folds[first].members.len(),
                        folds[second].members.len(),
                    ),
                })
            })
            .collect();
        let fold_identities: Vec<serde_json::Value> = folds
            .iter()
            .enumerate()
            .map(|(fold_index, fold)| {
                let (img_lo, img_hi) = territory[fold_index];
                let mut value = json!({
                    "members": fold.members.iter().map(|m| json!({
                        "path_name": m.path_name, "start": m.start, "end": m.end,
                    })).collect::<Vec<_>>(),
                    "length": fold.len,
                });
                if territory_normalized {
                    value["image"] = if img_hi > img_lo {
                        json!([img_lo, img_hi])
                    } else {
                        serde_json::Value::Null
                    };
                }
                value
            })
            .collect();
        let class_fold_pairs: Vec<[usize; 2]> = (0..pair_count)
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                [first, second]
            })
            .collect();
        let identical_pair_fold = folds.iter().position(|fold| {
            fold.members.iter().any(|m| m.path_name == "BTE#3#block28_contig1")
                && fold.members
                    .iter()
                    .any(|m| m.path_name == "BTE#4#block28_contig1")
        });
        // The log-gap floor attribution: how much of the winner-truth
        // gap is the all-mismatch floor (reads the pair can place on
        // neither homolog) versus per-base evidence on the reads both
        // place. The floor's magnitude is the model's derived minimum,
        // not the true unanchored likelihood (named).
        let pair_places =
            |i: usize, j: usize, u: usize| matrix.pinned(i, u) || matrix.pinned(j, u);
        let mut truth_unplaced_units = 0u64;
        let mut truth_unplaced_mass = 0.0f64;
        let mut winner_unplaced_units = 0u64;
        let mut winner_unplaced_mass = 0.0f64;
        let mut both_place_gap = 0.0f64;
        if let (Some((ti, tj)), Some((wi, wj))) = (truth_pair, Some(winner_pair)) {
            for (u, unit) in units.iter().enumerate() {
                let truth_places = pair_places(ti, tj, u);
                let winner_places = pair_places(wi, wj, u);
                if !truth_places {
                    truth_unplaced_units += 1;
                    truth_unplaced_mass += unit.count as f64;
                }
                if !winner_places {
                    winner_unplaced_units += 1;
                    winner_unplaced_mass += unit.count as f64;
                }
                if truth_places && winner_places {
                    both_place_gap += unit.count as f64
                        * (mix_logsumexp(matrix.ll(wi, u), matrix.ll(wj, u))
                            - mix_logsumexp(matrix.ll(ti, u), matrix.ll(tj, u)));
                }
            }
        }
        let locality_totals = {
            let mut in_axis = 0u64;
            let mut overhang = 0u64;
            let mut extension = 0u64;
            for &record in &locus_records {
                let data = &touching[&record];
                in_axis += data.in_axis;
                overhang += data.overhang;
                extension += data.extension;
            }
            (in_axis, overhang, extension)
        };
        let locus_seconds = locus_started.elapsed().as_secs_f64();
        let rss_kb = rss.probe(&format!("locus_{locus}"))?;
        serde_json::to_writer(&mut report_buf, &json!({
            "locus": locus,
            "partition": partition,
            "model": model_name,
            "skeleton": {
                "scheme": if stored_frame { "stored_walk" } else { "canonical_scheme" },
                "rows": skeleton_rows,
                "rows_affected": skeleton_affected,
                "added": skeleton_added,
                "replaced": skeleton_replaced,
                "dropped": skeleton_dropped,
                "kept_forward": skeleton_forward,
                "kept_reverse": skeleton_reverse,
                "note": "the canonical-scheme pin skeleton (slice D's frame \
                         repair: per position the frame that spells the \
                         k-mer's canonical form forward decides; every step \
                         sequence-verified against the fetched panel sequence); \
                         added = rc-frame-only anchors the stored walk lacks",
            },
            "scoring": {
                "phred": phred,
                "epsilon": epsilon,
                "match_log_prob": scoring.a,
                "mismatch_log_prob": scoring.b,
                "read_length": read_length,
                "uniform_qualities": true,
                "elsewhere_log_prob": elsewhere,
                "genome_diploid_bp": genome_diploid_bp,
                "panel_mean_haploid_bp": panel_mean_haploid_bp,
                "panel_strains": strain_totals.len(),
                "panel_strain_list": strain_list,
            },
            "record_count": locus_records.len(),
            "unit_count": n_units,
            "folds": n_folds,
            "physical_members": map.members.len(),
            "class_count": pair_count,
            "class_fold_pairs": class_fold_pairs,
            // null entries (the exhaustive form) = no valid anchored
            // placement for any unit of the class (the anchored
            // evidence model's floor); null entries (the labeled
            // bounded form, stated by the -classscan-bounded model
            // string and the gated classscan_bounded certificate
            // member) = SKIPPED by the derived bound — the pair's
            // exact sum sits below best + the derived exp-underflow
            // point (its every QUAL contribution is exactly +0.0),
            // its per-pair bound and the bound's derived
            // rounding-error slack are emitted in the certificate,
            // and the answers are the exhaustive run's own by the
            // bound-plus-slack's arithmetic and the final-value
            // tie-set rebuild (SS10.6.21, SS10.6.22).
            "class_log_likelihoods": class_lls,
            "best_log_likelihood": best_ll,
            "best_fold_indices": [winner_pair.0, winner_pair.1],
            "called_class_count": called.len(),
            "called_classes": called_classes_json,
            "fold_identities": fold_identities,
            "truth_folds": truth_pair.map(|(a, b)| [a, b]),
            "truth_pair_expressible": truth_pair.is_some(),
            "truth_log_likelihood": truth_log_likelihood,
            "truth_rank": truth_rank,
            "truth_tied_classes": truth_tied_classes,
            "truth_in_called_set": truth_in_called_set,
            "log_gap": log_gap,
            "identical_pair_fold": identical_pair_fold,
            "qual": qual,
            "qual_unbounded": qual_unbounded,
            "qual_spectrum_shape": cluster_shape,
            "qual_knee_distance": cluster_knee,
            "qual_cluster_size": cluster_size,
            "qual_cluster_k": cluster_k,
            "qual_spectrum_classes": called_score.map(|_| spectrum.len()),
            "qual_spectrum_zero_distance_classes": spectrum_zero_distance,
            "qual_nearest_rival_distance": nearest_rival_distance,
            "qual_excluded_class_count": cluster_excluded,
            "qual_distance_spectrum": spectrum.iter().map(|s| s.0).collect::<Vec<_>>(),
            "qual_distance_spectrum_scores": spectrum.iter().map(|s| s.1).collect::<Vec<_>>(),
            "qual_signature_distance_spectrum": spectrum.iter().map(|s| s.2).collect::<Vec<_>>(),
            "qual_alternative_similarity": cluster_alternative,
            "qual_alternative_log_gap": alternative_log_gap,
            "qual_alternative_classes": qual_alternative_classes,
            "qual_posterior_logsumexp": posterior_logsumexp,
            "qual_posterior_top_cluster": top_cluster_posterior,
            "observed_node_mass": observed_node_mass,
            "observed_edge_mass": observed_edge_mass,
            "observed_node_keys": observed_nodes.len(),
            "observed_edge_keys": observed_edges.len(),
            "log_gap_decomposition": {
                "truth_unplaced_units": truth_unplaced_units,
                "truth_unplaced_mass": truth_unplaced_mass,
                "winner_unplaced_units": winner_unplaced_units,
                "winner_unplaced_mass": winner_unplaced_mass,
                "both_placed_evidence_gap": both_place_gap,
                "note": "unplaced reads score at the derived ELSEWHERE \
                         branch E (the generated-elsewhere marginalization; \
                         the genome outside this locale at the measured \
                         per-base rate under the max-entropy uniform origin), \
                         not the all-mismatch 150*B; ranks are bounded by E, \
                         gap magnitudes are E-dominated where unplaced \
                         counts differ",
            },
            "locality": {
                "touching_in_axis": locality_totals.0,
                "touching_overhang": locality_totals.1,
                "touching_extension": locality_totals.2,
                "examples": locality_examples,
            },
            "validation": {
                "exactness_assert": exactness_assert,
                "factorized_checked": factorized_checked,
                "factorized_equal": factorized_equal,
                "edge_discarded_placements": edge_discarded,
                "floor_only_unit_fold": no_pin,
                "abstained_hull_bp": abstained_hull_bp,
                "exactness_sample_pairs": sampled.len(),
                "exactness_valid_pairs": valid_len,
                "exactness_mismatch_pairs": mismatch_len,
                "exactness_strata": stratum_counts,
            },
            "walls": {
                "folds_seconds": fold_started.elapsed().as_secs_f64() - scoring_seconds - class_seconds - qual_seconds - exactness_seconds - anatomy_seconds,
                "scoring_seconds": scoring_seconds,
                "classes_seconds": class_seconds,
                "qual_seconds": qual_seconds,
                "exactness_seconds": exactness_seconds,
                "anatomy_seconds": anatomy_seconds,
                "locus_seconds": locus_seconds,
            },
            "rss_kb": rss_kb,
        }))?;
        // (Stage 2: the territory field is emitted ONLY under the
        // normalized convention — the default receipt stays
        // byte-identical to the committed convention, the identity
        // gate. serde_json's default map is ordered, so the late
        // insertion lands alphabetically like every other key.)
        if territory_normalized {
            let mut report_value: serde_json::Value =
                serde_json::from_slice(&report_buf)?;
            report_value["territory"] = json!({
                "convention": "normalized",
                "axis_fold": territory_axis_fold,
                "images": territory.iter().map(|&(lo, hi)| {
                    if hi > lo { json!([lo, hi]) } else { serde_json::Value::Null }
                }).collect::<Vec<_>>(),
                "note": "the locus's normalized comparison domain: the \
                         window's axis row (the locus's own coordinate); per \
                         fold the image is the interval spanned by the \
                         fold's walk steps shared with the axis row's walk; \
                         placements and the uniform prior are restricted to \
                         the image; out-of-image mass enters at the \
                         candidate-independent E branch for both candidates \
                         equally (rule 8)",
            });
            // The tied-winner tie certificate (the cluster machinery's
            // union-find input for tied called classes — the committed
            // receipts never carried a tied winner, so the checker's
            // phase 6 could re-derive k with an empty band list; the
            // normalized convention's converted loci DO tie (the
            // identical-LL called set), so the certificate is emitted
            // for the checker's faithful re-derivation. Territory
            // receipts only — the default receipt schema stays
            // byte-identical. THE CERTIFICATE FORM (the per-pair
            // constant's fix, mirroring the committed degenerate
            // precedent): the compact CUT-INDEX list at EVERY tie
            // shape — each tied class's spectrum index and its band's
            // exact <=cut predicate set, the only part of the band
            // the cluster machinery consumes (the naive band matrix
            // at the fully degenerate monsters is 35,510 x 35,510
            // distances, a ~6GB emission that breached the 64GiB
            // guard; at the non-degenerate fold-heavy shapes it is
            // the |called| - 1 x |spectrum| per-pair constant the
            // pre-merge anatomy receipted). The pre-fix
            // non-degenerate receipts' qual_tied_bands member is
            // replaced by this member — the answers are unchanged
            // BY CONSTRUCTION (the oracle tests).)
            // (THE CERTIFICATE STREAMS - the receipts-phase memory
            // audit's Finding 1. This member is the one structure in a
            // locus whose cardinality is (|called| - 1) x |spectrum|;
            // at the structurally degenerate shape (called.len() ==
            // pair_count, the equality this code already relies on for
            // the certificate's own form) that is pair_count squared.
            // Materializing it as a serde_json value tree cost ~32 B per
            // entry ON TOP of the typed list's 8 B and the emitted
            // text's ~6.5 B - three simultaneous copies, tens of GB at
            // the pathological shape. The member is now spliced into the
            // streamed report straight from the typed lists, so the
            // Value never exists. A Null sentinel keeps the key in the
            // map, which is what pins its alphabetical emission slot.
            report_value["qual_tied_cut_indices"] = serde_json::Value::Null;
            report_buf.clear();
            write_report_with_cut_indices(
                &mut report_buf,
                &report_value,
                "qual_tied_cut_indices",
                &cut_index_lists,
            )?;
            // (The tie certificate's emission cost - the Value tree plus
            // the reserialized text - isolated against `_certificate`.)
            rss.probe(&format!("locus_{locus}_tie_emitted"))?;
        }
        // (The alignment-connected fold relaxation's report field —
        // emitted ONLY under the relaxed rule: the merge anatomy with
        // every group's pocket masses (zero at every admitted merge —
        // the safety property, checker-verifiable), the refusal
        // counters, and the truth-second correspondence receipt. The
        // default receipt schema stays byte-identical — the identity
        // gate.)
        if alignment_connected {
            let mut report_value: serde_json::Value =
                serde_json::from_slice(&report_buf)?;
            report_value["alignment_connected_folds"] = json!({
                "rule": "the pggb bubble (the connected alignment component 
                          containing the axis correspondence) + read-level 
                          indistinguishability: bit-identical per-unit 
                          likelihood columns and zero observed read mass on 
                          every variant pocket, or the fold does not happen",
                "merge_groups": merge_groups.iter().map(|group| {
                    let row = match &territory_frame {
                        Some(frame) => {
                            let (name, start, end) = frame
                                .local_row(
                                    &group.representative_row.path_name,
                                    group.representative_row.start,
                                    group.representative_row.end,
                                )
                                .unwrap_or_else(|_| (
                                    group.representative_row.path_name.clone(),
                                    group.representative_row.start,
                                    group.representative_row.end,
                                ));
                            (name, start, end)
                        }
                        None => (
                            group.representative_row.path_name.clone(),
                            group.representative_row.start,
                            group.representative_row.end,
                        ),
                    };
                    json!({
                    "fold": group.index,
                    "representative": group.representative,
                    "representative_row": {
                        "path_name": row.0,
                        "start": row.1,
                        "end": row.2,
                    },
                    "members": group.members,
                    "pocket_node_mass": group.pocket_node_mass,
                    "pocket_edge_mass": group.pocket_edge_mass,
                })}).collect::<Vec<_>>(),
                "counters": {
                    "equal_column_pairs": merge_counters.equal_column_pairs,
                    "refused_outside_bubble": merge_counters.refused_outside_bubble,
                    "refused_unshared": merge_counters.refused_unshared,
                    "refused_pocket_mass": merge_counters.refused_pocket_mass,
                    "merged_pairs": merge_counters.merged_pairs,
                },
                "truth_second_correspondence": truth_second_correspondence,
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
        }
        // (The call-core census emission — ONLY under the env
        // instrument: the four measured work classes beside the
        // committed walls, so the autopsy's arithmetic is
        // checker-readable from the receipt itself. The default
        // receipt schema stays byte-identical — the identity gate.)
        let census_seconds = |nanos: &AtomicU64| {
            nanos.load(Ordering::Relaxed) as f64 / 1.0e9
        };
        if locus_scoring_census {
            let anchor_seconds = census_seconds(&census_anchor_ns);
            let placement_seconds = census_seconds(&census_placement_ns);
            let assertion_seconds = census_seconds(&census_assertion_ns);
            let marginal_seconds = census_seconds(&census_marginal_ns);
            let anchor_calls = census_anchor_calls.load(Ordering::Relaxed);
            let placement_calls = census_placement_calls.load(Ordering::Relaxed);
            let assertion_calls = census_assertion_calls.load(Ordering::Relaxed);
            let mut report_value: serde_json::Value =
                serde_json::from_slice(&report_buf)?;
            report_value["scoring_census"] = json!({
                "note": "the pair loop's four CPU-work classes, measured \n
                         under the per-call Instant timers (env \n
                         IMPG_LOCUS_SCORING_CENSUS; the timer's own \n
                         overhead makes these upper bounds); anchor = the \n
                         per-orientation pin build (anchor_correspondence + \n
                         the merged backbone), placement = placement_score, \n
                         assertion = the direct_score exactness check (env \n
                         IMPG_SCORE_EXACTNESS_ASSERT; the production default \n
                         is OFF — the receipt/gate form re-verifies every \n
                         placement; 0.0s/0 calls when off), \n
                         marginal = the marginalization + the row build",
                "anchor_seconds": anchor_seconds,
                "placement_seconds": placement_seconds,
                "assertion_seconds": assertion_seconds,
                "marginal_seconds": marginal_seconds,
                "anchor_calls": anchor_calls,
                "placement_calls": placement_calls,
                "assertion_calls": assertion_calls,
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
            eprintln!(
                "[score] locus {locus} scoring census: anchor {anchor_seconds:.2}s \n
                 ({anchor_calls} calls), placement {placement_seconds:.2}s \n
                 ({placement_calls} calls), assertion {assertion_seconds:.2}s \n
                 ({assertion_calls} calls), marginal {marginal_seconds:.2}s \n
                 [scoring {scoring_seconds:.2}s]",
            );
        }
        // (THE LOCUS-PHASE CENSUS's own receipt member — ONLY under
        // the env: the per-phase walls beside the shape counters
        // that price them, so the per-pair constants read off the
        // receipt itself (the fold anatomy's OPEN term's own
        // instrument). The default receipt schema stays
        // byte-identical.)
        if locus_phase_census {
            let mut report_value: serde_json::Value = serde_json::from_slice(&report_buf)?;
            report_value["phase_census"] = json!({
                "note": "the per-locus phase walls and their shape \
                         counters (env IMPG_LOCUS_PHASE_CENSUS; the \
                         honest rung is the serial call): merge = the \
                         alignment-connected fold-relaxation pass; \
                         classes = the class scan (pairs x units is \
                         the unit-pair model's term); qual_spectrum = \
                         the winner's observed-mass spectrum (one \
                         merged-multiset pair of walks per class \
                         pair); qual_certificate = the tie \
                         certificate (the cut-index-list form: \
                         per-fold mass signatures once, O(1) \
                         two-sided bounds per tied-class x \
                         spectrum-entry, the exact walk as the \
                         undecided-band fallback) — the pre-fix \
                         receipts' qual_bands_seconds (the naive \
                         |called|-1 x |spectrum| band matrix at \
                         non-degenerate shapes) and \
                         qual_degenerate_seconds are this same wall; \
                         cluster = the cluster form; best_outside = \
                         the bounded scan's best-outside fallback",
                "merge_seconds": merge_seconds,
                "merge_equal_column_pairs": merge_counters.equal_column_pairs,
                "merge_admitted_pairs": merge_counters.merged_pairs,
                "class_pairs": pair_count,
                "class_units": n_units,
                "classes_seconds": class_seconds,
                "qual_spectrum_seconds": qual_spectrum_seconds,
                "qual_spectrum_entries": spectrum.len(),
                "qual_tied_list_count": phase_census_tied_list_count,
                "qual_tied_list_entries": phase_census_tied_list_entries,
                "qual_certificate_seconds": qual_certificate_seconds,
                "qual_cluster_seconds": qual_cluster_seconds,
                "best_outside_seconds": best_outside_seconds,
                "called_class_count": called.len(),
                "degenerate": degenerate,
                "fold_multiset_entries_total": phase_census_fold_multiset_entries,
                "fold_multiset_entries_max": phase_census_fold_multiset_max,
                "scan_memo": bounded_state
                    .as_ref()
                    .and_then(|state| state.scan_memo.as_ref())
                    .or(scan_memo.as_ref())
                    .map(|memo| {
                        json!({
                            "note": "the bit-equal-row memoization (the \
                                     fourth expressibility reduction; \
                                     IMPG_CLASS_SCAN_MEMO=off restores \
                                     the direct per-pair form): folds \
                                     -> distinct bit-equal row classes; \
                                     canonical_pairs = the distinct \
                                     row-class pair sums the scan \
                                     computes; dedup = pairs served from \
                                     the memo; build = the canonicalization \
                                     walk's own wall",
                            "folds": n_folds,
                            "distinct_rows": memo.distinct_rows,
                            "deduplicated_rows": n_folds - memo.distinct_rows,
                            "canonical_pairs": memo.canonical_pairs,
                            "pairs": pair_count,
                            "dedup_pairs": pair_count as u64 - memo.canonical_pairs,
                            "computed_sums": memo.computed_sums.load(Ordering::Relaxed),
                            "memo_hits": memo.memo_hits.load(Ordering::Relaxed),
                            "build_seconds": memo.build_seconds,
                        })
                    }),
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
            let memo_line = bounded_state
                .as_ref()
                .and_then(|state| state.scan_memo.as_ref())
                .or(scan_memo.as_ref())
                .map(|memo| {
                    format!(
                        ", scan memo {}folds -> {} rows ({} dedup), \
                         {} of {} pairs dedup ({:.1}%), {} sums computed / \
                         {} memo hits, build {:.2}s",
                        n_folds,
                        memo.distinct_rows,
                        n_folds - memo.distinct_rows,
                        pair_count as u64 - memo.canonical_pairs,
                        pair_count,
                        100.0 * (pair_count as u64 - memo.canonical_pairs) as f64
                            / pair_count as f64,
                        memo.computed_sums.load(Ordering::Relaxed),
                        memo.memo_hits.load(Ordering::Relaxed),
                        memo.build_seconds,
                    )
                })
                .unwrap_or_default();
            eprintln!(
                "[score] locus {locus} phase census: merge {merge_seconds:.2}s \
                 ({} eq-col pairs, {} admitted), classes {class_seconds:.2}s \
                 ({pair_count} pairs x {n_units} units), qual spectrum \
                 {qual_spectrum_seconds:.2}s ({} entries), tie certificate \
                 {qual_certificate_seconds:.2}s ({phase_census_tied_list_count} lists x \
                 {phase_census_tied_list_entries} predicate entries), cluster \
                 {qual_cluster_seconds:.2}s, \
                 best-outside {best_outside_seconds:.2}s, called {}, fold \
                 multiset entries {} total / {} max{memo_line}",
                merge_counters.equal_column_pairs,
                merge_counters.merged_pairs,
                spectrum.len(),
                called.len(),
                phase_census_fold_multiset_entries,
                phase_census_fold_multiset_max,
            );
        }
        // (THE FOLD-SET DIGEST's own receipt member — ONLY under the
        // env: FNV-1a over the FINAL fold space's ordered identities
        // (per fold, in fold order: the length as u64 LE bytes, then
        // each member row's path-name bytes, start and end as u64 LE
        // bytes) — exactly the fold_identities fields, in that order,
        // so the digest is recomputable from any committed receipt
        // and the fold-set bit-identity gate compares two digest
        // lines, never two fold trees on disk.)
        if fold_set_digest {
            let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
            let mut feed = |bytes: &[u8], digest: &mut u64| {
                for &byte in bytes {
                    *digest ^= byte as u64;
                    *digest = digest.wrapping_mul(0x100_0000_01b3);
                }
            };
            for fold in folds.iter() {
                feed(&fold.len.to_le_bytes(), &mut digest);
                for member in &fold.members {
                    feed(member.path_name.as_bytes(), &mut digest);
                    feed(&member.start.to_le_bytes(), &mut digest);
                    feed(&member.end.to_le_bytes(), &mut digest);
                }
            }
            let mut report_value: serde_json::Value = serde_json::from_slice(&report_buf)?;
            report_value["fold_set_digest"] = json!({
                "folds": folds.len(),
                "digest": format!("{digest:016x}"),
                "note": "FNV-1a over the final fold space's ordered \
                         identities: per fold, in fold order, the \
                         length then each member row's (path_name, \
                         start, end), all integers u64 little-endian — \
                         the fold_identities fields in that order",
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
        }
        // (The class-scan order's own receipt member — ONLY under the
        // fast mode or the compare instrument: the accumulation mode
        // and the measured float deltas, so a fast-mode receipt
        // states what it is. The default receipt schema stays
        // byte-identical.)
        if sparse_skip || class_scan_compare {
            let mut report_value: serde_json::Value = serde_json::from_slice(&report_buf)?;
            report_value["class_scan"] = json!({
                "order": if sparse_skip { "sparse-skip" } else { "exact" },
                "note": "the accumulation-order toggle (the owner's Ruling 3): \
                         exact = the committed unit-order accumulation; \
                         sparse-skip = the constant entries (units unpinned \
                         on both folds) return as one elsewhere*mass term — \
                         the float accumulation CHANGES, deltas named here",
                "compare_pairs": compare_pairs.load(Ordering::Relaxed),
                "compare_pairs_differing": compare_diff.load(Ordering::Relaxed),
                "compare_max_ulp": compare_max_ulp.load(Ordering::Relaxed),
                "constant_unit_pairs_skipped": compare_skipped.load(Ordering::Relaxed),
                "classes_seconds": class_seconds,
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
        }
        // (THE SANCTIONED BOUNDED CLASS SCAN's own receipt member —
        // ONLY under the bounded mode (the sparse-skip precedent:
        // the exhaustive receipt schema stays byte-identical; the
        // labeled bounded form states its own skip certificate —
        // the kept/skipped counts, stage one's wall, the derived
        // threshold, and the per-pair bounds E*C + A_f + A_g, so a
        // bounded receipt can never be mistaken for an exhaustive
        // one and every skip is auditable offline.)
        if let Some(bounded) = bounded_state.as_ref() {
            let mut report_value: serde_json::Value = serde_json::from_slice(&report_buf)?;
            let known = bounded
                .computed
                .iter()
                .filter(|value| value.is_some())
                .count();
            report_value["classscan_bounded"] = json!({
                "kept_pairs": known,
                "skipped_pairs": pair_count - known,
                "stage_one_seconds": bounded.stage_one_seconds,
                "stage_one_best": bounded.best,
                "exp_underflow_point": exp_underflow_point(),
                "keep_threshold": bounded_keep_threshold(
                    bounded.best,
                    truth_log_likelihood,
                ),
                "pair_bounds": bounded.pair_bounds.clone(),
                "pair_slack": bounded.pair_slack.clone(),
                "note": "the labeled bounded class scan (the sanctioned \
                         frontier-1 mode): class_log_likelihoods is the exact \
                         committed sum at kept pairs and null at skipped ones; \
                         pair_bounds carries E*C + A_f + A_g per pair and \
                         pair_slack each bound's own derived rounding-error \
                         term (every skip comparison consumes bound + slack — \
                         the bound's f64 evaluation can round below the exact \
                         unit-order sum, the receipted 12-ulp class at T03226 \
                         locus 16); the emitted called set is rebuilt from the \
                         final class_lls (the tie-set fix); every answer (the \
                         MAP pair, the called set, the QUAL terms, the truth \
                         fields, the best-outside rival) is the exhaustive \
                         run's own (SS10.6.21, SS10.6.22)",
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
        }
        // (THE UNIT-FOLD REDUCTION's own receipt member — ONLY when
        // the plan is attached (env IMPG_UNIT_FOLD): the shape counters
        // the `[unit-fold]` stderr line states, so a folded receipt can
        // never be mistaken for the committed walk's. The default
        // receipt schema stays byte-identical, the labeled-toggle
        // doctrine's own form.)
        if let Some(plan) = matrix.unit_fold() {
            let mut report_value: serde_json::Value = serde_json::from_slice(&report_buf)?;
            report_value["unit_fold"] = json!({
                "note": "the unit-fold reduction (env IMPG_UNIT_FOLD, \
                         DEFAULT OFF): the class scan's unit dimension \
                         walked over the distinct bit-identical ll \
                         COLUMNS instead of the units, grouped by the \
                         exact f64 bits of each unit's whole column and \
                         never by a rounded value. A constant column \
                         (every fold bit-equal) contributes count_total \
                         * c to every pair identically (mix(c,c) = c \
                         exactly) and enters each pair sum as ONE scalar \
                         term — kept in the class LL, so the values stay \
                         comparable with the committed walk's; units with \
                         bit-identical columns share one term with their \
                         summed count mass, at the group's smallest unit \
                         index, so the surviving terms keep the original \
                         unit order. Exact in real arithmetic (the same \
                         addends), NOT bit-identical (a different \
                         association of the fixed left-to-right \
                         accumulation: T00048 locus 8, 65,520 of 65,535 \
                         kept pairs change bits, max 1,362 ulp, argmax \
                         unchanged by 1.5e7x)",
                "units": plan.units,
                "columns": plan.distinct_columns,
                "constant_columns": plan.constant_columns,
                "constant_units": plan.constant_units,
                "merged_units": plan.merged_units(),
                "merged_columns": plan.merged_columns(),
                "pair_terms_before": plan.units,
                "pair_terms_after": plan.pair_terms,
                "build_seconds": plan.build_seconds,
                "classes_seconds": class_seconds,
            });
            report_buf.clear();
            serde_json::to_writer(&mut report_buf, &report_value)?;
        }
        report_buf.push(b'\n');
        eprintln!(
            "[score] locus {locus}: receipts/IO (named classes + ingredients + \
             the main receipt) [{:.1}s]",
            receipts_started.elapsed().as_secs_f64(),
        );
        eprintln!(
            "[score] locus {locus} done: {n_folds} folds, {n_units} units, {pair_count} classes, \
             truth_rank {:?}, log_gap {log_gap:?}, qual {qual:?} [{locus_seconds:.1}s]",
            truth_rank,
        );
        // The strain roster of a fold (the committed QV receipts'
        // fold-identity form: the sorted distinct strain names of the
        // member rows).
        let fold_strains = |fold: &Fold| -> Vec<String> {
            let mut strains: Vec<String> = fold
                .members
                .iter()
                .map(|m| {
                    m.path_name
                        .split('#')
                        .next()
                        .unwrap_or_default()
                        .to_string()
                })
                .collect();
            strains.sort();
            strains.dedup();
            strains
        };
        // --------------------------------------- THE PRODUCT EMISSION
        // (the owner's product swap: the per-locus diplotype call —
        // the two called paths, the class identity, the model-internal
        // QUAL and the E/alternative-gap fields — truth-free; the
        // instrument receipt above is unchanged, and this record is
        // emitted ONLY when --product-calls is given).
        // (THE FRAME BACK-TRANSLATION, the global-panel call mode: the
        // product rows are emitted in the PARTITION'S OWN frame — the
        // committed reference-anchored output shape — via the inverse
        // record table; the machinery's global rows translate back
        // fail-closed. Under the default frame this is the identity.)
        if product_emission {
            let product_fold = |index: usize| -> io::Result<serde_json::Value> {
                let (image_lo, image_hi) = territory[index];
                let members: Vec<serde_json::Value> = folds[index]
                    .members
                    .iter()
                    .map(|m| match &territory_frame {
                        Some(frame) => {
                            let (name, start, end) =
                                frame.local_row(&m.path_name, m.start, m.end)?;
                            Ok(json!({
                                "path_name": name, "start": start, "end": end,
                            }))
                        }
                        None => Ok(json!({
                            "path_name": m.path_name, "start": m.start, "end": m.end,
                        })),
                    })
                    .collect::<io::Result<Vec<_>>>()?;
                // (The territory image is FOLD-RELATIVE — the same
                // extent in either frame; no translation.)
                Ok(json!({
                    "fold_index": index,
                    "members": members,
                    "length": folds[index].len,
                    "territory_image": if image_hi > image_lo {
                        json!([image_lo, image_hi])
                    } else {
                        serde_json::Value::Null
                    },
                }))
            };
            let (winner_first, winner_second) = winner_pair;
            // THE MODEL-INTERNAL QUAL (the owner's cluster semantics,
            // derived and stated): k = the cluster count (the winner's
            // connected component under the knee cut plus every
            // exactly-tied called class's own component); the
            // alternative = the best relative likelihood e^{-Δ} over
            // the classes OUTSIDE the excluded band (Δ = the
            // alternative log-gap); the two-hypothesis flat-prior
            // posterior p = 1/(k + e^{-Δ}); QUAL = -10 log10(1 - p).
            // Where the f64 posterior is representable (p < 1) the
            // emitted value is the committed cluster form verbatim
            // (field-identical to the receipts); where it saturated
            // to exactly 1.0 (the receipts' qual null +
            // qual_unbounded) the SAME form is evaluated in log
            // space — QUAL = 10/ln(10) * ln((k + e^{-Δ}) /
            // (k - 1 + e^{-Δ})) — so a realignment separation of Δ
            // yields QUAL = 4.343 * Δ Phred, growing linearly with
            // the log-gap and never saturating meaninglessly (k > 1
            // saturates only at the honest cap 10*log10(k/(k-1)) of
            // exactly-tied rivals).
            let product_qual = match (qual, qual_unbounded) {
                (Some(value), _) => json!(value),
                (None, true) => match (cluster_k, alternative_log_gap) {
                    // k = 1 (the only unbounded case): QUAL =
                    // 10/ln(10) * (Δ + ln(1 + e^{-Δ})) — the stable
                    // log-space form (the naive ln(e^{-Δ}) underflows
                    // to -inf at deep separations; the exact value is
                    // 4.343 * Δ Phred + the vanishing tie bonus).
                    (Some(1), Some(delta)) => json!(
                        10.0 * (delta + (-delta).exp().ln_1p())
                            / std::f64::consts::LN_10
                    ),
                    (Some(k), Some(delta)) if k >= 2 => {
                        let alternative = (-delta).exp();
                        json!(10.0
                            * ((k as f64 + alternative).ln()
                                - (k as f64 - 1.0 + alternative).ln())
                            / std::f64::consts::LN_10)
                    }
                    _ => serde_json::Value::Null,
                },
                (None, false) => serde_json::Value::Null,
            };
            let mut product_value = json!({
                    "locus": locus,
                    "partition": partition,
                    "component": options.component,
                    "model": model_name,
                    "diplotype": {
                        "fold_indices": [winner_first, winner_second],
                        "paths": [product_fold(winner_first)?, product_fold(winner_second)?],
                    },
                    "class_count": pair_count,
                    "called_class_count": called.len(),
                    "called_classes": called.iter()
                        .map(|&(first, second)| json!([first, second]))
                        .collect::<Vec<_>>(),
                    "unit_count": n_units,
                    "record_count": locus_records.len(),
                    "best_log_likelihood": best_ll,
                    "elsewhere_log_prob": elsewhere,
                    "alternative_log_gap": alternative_log_gap,
                    "qual": product_qual,
                    "qual_unbounded": qual_unbounded,
                    "qual_spectrum_shape": cluster_shape,
                    "qual_knee_distance": cluster_knee,
                    "qual_cluster_size": cluster_size,
                    "qual_cluster_k": cluster_k,
                    "qual_alternative_similarity": cluster_alternative,
                    "qual_alternative_fold_pairs": qual_alternative_classes.iter()
                        .map(|value| value["fold_indices"].clone())
                        .collect::<Vec<_>>(),
                    "nearest_rival_distance": nearest_rival_distance,
            });
            if excluded_rows > 0 {
                if let Some(object) = product_value.as_object_mut() {
                    object.insert("excluded_rows".to_string(), json!(excluded_rows));
                    object.insert(
                        "exclude_paths".to_string(),
                        json!(options.exclude_paths),
                    );
                }
                eprintln!(
                    "[score] locus {locus}: {excluded_rows} candidate rows excluded by \
                     --exclude-paths (the stranger filter; the axis row kept as the frame)"
                );
            }
            if flank_excluded_rows > 0 {
                if let Some(object) = product_value.as_object_mut() {
                    object.insert(
                        "flank_excluded_rows".to_string(),
                        json!(flank_excluded_rows),
                    );
                    object.insert(
                        "flank_overlap".to_string(),
                        json!(options
                            .flank_overlap
                            .map(|f| f.setting())
                            .unwrap_or_default()),
                    );
                }
                eprintln!(
                    "[score] locus {locus} (partition {partition}): {flank_excluded_rows} \
                     fragmentary candidate rows excluded by --flank-overlap {} \
                     (missing a flank of [{axis_start}, {axis_end}); the axis row kept \
                     as the frame)",
                    options
                        .flank_overlap
                        .map(|f| f.setting())
                        .unwrap_or_default()
                );
            }
            serde_json::to_writer(&mut product_buf, &product_value)?;
            // The receipt's own serialization convention, mirrored
            // exactly: under the territory/locality conventions the
            // receipt is reparsed through serde_json::from_slice
            // before its final emission (the committed double-buffer
            // discipline), and serde_json's float parse is up to 1
            // ULP from the correctly-rounded value at some decimals —
            // so the product record passes through the IDENTICAL
            // round trip under the same conventions, keeping its
            // f64 fields bit-identical to the receipt's (the
            // translation gate's field parity). Under the default
            // convention neither record is reparsed (raw values
            // both).
            if territory_normalized {
                let mut value: serde_json::Value = serde_json::from_slice(&product_buf)?;
                product_buf.clear();
                serde_json::to_writer(&mut product_buf, &value)?;
            }
            product_buf.push(b'\n');
            // ------------------ THE DIPLOTYPE SEQUENCE EMISSION
            // (the write-up codification, the owner's ruling: the
            // diplotype sequence is a PRODUCT OUTPUT at every called
            // locus): the chosen pair's walks from the local graph
            // spelled into per-locus hap1/hap2 records - the fold's
            // own spelled sequence, the SAME material the truth-QV
            // alignments consume - each record carrying the
            // model-internal QUAL, the class identity, and the tied
            // classes' coalesced called set (every called fold's
            // sequence, the owner's ambiguity ruling). Truth-free:
            // no truth name, alignment or QV field enters here.
            if options.product_sequences {
                let fold_record = |index: usize| -> io::Result<serde_json::Value> {
                    let members: Vec<serde_json::Value> = folds[index]
                        .members
                        .iter()
                        .map(|m| match &territory_frame {
                            Some(frame) => {
                                let (name, start, end) =
                                    frame.local_row(&m.path_name, m.start, m.end)?;
                                Ok(json!({
                                    "path_name": name, "start": start, "end": end,
                                }))
                            }
                            None => Ok(json!({
                                "path_name": m.path_name, "start": m.start, "end": m.end,
                            })),
                        })
                        .collect::<io::Result<Vec<_>>>()?;
                    Ok(json!({
                        "fold_index": index,
                        "members": members,
                        "length": folds[index].len,
                        "sequence": fold_sequence_string(&folds[index].seq)?,
                    }))
                };
                let mut sequences_value = json!({
                    "locus": locus,
                    "partition": partition,
                    "component": options.component,
                    "qual": product_qual,
                    "qual_unbounded": qual_unbounded,
                    "class_count": pair_count,
                    "called_class_count": called.len(),
                    "hap1": fold_record(winner_first)?,
                    "hap2": fold_record(winner_second)?,
                });
                if let Some(object) = sequences_value.as_object_mut() {
                    object.insert(
                        "called_set".to_string(),
                        coalesced_called_set(&called)
                            .into_iter()
                            .map(fold_record)
                            .collect::<io::Result<Vec<_>>>()?
                            .into(),
                    );
                }
                serde_json::to_writer(&mut product_sequences_buf, &sequences_value)?;
                // The receipt's own serialization convention, mirrored
                // from the product record: under the territory/locality
                // conventions the record is reparsed through
                // serde_json::from_slice before its final emission (the
                // committed double-buffer discipline), so the emitted
                // QUAL floats are bit-identical to the product call's.
                if territory_normalized {
                    let mut value: serde_json::Value =
                        serde_json::from_slice(&product_sequences_buf)?;
                    product_sequences_buf.clear();
                    serde_json::to_writer(&mut product_sequences_buf, &value)?;
                }
                product_sequences_buf.push(b'\n');
            }
            // ------------------------------- THE TRUTH-QV TEST MODE
            // (the owner's ruling: the called-vs-truth sequence QV
            // ONLY here, behind --truth-qv with the truth haplotypes
            // as input — a separate test artifact, never mixed into
            // the product output. The committed convention verbatim:
            // the best no-replacement assignment of the two called
            // homologs to the two truth homologs, biWFA gap-affine
            // End2End 0/4/6/2, error = (mismatches + gap columns) /
            // total columns, QV = -10 log10(error) with the QV-60 cap
            // at error <= 1e-6.)
            let truth_qv_value = match truth_pair_test {
                    None => json!({
                        "locus": locus,
                        "partition": partition,
                        "component": options.component,
                        "truth_pair_expressible": false,
                        "truth_second_correspondence":
                            truth_second_correspondence_test,
                        "called_folds": [
                            {"index": winner_first,
                             "strains": fold_strains(&folds[winner_first])},
                            {"index": winner_second,
                             "strains": fold_strains(&folds[winner_second])},
                        ],
                        "qv": serde_json::Value::Null,
                        "perfect": false,
                    }),
                    Some((truth_first_index, truth_second_index)) => {
                        let pair_counts = |called_index: usize, truth_index: usize| {
                            let Some(aligner) = truth_qv_aligner.as_ref() else {
                                return Err(invalid(
                                    "the truth-QV aligner is not available",
                                ));
                            };
                            let [matches, mismatches, ins, dels] = aligner
                                .lock()
                                .map_err(|_| invalid("the truth-QV aligner lock poisoned"))?
                                .align(&folds[called_index].seq, &folds[truth_index].seq)?;
                            let (edits, columns) = (mismatches + ins + dels,
                                matches + mismatches + ins + dels);
                            Ok(json!({
                                "matches": matches, "mismatches": mismatches,
                                "ins": ins, "dels": dels,
                                "columns": columns, "edits": edits,
                            }))
                        };
                        let p00 = pair_counts(winner_first, truth_first_index)?;
                        let p11 = pair_counts(winner_second, truth_second_index)?;
                        let p01 = pair_counts(winner_first, truth_second_index)?;
                        let p10 = pair_counts(winner_second, truth_first_index)?;
                        let order_score = |first: &serde_json::Value,
                                           second: &serde_json::Value| {
                            let edits = first["edits"].as_u64().unwrap_or(0)
                                + second["edits"].as_u64().unwrap_or(0);
                            let columns = first["columns"].as_u64().unwrap_or(0)
                                + second["columns"].as_u64().unwrap_or(0);
                            let error = if columns > 0 {
                                edits as f64 / columns as f64
                            } else {
                                0.0
                            };
                            (error, edits, columns)
                        };
                        let direct = order_score(&p00, &p11);
                        let crossed = order_score(&p01, &p10);
                        // The owner's yardstick: the lower rate wins;
                        // ties to the lower total edits, then direct.
                        let (assignment, (error, edits, columns)) =
                            if (direct.0, direct.1, 0)
                                <= (crossed.0, crossed.1, 1)
                            {
                                ("direct", direct)
                            } else {
                                ("crossed", crossed)
                            };
                        // likegt's own rule (the higher mean pair
                        // identity) — disagreement reported, never
                        // hidden.
                        let mean_identity = |first: &serde_json::Value,
                                             second: &serde_json::Value| {
                            let one = |pair: &serde_json::Value| {
                                let (edits, columns) = (
                                    pair["edits"].as_u64().unwrap_or(0),
                                    pair["columns"].as_u64().unwrap_or(0),
                                );
                                    if columns > 0 {
                                        1.0 - edits as f64 / columns as f64
                                    } else {
                                        1.0
                                    }
                            };
                            (one(first) + one(second)) / 2.0
                        };
                        let likegt_assignment = if mean_identity(&p00, &p11)
                            >= mean_identity(&p01, &p10)
                        {
                            "direct"
                        } else {
                            "crossed"
                        };
                        let (chosen_first, chosen_second) =
                            if assignment == "direct" { (&p00, &p11) } else { (&p01, &p10) };
                        let qv = if error > 0.0 {
                            (-10.0 * error.max(1e-6).log10()).min(60.0)
                        } else {
                            60.0
                        };
                        let mut record = json!({
                            "locus": locus,
                            "partition": partition,
                            "component": options.component,
                            "truth_pair_expressible": true,
                            "truth_second_correspondence":
                                truth_second_correspondence_test,
                            "truth_rank": truth_rank_test,
                            "rank1": truth_rank_test == Some(1),
                            "called_folds": [
                                {"index": winner_first,
                                 "strains": fold_strains(&folds[winner_first])},
                                {"index": winner_second,
                                 "strains": fold_strains(&folds[winner_second])},
                            ],
                            "truth_folds": [
                                {"index": truth_first_index,
                                 "strains": fold_strains(&folds[truth_first_index])},
                                {"index": truth_second_index,
                                 "strains": fold_strains(&folds[truth_second_index])},
                            ],
                            "assignment": assignment,
                            "likegt_assignment": likegt_assignment,
                            "assignment_disagreement": assignment != likegt_assignment,
                            "pairs": {
                                "called0_truth0": p00,
                                "called1_truth1": p11,
                                "called0_truth1": p01,
                                "called1_truth0": p10,
                            },
                            "chosen_pairs": [chosen_first, chosen_second],
                            "total_edits": edits,
                            "total_columns": columns,
                            "mismatches": chosen_first["mismatches"].as_u64().unwrap_or(0)
                                + chosen_second["mismatches"].as_u64().unwrap_or(0),
                            "gap_columns": chosen_first["ins"].as_u64().unwrap_or(0)
                                + chosen_first["dels"].as_u64().unwrap_or(0)
                                + chosen_second["ins"].as_u64().unwrap_or(0)
                                + chosen_second["dels"].as_u64().unwrap_or(0),
                            "error": error,
                            "identity": 1.0 - error,
                            "qv": qv,
                            "perfect": edits == 0,
                        });
                        // THE OFFICIAL SEQUENCE-QV CONVENTION (the
                        // owner's ruling, the likegt max-qv form),
                        // env-gated beside the committed biWFA fields:
                        // allwave aligns the two truth haplotypes
                        // against the two called homologs, per-hap
                        // best selection, QV from the best
                        // alignments' per-base error (identity =
                        // matches / query_span, the QV-50 cap), the
                        // diplotype QV the mean of the two per-hap
                        // QVs. Additive only: the committed fields
                        // above are byte-identical.
                        if let Some(allwave) = truth_qv_allwave.as_ref() {
                            let official = allwave.locus_qv(
                                locus,
                                &folds[winner_first].seq,
                                &folds[winner_second].seq,
                                &folds[truth_first_index].seq,
                                &folds[truth_second_index].seq,
                            )?;
                            if let Some(object) = record.as_object_mut() {
                                object.insert("official_allwave".to_string(), official);
                            }
                        }
                        record
                    }
                };
            serde_json::to_writer(&mut truth_qv_buf, &truth_qv_value)?;
            truth_qv_buf.push(b'\n');
        }
        Ok(LocusOutput {
            skeleton: skeleton_buf,
            exactness: exactness_buf,
            records: records_buf,
            ingredients_spill: ingredients_spill_path,
            report: report_buf,
            anatomy: anatomy_buf,
            product: product_buf,
            product_sequences: product_sequences_buf,
            truth_qv: truth_qv_buf,
            factorized_checked,
            factorized_equal,
            locus_seconds,
        })
    };
    let outputs: Vec<LocusOutput> = match if seam_parallel {
        seam_map(&loci, seam_width, &|&locus| compute_locus(locus))
    } else {
        loci
            .iter()
            .map(|&locus| compute_locus(locus))
            .collect::<io::Result<Vec<_>>>()
    } {
        Ok(outputs) => outputs,
        Err(error) => {
            // (A failed locus may have left its ingredients spill
            // behind; the run is failing — remove the spills.)
            for &locus in &loci {
                let _ =
                    std::fs::remove_file(ingredients_spill_path(&options.out, locus));
            }
            return Err(error);
        }
    };
    // The ordered emission: the per-locus buffers are written in LOCUS
    // ORDER — the receipt line order of the committed serial runs.
    for output in outputs {
        if let Some(writer) = skeleton.as_mut() {
            writer.write_all(&output.skeleton)?;
        }
        exactness.write_all(&output.exactness)?;
        records_file.write_all(&output.records)?;
        // (Candidate C: stream the locus's ingredients spill into the
        // ordered file — identical bytes to the committed buffered
        // write — then remove the spill.)
        let spill = File::open(&output.ingredients_spill)?;
        std::io::copy(&mut BufReader::new(spill), &mut ingredients)?;
        std::fs::remove_file(&output.ingredients_spill)?;
        report.write_all(&output.report)?;
        if let Some(writer) = product_calls.as_mut() {
            writer.write_all(&output.product)?;
        }
        if let Some(writer) = product_sequences_out.as_mut() {
            writer.write_all(&output.product_sequences)?;
        }
        if let Some(writer) = truth_qv_out.as_mut() {
            writer.write_all(&output.truth_qv)?;
        }
        if let Some(writer) = anatomy.as_mut() {
            writer.write_all(&output.anatomy)?;
        }
        total_factorized_checked += output.factorized_checked;
        total_factorized_equal += output.factorized_equal;
        loci_seconds_total += output.locus_seconds;
    }

    report.flush()?;
    exactness.flush()?;
    ingredients.flush()?;
    records_file.flush()?;
    if let Some(writer) = product_calls.as_mut() {
        writer.flush()?;
    }
    if let Some(writer) = product_sequences_out.as_mut() {
        writer.flush()?;
    }
    if let Some(writer) = truth_qv_out.as_mut() {
        writer.flush()?;
    }
    if let Some(writer) = anatomy.as_mut() {
        writer.flush()?;
    }
    ensure(
        total_factorized_checked == total_factorized_equal,
        "the factorization gate failed somewhere",
    )?;
    eprintln!(
        "[score] complete: loci {loci:?}, factorized placements {total_factorized_equal}/\
         {total_factorized_checked} exact, total {:.1}s",
        started.elapsed().as_secs_f64(),
    );
    eprintln!(
        "[score] lane cache: {} lanes loaded ({:.1} MB, {:.1}s of load time)",
        lane_cache.loads.load(Ordering::Relaxed),
        lane_cache.bytes.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
        lane_cache.nanos.load(Ordering::Relaxed) as f64 / 1e9,
    );
    eprintln!(
        "[score] phase summary: inputs {inputs_seconds:.1}s, quality (FASTQ scan) \
         {quality_seconds:.1}s, elsewhere {elsewhere_seconds:.3}s, derive cache \
         {derive_seconds:.1}s, binding {binding_seconds:.1}s, variants \
         {variants_seconds:.1}s, rows/context assembly {rows_seconds:.1}s, \
         loci (folds+units+scoring+classes+qual+exactness+receipts) \
         {loci_seconds_total:.1}s; AGC fetches {} total ({:.1} MB)",
        fetch_calls.load(Ordering::Relaxed),
        fetch_bytes.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
    );
    if let Some(aligner) = truth_qv_aligner.as_ref() {
        aligner
            .lock()
            .map_err(|_| invalid("the truth-QV aligner lock poisoned"))?
            .finish()?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests: the factorized exact form against the direct per-base walk,
// the strand symmetry, the projection equations, and the machinery.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const K: u64 = 8;

    fn fold_seq() -> Vec<u8> {
        // 60 bp, deterministic and non-repetitive enough for the tests.
        (0..60)
            .map(|i| match i % 5 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b"ACGT"[(i / 5) % 4],
            })
            .collect()
    }

    /// The canonical pattern: anchors (node 5 @ 0), (node 7 @ 10);
    /// span 18; the fold carries node 5 forward at rel 4 and node 7
    /// forward at rel 14. The forward read F (30 bp, w_lo 2) spells the
    /// fold from rel 2 to 32 with ONE deliberate mismatch in the
    /// 2 bp interior hull gap.
    fn forward_read(seq: &[u8]) -> (Vec<u8>, usize) {
        let mut read = vec![0u8; 30];
        for i in 0..30usize {
            read[i] = seq[i + 2];
        }
        read[11] = match read[11] {
            b'A' => b'T',
            other => b'A' + (other == b'A') as u8,
        };
        (read, 13) // the mismatching fold coordinate (rel 13)
    }

    fn pin() -> Pin {
        Pin {
            orientation: 0,
            pinned: vec![Some(4), Some(14)],
            backbone_bp: 16,
        }
    }

    fn canonical() -> Vec<(i32, u64)> {
        vec![(5, 0), (7, 10)]
    }


    #[test]
    fn factorized_equals_direct_on_the_forward_placement() {
        let seq = fold_seq();
        let (read, mismatch_coord) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 60)
            .expect("valid placement");
        // Backbone 16 matches + 14 skipped votes, one of them the
        // deliberate mismatch at fold rel 13.
        assert_eq!((m, c), (29, 1));
        let (dm, dc) = direct_score(&unit, &canonical(), &pin, &seq, 60, &reads, K).unwrap();
        assert_eq!((m, c), (dm, dc));
        // The mismatch is where we placed it.
        let canonical = canonical();
        let own_positions = own_positions_of(&canonical, false, 2);
        let skipped = skipped_positions(&own_positions, K, 30);
        assert_eq!(skipped, vec![(0, 2), (10, 12), (20, 30)]);
        let fold_base = seq[mismatch_coord];
        let read_base = reads[0][11];
        assert_ne!(read_base, fold_base);
    }

    #[test]
    fn the_reverse_strand_placement_is_strand_symmetric() {
        // M = revcomp(F): the same physical read; its own orientation
        // extraction sees the pattern mirrored (w_lo 10), and at the
        // same orientation-0 pin the placement is physically reverse.
        // The per-base evidence must be IDENTICAL: (m, c) equal to the
        // forward case, every voted base the same physical comparison.
        let seq = fold_seq();
        let (forward, _) = forward_read(&seq);
        let mirrored = revcomp(&forward);
        let reads = vec![forward, mirrored];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: true,
            read: 1,
            w_lo: 10,
            count: 1,
        };
        // Mirror own positions: canonical anchor 0 at hull 10, anchor 1
        // at hull 0 — the correspondence is canonical-order-based and
        // mirror-blind, so the pin is the same.
        let own = own_positions_of(&canonical(), true, 10);
        assert_eq!(own, vec![20, 10]);
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 60)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (dm, dc) = direct_score(&unit, &canonical(), &pin, &seq, 60, &reads, K).unwrap();
        assert_eq!((m, c), (dm, dc));
    }

    #[test]
    fn placements_with_out_of_extent_votes_are_invalid() {
        // A fold whose extent ends inside the read's right flank: the
        // voted bases do not all project inside, so the placement is
        // not a placement of the model (the read must lie fully
        // inside).
        let seq = fold_seq();
        let (read, _) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = Pin {
            orientation: 0,
            pinned: vec![Some(4), Some(14)],
            backbone_bp: 16,
        };
        let placed = placement_score(&unit, &canonical(), &pin, &seq, 24, &reads, K, 0, 24);
        assert!(placed.is_none());
        // Backbone windows partly outside the extent are not
        // placements either (the windows_inside check at pin
        // construction discards them).
        let pin_outside = [Some(196i64), Some(206)];
        assert!(pin_outside.iter().flatten().any(|&r| r + K as i64 > 24));
    }

    #[test]
    fn the_territory_image_restricts_placements_symmetrically() {
        // Stage 2, rule 8: under the normalized convention the
        // comparison domain is the fold's TERRITORY IMAGE — a
        // placement whose read span exits the image is not a local
        // placement of the locus, exactly as one exiting the fold's
        // extent is not a placement under the committed convention.
        // The SAME placement that is valid at the full extent
        // [0, 60) must score IDENTICALLY inside a containing image
        // (the degenerate default), and must be INVALID once the
        // image ends before the read's voted bases: the rule
        // restricts every candidate by its OWN image of the SAME
        // territory — it cannot collect out-of-image mass.
        let seq = fold_seq();
        let (read, _) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = pin();
        // The containing image: identical score (the degenerate
        // default convention, bit-for-bit).
        let (m, c) =
            placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 60)
                .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (m2, c2) =
            placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 34)
                .expect("the read span [2, 32) lies inside the image");
        assert_eq!((m2, c2), (29, 1));
        // The image ending at 31 clips the read's last voted base
        // (fold rel 32): not a placement — the read is not generated
        // wholly within this fold's image of the territory.
        assert!(
            placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 31).is_none()
        );
        // An empty image admits no placement at all (a fold sharing
        // no anchor with the axis row spells no territory material).
        assert!(
            placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 0).is_none()
        );
        // The prior counts the IMAGE's placements under the
        // normalized convention (the support length, not the fold's
        // full extent) — and an image shorter than the read has no
        // support at all.
        let a = 0.0f64;
        let full = unit_log_likelihood(&[a], 200, READ_LENGTH_FRESH);
        let image = unit_log_likelihood(&[a], 150, READ_LENGTH_FRESH);
        assert!(full < image); // the smaller support carries the larger prior
        assert_eq!(
            unit_log_likelihood(&[a], 149, READ_LENGTH_FRESH),
            f64::NEG_INFINITY,
        );
    }

    #[test]
    fn the_backbone_is_mirror_invariant() {
        // The merged pinned-window coverage in hull coordinates has
        // the same total length under both mirror states (the mirrored
        // window set is the reflection of the forward set).
        let canonical = vec![(5, 0), (7, 10), (-9, 14)];
        for mirror in [false, true] {
            let own = own_positions_of(&canonical, mirror, 0);
            let windows: Vec<(u64, u64)> = own.iter().map(|&p| (p, p + K)).collect();
            let total: u64 = merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
            assert_eq!(total, 20); // windows [0,8),[10,18),[14,22) -> merged [0,8)+[10,22)
        }
    }

    #[test]
    fn unit_likelihood_is_prior_plus_logsumexp() {
        // No placements: the anchored floor.
        assert_eq!(
            unit_log_likelihood(&[], 200, READ_LENGTH_FRESH),
            f64::NEG_INFINITY,
        );
        assert_eq!(
            unit_log_likelihood(&[0.0], 100, READ_LENGTH_FRESH),
            f64::NEG_INFINITY,
        );
        let a = -3.0f64;
        let b = -5.0f64;
        let ll = unit_log_likelihood(&[a, b], 200, READ_LENGTH_FRESH);
        let expect = -(2.0 * (200 - 149) as f64).ln() + (a.exp() + b.exp()).ln();
        assert!((ll - expect).abs() < 1e-12);
        // Bit-exact when the placements are equal (the identical-fold
        // case): the sum is the placement plus ln 2.
        let ll = unit_log_likelihood(&[a, a], 200, READ_LENGTH_FRESH);
        let expect = -(2.0 * (200 - 149) as f64).ln() + (2.0 * a.exp()).ln();
        assert!((ll - expect).abs() < 1e-12);
    }

    #[test]
    fn the_serving_anchor_brackets_the_skipped_position() {
        // Pinned anchors at hull positions 0 and 10: the left flank is
        // served by the LEAST, the interior gap and the right flank by
        // the GREATEST pinned anchor at or before the base.
        let pinned = vec![Some(4i64), Some(14)];
        let own_hull = vec![0u64, 10];
        // No pinned anchor at or before d: the LEAST serves.
        assert_eq!(serving_anchor(&pinned, &own_hull, -2), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 0), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 8), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 10), Some((1, 10)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 25), Some((1, 10)));
    }

    #[test]
    fn unrank_pair_is_the_inverse_of_flat() {
        for i in 0..40usize {
            for j in i..40usize {
                let flat = j * (j + 1) / 2 + i;
                assert_eq!(unrank_pair(flat), (i, j));
            }
        }
    }

    #[test]
    fn collinearity_holds_for_true_placements_both_strands() {
        // Forward: r_j - own_hull_j invariant.
        let pinned = vec![Some(4i64), Some(14)];
        let own_forward = vec![0u64, 10];
        let (collinear, offset) = placement_offset(&pinned, &own_forward, true);
        assert!(collinear);
        assert_eq!(offset, Some(4));
        // Physically reverse: the reflection r_j + own_hull_j invariant
        // (the mirrored read's own hull positions are reflected).
        let own_mirror = vec![10u64, 0];
        let (collinear, offset) = placement_offset(&pinned, &own_mirror, false);
        assert!(collinear);
        assert_eq!(offset, Some(14));
        // An indel between read and row pins the anchors at mutually
        // incompatible offsets: NOT one continuous path.
        let pinned_indel = vec![Some(4i64), Some(20)];
        let (collinear, _) = placement_offset(&pinned_indel, &own_forward, true);
        assert!(!collinear);
    }

    #[test]
    fn the_one_continuous_path_score_equals_the_placement_score_when_collinear() {
        // The slice-B piecewise (m, c) and the whole-read single-offset
        // score agree exactly on a collinear placement where every
        // anchor pins: the piecewise serving projects the flanks at
        // the same offset, and the previously-abstaining within-window
        // bases agree because the read is a contiguous segment.
        let seq = fold_seq();
        let (read, _) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 60)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (sm, sc, oob) = single_offset_score(&unit, 4, true, &seq, 60, &reads, K);
        assert_eq!((sm, sc, oob), (29, 1, 0));
        // The mirrored read at the same pin: the reflected offset
        // const 14 gives the identical per-base evidence.
        let (forward_read_bytes, _) = forward_read(&seq);
        let mirrored = revcomp(&forward_read_bytes);
        let reads = vec![forward_read_bytes, mirrored];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: true,
            read: 1,
            w_lo: 10,
            count: 1,
        };
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K, 0, 60)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (sm, sc, oob) = single_offset_score(&unit, 14, false, &seq, 60, &reads, K);
        assert_eq!((sm, sc, oob), (29, 1, 0));
    }

    #[test]
    fn the_diploid_mixture_handles_the_floor() {
        assert_eq!(mix_logsumexp(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
        let a = -2.0f64;
        // The homozygous class: the mixture of a fold with itself is
        // the fold's own likelihood (0.5 P + 0.5 P = P).
        assert!((mix_logsumexp(a, a) - a).abs() < 1e-12);
        assert!((mix_logsumexp(a, f64::NEG_INFINITY) - (a - (2f64).ln())).abs() < 1e-12);
        let mixed = mix_logsumexp(-2.0, -3.0);
        assert!((mixed - ((-2.0f64).exp() * 0.5 + (-3.0f64).exp() * 0.5).ln()).abs() < 1e-12);
    }

    #[test]
    fn the_knee_cut_is_the_max_relative_jump_above_the_typical() {
        // A spectrum with a clean knee: distances 1,1,2,2,3,3,90.
        let sorted = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 90.0];
        let knee = spectrum_knee(&sorted);
        assert!(knee.has_knee);
        assert_eq!(knee.cut, 3.0);
        // No knee in a geometric spectrum (every relative jump equal
        // to the typical): the derived no-knee case.
        let geometric = vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];
        assert!(!spectrum_knee(&geometric).has_knee);
    }

    #[test]
    fn the_elsewhere_branch_is_derived_and_bounds_the_floor() {
        // E(read) = len*A - ln(2*(G - len + 1)): the measured per-base
        // match rate under the max-entropy uniform origin over the
        // genome's both-strand full-read placements.
        let a = (1.0 - 1e-4f64).ln();
        let e = elsewhere_log_prob(150, a, 24_400_000.0);
        let expect = 150.0 * a - (2.0f64 * (24_400_000.0 - 150.0 + 1.0)).ln();
        assert!((e - expect).abs() < 1e-9);
        // The floor for an unplaceable read is E — far above the
        // all-mismatch minimum, and the mixture with a missing local
        // branch (NEG_INFINITY) is exactly E.
        assert_eq!(logsumexp2(f64::NEG_INFINITY, e), e);
        // A local branch scoring below E is absorbed by it: the read
        // does not vote beyond E (a full-mismatch or heavily-mismatched
        // placement carries no credit).
        let bad_local = -(2.0f64 * (10_000.0 - 149.0)).ln() + 140.0 * a + 10.0 * (1e-4f64 / 3.0).ln();
        assert!(logsumexp2(bad_local, e) < e + 1e-6);
        // A full-match local placement on a 10kb row survives and
        // beats E by the prior ratio (the locale's concentration).
        let good_local = -(2.0f64 * (10_000.0 - 149.0)).ln() + 150.0 * a;
        assert!(good_local > e);
        assert!(logsumexp2(good_local, e) > e);
    }

    #[test]
    fn the_marginalization_is_monotone_and_preserves_ties() {
        // The L4-control preservation property: logsumexp2 is monotone
        // in the local branch, so a candidate whose local branch is
        // at least every rival's on every unit keeps its per-unit
        // ordering after the marginalization; identical local branches
        // stay identical (the identical-fold non-votes stay non-votes).
        let e = elsewhere_log_prob(150, (1.0 - 1e-4f64).ln(), 24_400_000.0);
        let a_local = -9.9f64;
        let b_local = -25.0f64;
        assert!(logsumexp2(a_local, e) > logsumexp2(b_local, e));
        assert_eq!(logsumexp2(a_local, e), logsumexp2(a_local, e));
        // Deep below E the branches are indistinguishable: the read
        // does not vote.
        let x = logsumexp2(-1546.0f64, e);
        let y = logsumexp2(-2000.0f64, e);
        assert!((x - y).abs() < 1e-12);
        assert!((x - e).abs() < 1e-12);
    }

    #[test]
    fn qual_is_the_cluster_form_over_the_excluded_band() {
        // One divergent rival far outside a tight band: the knee cuts
        // the band, k = 2, p = 1/(2 + a).
        let spectrum: Vec<(f64, f64)> = vec![
            (0.0, 1.0),
            (1.0, 0.5),
            (2.0, 0.1),
            (100.0, 1e-40),
        ];
        let state = cluster_form_qual(1.0, &spectrum, &TieCertificate::Bands(&[]));
        assert!(state.knee.is_some());
        // A single called class: k = 1, and p = 1/(1 + a) is capped by
        // the alternative's likelihood mass.
        assert_eq!(state.k, 1);
        assert_eq!(state.alternative, Some(1e-40));
        let p = state.p.expect("finite p");
        assert!((p - 1.0 / (1.0 + 1e-40)).abs() < 1e-12);
    }

    #[test]
    fn the_observed_mass_arrays_match_the_hash_probe_walk() {
        // THE ARRAY WALK'S EXACTNESS LEMMA (the residual fix's gate):
        // `differing_observed_mass_from_units` must be the committed
        // `differing_observed_mass` BIT FOR BIT on the same inputs,
        // including the equal branch, whose product multiplies the unit
        // observed mass by `delta` and not by `count` — the branch a
        // stored product would break, and the branch most undecided
        // ties sit on (near-identical shapes share most of their keys).
        // The observed map carries hostile values (non-representable
        // decimals, 1/3, a magnitude that annihilates small addends) so
        // a reordered accumulation or a stored product would show up in
        // the bits.
        let mut observed: HashMap<u64, f64> = HashMap::new();
        for key in 0..40u64 {
            observed.insert(key, 0.1 * (key as f64 + 1.0) + key as f64 / 3.0);
        }
        observed.insert(9, 1e16);
        observed.insert(17, 1.0 / 3.0);
        let lookup = |key: u64| observed.get(&key).copied().unwrap_or(0.0);

        // the plain lemma, hand-written: key 2 equal at DIFFERENT counts
        // (the delta product), key 7 equal at equal counts (skipped),
        // key 11 left-only and key 13 right-only (the two tails)
        let left = vec![(2u64, 5u32), (7, 3), (11, 2)];
        let right = vec![(2u64, 1u32), (7, 3), (13, 4)];
        let left_observed: Vec<f64> = left.iter().map(|&(key, _)| lookup(key)).collect();
        let right_observed: Vec<f64> =
            right.iter().map(|&(key, _)| lookup(key)).collect();
        let expected = differing_observed_mass(&left, &right, &lookup);
        // (the accumulation is the committed one: observed(2) * delta,
        // then the left tail, then the right tail — no term carries
        // count 5 or count 1, which is what the equal branch forbids)
        assert_eq!(expected, lookup(2) * 4.0 + lookup(11) * 2.0 + lookup(13) * 4.0);
        assert_eq!(
            differing_observed_mass_from_units(
                &left,
                &left_observed,
                &right,
                &right_observed,
            )
            .to_bits(),
            expected.to_bits()
        );

        // the sweep: every branch fires, on both key lists, at every
        // shape; the array walk's bits are the probe walk's bits
        let mut seed = 0x0b5e_c0deu64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        let mut equal_with_delta = 0usize;
        let mut left_only = 0usize;
        let mut right_only = 0usize;
        for _ in 0..2000 {
            let span = 1 + next() % 12;
            let mut left: Vec<(u64, u32)> = Vec::new();
            let mut right: Vec<(u64, u32)> = Vec::new();
            for key in 0..span as u64 {
                if next() % 4 != 0 {
                    left.push((key * 3, 1 + next() % 4));
                }
                if next() % 4 != 0 {
                    right.push((key * 3, 1 + next() % 4));
                }
            }
            // (the branch census: the equal-branch delta product and both
            // tails must actually occur, so the gate cannot pass vacuously)
            let (mut i, mut j) = (0usize, 0usize);
            while i < left.len() && j < right.len() {
                match left[i].0.cmp(&right[j].0) {
                    std::cmp::Ordering::Equal => {
                        if left[i].1 != right[j].1 {
                            equal_with_delta += 1;
                        }
                        i += 1;
                        j += 1;
                    }
                    std::cmp::Ordering::Less => {
                        left_only += 1;
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        right_only += 1;
                        j += 1;
                    }
                }
            }
            left_only += left.len() - i;
            right_only += right.len() - j;
            let left_observed: Vec<f64> =
                left.iter().map(|&(key, _)| lookup(key)).collect();
            let right_observed: Vec<f64> =
                right.iter().map(|&(key, _)| lookup(key)).collect();
            let expected = differing_observed_mass(&left, &right, &lookup);
            let got = differing_observed_mass_from_units(
                &left,
                &left_observed,
                &right,
                &right_observed,
            );
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "left {left:?} right {right:?}"
            );
        }
        assert!(equal_with_delta > 100, "equal-branch deltas {equal_with_delta}");
        assert!(left_only > 100, "left-only keys {left_only}");
        assert!(right_only > 100, "right-only keys {right_only}");
    }

    #[test]
    fn the_degenerate_tie_certificate_matches_the_naive_band_oracle() {
        // THE WALL FIX'S EXACTNESS GATE (the chrIX oracle discipline):
        // the degenerate-locus cut-index path against a VERBATIM copy
        // of the pre-fix naive band matrix. Randomized folds and
        // observed maps; the naive path materializes every band and
        // derives its <=cut predicate set; the tiered path must emit
        // the IDENTICAL predicate set and drive cluster_form_qual to
        // the IDENTICAL state (k, excluded, alternative, cluster_size,
        // cut, shape, p) — the bounds are exact-deciding, not
        // approximating, and the fallback computes the very distance
        // the band would have carried.
        let mut seed = 0x5eed_1234u64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for case in 0..3000 {
            let n_folds = 2 + (next() % 9) as usize; // 2..=10 folds
            let pair_count = n_folds * (n_folds + 1) / 2;
            let folds: Vec<Fold> = (0..n_folds)
                .map(|_| {
                    let mut nodes: Vec<(u64, u32)> = (0..(next() % 7))
                        .map(|_| (1 + (next() % 24) as u64, 1 + next() % 3))
                        .collect();
                    nodes.sort_unstable();
                    nodes.dedup_by(|a, b| {
                        if a.0 == b.0 {
                            b.1 += a.1;
                            true
                        } else {
                            false
                        }
                    });
                    let mut edges: Vec<(u64, u32)> = (0..(next() % 5))
                        .map(|_| (1000 + (next() % 24) as u64, 1 + next() % 2))
                        .collect();
                    edges.sort_unstable();
                    edges.dedup_by(|a, b| {
                        if a.0 == b.0 {
                            b.1 += a.1;
                            true
                        } else {
                            false
                        }
                    });
                    Fold {
                        seq: vec![],
                        len: 0,
                        walk: vec![],
                        node_positions: HashMap::new(),
                        members: vec![],
                        nodes,
                        edges,
                    }
                })
                .collect();
            let observed_node: HashMap<u64, f64> = (0..8)
                .map(|_| (1 + (next() % 30) as u64, (next() % 97) as f64))
                .collect();
            let observed_edge: HashMap<u64, f64> = (0..6)
                .map(|_| (1000 + (next() % 30) as u64, (next() % 89) as f64))
                .collect();
            let node_fn = |key: u64| observed_node.get(&key).copied().unwrap_or(0.0);
            let edge_fn = |key: u64| observed_edge.get(&key).copied().unwrap_or(0.0);
            // The degenerate locus: every class ties at the pure-E
            // score, the winner is (0, 0) by flat order, and the
            // spectrum is every other pair — the instrument's own
            // construction, all scores 1.0.
            let called: Vec<(usize, usize)> = (0..pair_count).map(unrank_pair).collect();
            let (winner_first, winner_second) = called[0];
            let winner_nodes =
                merged_multiset(&folds[winner_first].nodes, &folds[winner_second].nodes);
            let winner_edges =
                merged_multiset(&folds[winner_first].edges, &folds[winner_second].edges);
            let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
                .filter_map(|flat| {
                    let (first, second) = unrank_pair(flat);
                    if first == winner_first && second == winner_second {
                        return None;
                    }
                    let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
                    let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
                    let distance =
                        differing_observed_mass(&nodes, &winner_nodes, &node_fn)
                            + differing_observed_mass(&edges, &winner_edges, &edge_fn);
                    let signature = signature_cosine_distance(
                        &nodes, &edges, &winner_nodes, &winner_edges,
                    );
                    Some((distance, 1.0, signature, first, second, flat))
                })
                .collect();
            spectrum
                .sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
            let mut sorted: Vec<f64> =
                spectrum.iter().map(|&(distance, ..)| distance).collect();
            sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
            let cut = spectrum_knee(&sorted).cut;
            // The VERBATIM naive oracle (the pre-fix band
            // construction, kept here as the test's oracle).
            let mut oracle_bands: Vec<(usize, Vec<f64>)> = Vec::new();
            for &(first, second) in called.iter().skip(1) {
                let index = spectrum
                    .iter()
                    .position(|&(_, _, _, a, b, _)| a == first && b == second)
                    .expect("called class missing from the winner's spectrum");
                let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
                let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
                let band: Vec<f64> = spectrum
                    .iter()
                    .map(|&(_, _, _, other_first, other_second, _)| {
                        let other_nodes =
                            merged_multiset(&folds[other_first].nodes, &folds[other_second].nodes);
                        let other_edges =
                            merged_multiset(&folds[other_first].edges, &folds[other_second].edges);
                        differing_observed_mass(&nodes, &other_nodes, &node_fn)
                            + differing_observed_mass(&edges, &other_edges, &edge_fn)
                    })
                    .collect();
                oracle_bands.push((index, band));
            }
            let oracle_cut_lists: Vec<(usize, Vec<usize>)> = oracle_bands
                .iter()
                .map(|(index, band)| {
                    (
                        *index,
                        band.iter()
                            .enumerate()
                            .filter(|&(_, &distance)| distance <= cut)
                            .map(|(position, _)| position)
                            .collect(),
                    )
                })
                .collect();
            let tiered = tied_class_cut_lists(
                &folds, &called, &spectrum, &node_fn, &edge_fn, cut,
            );
            // The oracle stays fully independent (it derives each band's
            // predicate set from the naive band vectors); the packed form
            // is unpacked back into the oracle's own representation, so
            // this also proves the bit order and the ASCENDING iteration.
            let tiered_indices: Vec<(usize, Vec<usize>)> = tiered
                .iter()
                .map(|(index, set)| (*index, set.iter().collect()))
                .collect();
            assert_eq!(
                tiered_indices, oracle_cut_lists,
                "case {case}: predicate sets differ"
            );
            let pairs: Vec<(f64, f64)> = spectrum
                .iter()
                .map(|&(distance, score, ..)| (distance, score))
                .collect();
            let state_bands =
                cluster_form_qual(1.0, &pairs, &TieCertificate::Bands(&oracle_bands));
            let state_lists =
                cluster_form_qual(1.0, &pairs, &TieCertificate::CutIndices(&tiered));
            // THE FAST PATH (`tie_cluster_fast`: k as connected components,
            // excluded from the shrinking candidate set) must answer the BAND
            // ORACLE identically - its safety net, at both tie shapes.
            let fast = tie_cluster_fast(&folds, &called, &spectrum, &node_fn, &edge_fn);
            assert_eq!(fast.k, state_bands.k, "case {case}: FAST k differs");
            assert_eq!(
                fast.excluded, state_bands.excluded,
                "case {case}: FAST excluded differs"
            );
            assert_eq!(
                fast.alternative, state_bands.alternative,
                "case {case}: FAST alternative differs"
            );
            assert_eq!(
                fast.cluster_size, state_bands.cluster_size,
                "case {case}: FAST cluster_size differs"
            );
            assert_eq!(fast.cut, state_bands.cut, "case {case}: FAST cut differs");
            assert_eq!(state_lists.k, state_bands.k, "case {case}: k differs");
            assert_eq!(state_lists.excluded, state_bands.excluded, "case {case}: excluded differs");
            assert_eq!(state_lists.alternative, state_bands.alternative, "case {case}: alternative differs");
            assert_eq!(state_lists.cluster_size, state_bands.cluster_size, "case {case}: cluster_size differs");
            assert_eq!(state_lists.cut, state_bands.cut, "case {case}: cut differs");
            assert_eq!(state_lists.shape, state_bands.shape, "case {case}: shape differs");
            assert_eq!(state_lists.p, state_bands.p, "case {case}: p differs");
        }
    }

    #[test]
    fn the_degenerate_certificate_handles_the_monster_shape_fast() {
        // The monster-locus shape at test scale (the chrXII L86
        // class): MANY folds, EVERY class tied at the pure-E score —
        // the naive certificate is (pair_count - 1) squared distances
        // (at the real locus 35,510 x 35,510 = 1.26e9, a 4,802s grind
        // and a ~6GB emission that breached the 64GiB guard); the
        // tiered path must reach the IDENTICAL cluster state well
        // inside the liveness bound, with the oracle proving the
        // equality at a size the oracle can afford (40 folds -> 820
        // classes, 819 x 819 naive distances).
        let started = Instant::now();
        let mut seed = 0x1234_abcdu64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        let n_folds = 40usize;
        let pair_count = n_folds * (n_folds + 1) / 2;
        let folds: Vec<Fold> = (0..n_folds)
            .map(|fold| {
                // Repeat-domain-shaped multisets: every fold carries a
                // few HEAVY observed nodes at varying copy counts plus
                // disjoint unobserved material.
                let mut nodes: Vec<(u64, u32)> = vec![(1, 1 + next() % 9), (2, 1 + next() % 9)];
                nodes.extend((0..(next() % 6)).map(|_| {
                    (50 + (next() % 200) as u64, 1 + next() % 3)
                }));
                nodes.sort_unstable();
                nodes.dedup_by(|a, b| {
                    if a.0 == b.0 {
                        b.1 += a.1;
                        true
                    } else {
                        false
                    }
                });
                let mut edges: Vec<(u64, u32)> = vec![(1001, 1 + next() % 5)];
                edges.extend((0..(next() % 4)).map(|_| {
                    (1100 + (next() % 100) as u64, 1 + next() % 2)
                }));
                edges.sort_unstable();
                edges.dedup_by(|a, b| {
                    if a.0 == b.0 {
                        b.1 += a.1;
                        true
                    } else {
                        false
                    }
                });
                let _ = fold;
                Fold {
                    seq: vec![],
                    len: 0,
                    walk: vec![],
                    node_positions: HashMap::new(),
                    members: vec![],
                    nodes,
                    edges,
                }
            })
            .collect();
        let observed_node: HashMap<u64, f64> =
            [(1u64, 37.0), (2u64, 11.0), (3u64, 53.0)].into_iter().collect();
        let observed_edge: HashMap<u64, f64> =
            [(1001u64, 29.0), (1002u64, 7.0)].into_iter().collect();
        let node_fn = |key: u64| observed_node.get(&key).copied().unwrap_or(0.0);
        let edge_fn = |key: u64| observed_edge.get(&key).copied().unwrap_or(0.0);
        let called: Vec<(usize, usize)> = (0..pair_count).map(unrank_pair).collect();
        let (winner_first, winner_second) = called[0];
        let winner_nodes =
            merged_multiset(&folds[winner_first].nodes, &folds[winner_second].nodes);
        let winner_edges =
            merged_multiset(&folds[winner_first].edges, &folds[winner_second].edges);
        let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
            .filter_map(|flat| {
                let (first, second) = unrank_pair(flat);
                if first == winner_first && second == winner_second {
                    return None;
                }
                let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
                let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
                let distance =
                    differing_observed_mass(&nodes, &winner_nodes, &node_fn)
                        + differing_observed_mass(&edges, &winner_edges, &edge_fn);
                let signature =
                    signature_cosine_distance(&nodes, &edges, &winner_nodes, &winner_edges);
                Some((distance, 1.0, signature, first, second, flat))
            })
            .collect();
        spectrum.sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
        let mut sorted: Vec<f64> = spectrum.iter().map(|&(distance, ..)| distance).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
        let cut = spectrum_knee(&sorted).cut;
        let mut oracle_bands: Vec<(usize, Vec<f64>)> = Vec::new();
        for &(first, second) in called.iter().skip(1) {
            let index = spectrum
                .iter()
                .position(|&(_, _, _, a, b, _)| a == first && b == second)
                .expect("called class missing from the winner's spectrum");
            let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
            let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
            let band: Vec<f64> = spectrum
                .iter()
                .map(|&(_, _, _, other_first, other_second, _)| {
                    let other_nodes =
                        merged_multiset(&folds[other_first].nodes, &folds[other_second].nodes);
                    let other_edges =
                        merged_multiset(&folds[other_first].edges, &folds[other_second].edges);
                    differing_observed_mass(&nodes, &other_nodes, &node_fn)
                        + differing_observed_mass(&edges, &other_edges, &edge_fn)
                })
                .collect();
            oracle_bands.push((index, band));
        }
        let oracle_cut_lists: Vec<(usize, Vec<usize>)> = oracle_bands
            .iter()
            .map(|(index, band)| {
                (
                    *index,
                    band.iter()
                        .enumerate()
                        .filter(|&(_, &distance)| distance <= cut)
                        .map(|(position, _)| position)
                        .collect(),
                )
            })
            .collect();
        let tiered = tied_class_cut_lists(&folds, &called, &spectrum, &node_fn, &edge_fn, cut);
        let tiered_indices: Vec<(usize, Vec<usize>)> = tiered
            .iter()
            .map(|(index, set)| (*index, set.iter().collect()))
            .collect();
        assert_eq!(
            tiered_indices, oracle_cut_lists,
            "monster shape: predicate sets differ"
        );
        let pairs: Vec<(f64, f64)> = spectrum
            .iter()
            .map(|&(distance, score, ..)| (distance, score))
            .collect();
        let state_bands =
            cluster_form_qual(1.0, &pairs, &TieCertificate::Bands(&oracle_bands));
        let state_lists =
            cluster_form_qual(1.0, &pairs, &TieCertificate::CutIndices(&tiered));
        assert_eq!(state_lists.k, state_bands.k);
        assert_eq!(state_lists.excluded, state_bands.excluded);
        assert_eq!(state_lists.p, state_bands.p);
        assert!(started.elapsed().as_secs() < 60, "the liveness bound");
    }
    #[test]
    fn the_non_degenerate_tie_certificate_matches_the_naive_band_oracle() {
        // THE PER-PAIR CONSTANT'S OWN EXACTNESS GATE (the fold-merge
        // open term's fix): the NON-degenerate tie certificate — a
        // STRICT tie subset (called.len() > 1 but called.len() <
        // pair_count, the fold-heavy shapes' uncollapsed tie sets;
        // the pre-fix naive path materialized |called| - 1 bands x
        // |spectrum| distances there) — against a VERBATIM copy of
        // the pre-fix naive band matrix. Randomized folds and observed
        // maps; the naive path materializes every band and derives
        // its <=cut predicate set; the cut-index path (per-fold mass
        // signatures computed once, O(1) two-sided bounds per entry,
        // the exact walk as the undecided-band fallback) must emit the
        // IDENTICAL predicate set and drive cluster_form_qual to the
        // IDENTICAL state (k, excluded, alternative, cluster_size,
        // cut, shape, p) — the bounds are exact-deciding, never
        // approximating, and the fallback computes the very distance
        // the band would have carried.
        let mut seed = 0x5eed_5678u64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for case in 0..4000 {
            let n_folds = 3 + (next() % 10) as usize; // 3..=12 folds
            let pair_count = n_folds * (n_folds + 1) / 2;
            let folds: Vec<Fold> = (0..n_folds)
                .map(|_| {
                    let mut nodes: Vec<(u64, u32)> = (0..(next() % 8))
                        .map(|_| (1 + (next() % 30) as u64, 1 + next() % 3))
                        .collect();
                    nodes.sort_unstable();
                    nodes.dedup_by(|a, b| {
                        if a.0 == b.0 {
                            b.1 += a.1;
                            true
                        } else {
                            false
                        }
                    });
                    let mut edges: Vec<(u64, u32)> = (0..(next() % 6))
                        .map(|_| (1000 + (next() % 30) as u64, 1 + next() % 2))
                        .collect();
                    edges.sort_unstable();
                    edges.dedup_by(|a, b| {
                        if a.0 == b.0 {
                            b.1 += a.1;
                            true
                        } else {
                            false
                        }
                    });
                    Fold {
                        seq: vec![],
                        len: 0,
                        walk: vec![],
                        node_positions: HashMap::new(),
                        members: vec![],
                        nodes,
                        edges,
                    }
                })
                .collect();
            let observed_node: HashMap<u64, f64> = (0..9)
                .map(|_| (1 + (next() % 34) as u64, (next() % 97) as f64))
                .collect();
            let observed_edge: HashMap<u64, f64> = (0..7)
                .map(|_| (1000 + (next() % 34) as u64, (next() % 89) as f64))
                .collect();
            let node_fn = |key: u64| observed_node.get(&key).copied().unwrap_or(0.0);
            let edge_fn = |key: u64| observed_edge.get(&key).copied().unwrap_or(0.0);
            // THE NON-DEGENERATE SHAPE: the winner is the flat-order
            // first pair (the exhaustive scan's own convention) and
            // the tie set is a STRICT random subset of the remaining
            // pairs (2..=pair_count - 1 members — always a proper
            // subset, never the full degenerate class).
            let winner = unrank_pair(0);
            let remaining: Vec<(usize, usize)> = (1..pair_count).map(unrank_pair).collect();
            let tie_count = 1 + (next() as usize % (remaining.len() - 1));
            let mut called: Vec<(usize, usize)> = vec![winner];
            let mut pool = remaining.clone();
            for _ in 0..tie_count {
                let pick = (next() as usize) % pool.len();
                called.push(pool.swap_remove(pick));
            }
            assert!(called.len() > 1 && called.len() < pair_count);
            let (winner_first, winner_second) = called[0];
            let winner_nodes =
                merged_multiset(&folds[winner_first].nodes, &folds[winner_second].nodes);
            let winner_edges =
                merged_multiset(&folds[winner_first].edges, &folds[winner_second].edges);
            let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
                .filter_map(|flat| {
                    let (first, second) = unrank_pair(flat);
                    if first == winner_first && second == winner_second {
                        return None;
                    }
                    let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
                    let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
                    let distance =
                        differing_observed_mass(&nodes, &winner_nodes, &node_fn)
                            + differing_observed_mass(&edges, &winner_edges, &edge_fn);
                    let signature = signature_cosine_distance(
                        &nodes, &edges, &winner_nodes, &winner_edges,
                    );
                    Some((distance, 1.0, signature, first, second, flat))
                })
                .collect();
            spectrum
                .sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
            let mut sorted: Vec<f64> =
                spectrum.iter().map(|&(distance, ..)| distance).collect();
            sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
            let cut = spectrum_knee(&sorted).cut;
            // The VERBATIM naive oracle (the pre-fix NON-degenerate
            // band construction, kept here as the test's oracle).
            let mut oracle_bands: Vec<(usize, Vec<f64>)> = Vec::new();
            for &(first, second) in called.iter().skip(1) {
                let index = spectrum
                    .iter()
                    .position(|&(_, _, _, a, b, _)| a == first && b == second)
                    .expect("called class missing from the winner's spectrum");
                let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
                let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
                let band: Vec<f64> = spectrum
                    .iter()
                    .map(|&(_, _, _, other_first, other_second, _)| {
                        let other_nodes =
                            merged_multiset(&folds[other_first].nodes, &folds[other_second].nodes);
                        let other_edges =
                            merged_multiset(&folds[other_first].edges, &folds[other_second].edges);
                        differing_observed_mass(&nodes, &other_nodes, &node_fn)
                            + differing_observed_mass(&edges, &other_edges, &edge_fn)
                    })
                    .collect();
                oracle_bands.push((index, band));
            }
            let oracle_cut_lists: Vec<(usize, Vec<usize>)> = oracle_bands
                .iter()
                .map(|(index, band)| {
                    (
                        *index,
                        band.iter()
                            .enumerate()
                            .filter(|&(_, &distance)| distance <= cut)
                            .map(|(position, _)| position)
                            .collect(),
                    )
                })
                .collect();
            let tiered = tied_class_cut_lists(
                &folds, &called, &spectrum, &node_fn, &edge_fn, cut,
            );
            // The oracle stays fully independent (it derives each band's
            // predicate set from the naive band vectors); the packed form
            // is unpacked back into the oracle's own representation, so
            // this also proves the bit order and the ASCENDING iteration.
            let tiered_indices: Vec<(usize, Vec<usize>)> = tiered
                .iter()
                .map(|(index, set)| (*index, set.iter().collect()))
                .collect();
            assert_eq!(
                tiered_indices, oracle_cut_lists,
                "case {case}: predicate sets differ"
            );
            let pairs: Vec<(f64, f64)> = spectrum
                .iter()
                .map(|&(distance, score, ..)| (distance, score))
                .collect();
            let state_bands =
                cluster_form_qual(1.0, &pairs, &TieCertificate::Bands(&oracle_bands));
            let state_lists =
                cluster_form_qual(1.0, &pairs, &TieCertificate::CutIndices(&tiered));
            // THE FAST PATH (`tie_cluster_fast`: k as connected components,
            // excluded from the shrinking candidate set) must answer the BAND
            // ORACLE identically - its safety net, at both tie shapes.
            let fast = tie_cluster_fast(&folds, &called, &spectrum, &node_fn, &edge_fn);
            assert_eq!(fast.k, state_bands.k, "case {case}: FAST k differs");
            assert_eq!(
                fast.excluded, state_bands.excluded,
                "case {case}: FAST excluded differs"
            );
            assert_eq!(
                fast.alternative, state_bands.alternative,
                "case {case}: FAST alternative differs"
            );
            assert_eq!(
                fast.cluster_size, state_bands.cluster_size,
                "case {case}: FAST cluster_size differs"
            );
            assert_eq!(fast.cut, state_bands.cut, "case {case}: FAST cut differs");
            assert_eq!(state_lists.k, state_bands.k, "case {case}: k differs");
            assert_eq!(state_lists.excluded, state_bands.excluded, "case {case}: excluded differs");
            assert_eq!(state_lists.alternative, state_bands.alternative, "case {case}: alternative differs");
            assert_eq!(state_lists.cluster_size, state_bands.cluster_size, "case {case}: cluster_size differs");
            assert_eq!(state_lists.cut, state_bands.cut, "case {case}: cut differs");
            assert_eq!(state_lists.shape, state_bands.shape, "case {case}: shape differs");
            assert_eq!(state_lists.p, state_bands.p, "case {case}: p differs");
        }
    }
    #[test]
    fn canonical_scheme_selects_per_window_form() {
        // The canonical-scheme selection (slice D's frame repair): per
        // position the frame that spells the k-mer's canonical form
        // forward decides; a position whose canonical frame did not
        // qualify carries NO step; the chosen frame's node is kept.
        let k = 4u64;
        // A 16 bp sequence: window at 0 is canonical-forward (ACGT),
        // the windows at 6 and 12 are rc-canonical (T leads and the
        // complement of the trailing base loses the comparison).
        let seq = b"ACGTGTTAACGTTTGAC".to_vec();
        let seq_lo = 100u64;
        assert!(crate::genome_inference::mem_records::window_is_canonical_forward(
            &seq[0..4]
        ));
        assert!(!crate::genome_inference::mem_records::window_is_canonical_forward(
            &seq[6..10]
        ));
        let mut forward = BTreeMap::new();
        forward.insert(seq_lo + 0, 11);
        forward.insert(seq_lo + 6, 13); // forward-qualified, rc-canonical
        let mut reverse = BTreeMap::new();
        reverse.insert(seq_lo + 6, -13); // the rc frame's node at the same position
        reverse.insert(seq_lo + 12, -17); // rc-only position, rc-canonical window
        assert!(!crate::genome_inference::mem_records::window_is_canonical_forward(
            &seq[12..16]
        ));
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        // The canonical-forward position keeps the forward node; the
        // rc-canonical positions take the reverse frame's node (or
        // drop when only the forward frame qualified); the forward-only
        // rc-canonical position is REPLACED by the reverse node here.
        assert_eq!(
            steps,
            vec![(seq_lo + 0, 11, 0u8), (seq_lo + 6, -13, 1u8), (seq_lo + 12, -17, 1u8)]
        );
    }

    #[test]
    fn canonical_scheme_drops_forward_only_rc_canonical_positions() {
        // A position where ONLY the forward frame qualified and the
        // window is rc-canonical carries NO step: no census record
        // (derived under the same scheme) can anchor there.
        let k = 4u64;
        let seq = b"ACGTGTTAACGTTTGAC".to_vec();
        let seq_lo = 0u64;
        let mut forward = BTreeMap::new();
        forward.insert(6, 13);
        let reverse: BTreeMap<u64, i32> = BTreeMap::new();
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        assert!(steps.is_empty());
        // With the rc frame qualifying there too, the position keeps
        // the rc frame's node.
        let mut reverse = BTreeMap::new();
        reverse.insert(6, -13);
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        assert_eq!(steps, vec![(6, -13, 1u8)]);
    }

    /// The pre-rebuild per-shift verification (phase 0's code, the
    /// reference for the hoisted union-range verify): one fetch per
    /// (occurrence x shift), with the fetch closure carrying the
    /// fetch_seq/Sources clamping semantics (high bound clamped to the
    /// path length; empty crop when the low bound reaches it).
    fn old_verify_key_at(
        shape: &KeyShape,
        occ: &CensusOccurrence,
        reads: &[Vec<u8>],
        fetch: &dyn Fn(u64, u64) -> Vec<u8>,
        k: u64,
        shift: i64,
    ) -> bool {
        let lo = (occ.start as i64 + shift - 1).max(0) as u64;
        let hi = occ.start + shape.span + 1;
        let seq = fetch(lo, hi);
        let origin = occ.start as i64 + shift;
        let read = &reads[shape.read];
        let forward = physical_forward(occ.orientation, shape.mirrored);
        for (j, &(node, canonical_pos)) in shape.canonical.iter().enumerate() {
            let _ = node;
            let window_lo =
                anchor_window_start(canonical_pos, origin as u64, shape.span, k, occ.orientation);
            let rel_lo = window_lo as i64 - lo as i64;
            if rel_lo < 0 || rel_lo + k as i64 > seq.len() as i64 {
                return false;
            }
            let window = &seq[rel_lo as usize..(rel_lo + k as i64) as usize];
            let own_pos = shape.own_positions[j] as usize;
            let kmer = &read[own_pos..own_pos + k as usize];
            let ok = if forward {
                window == kmer
            } else {
                revcomp(window) == kmer
            };
            if !ok {
                return false;
            }
        }
        true
    }

    #[test]
    fn hoisted_verify_reproduces_the_per_shift_fetch_verification() {
        // Phase 1, lever 1: the once-per-(occurrence, candidate) union
        // range fetch must reproduce, for every origin shift, exactly
        // the verification the three per-shift fetches performed —
        // including the edge cases where the clamped crop is short or
        // empty (occurrence start near 0; range overhanging the path
        // end), both orientations, both mirror states, and verdicts of
        // both true and false.
        let k = 4u64;
        let path: Vec<u8> = b"ACGTACGTACGTTTGACACGTGCATATCGGATTCGCA".to_vec();
        let path_len = path.len() as u64;
        let fetch = |lo: u64, hi: u64| -> Vec<u8> {
            let hi = hi.min(path_len);
            if lo >= hi {
                Vec::new()
            } else {
                path[lo as usize..hi as usize].to_vec()
            }
        };
        // anchors at canonical positions 0 and 8, own k-mers at read
        // positions 2 and 10; span 8 - 0 + 4 = 12.
        let make_shape = |mirrored: bool| KeyShape {
            canonical: vec![(7, 0), (9, 8)],
            span: 12,
            read: 0,
            mirrored,
            own_positions: vec![2, 10],
        };
        let mut read = path[2..16].to_vec();
        let reads_match = vec![read.clone()];
        read[3] = match read[3] {
            b'A' => b'T',
            other => b'A' + (other == b'A') as u8,
        };
        let reads_mismatch = vec![read];
        // Interior starts verify true (orientation 0, unmirrored);
        // near-zero starts exercise the low-bound clamp, path-end
        // starts exercise the high-bound clamp and the empty crop.
        for start in [0u64, 1, 2, 3, 5, 20, 30, 33, 34, 36, 38, 39, 41] {
            for orientation in [0u8, 1] {
                for mirrored in [false, true] {
                    let occ = CensusOccurrence {
                        path: 0,
                        start,
                        orientation,
                        partitions: vec![],
                        intervals: vec![],
                    };
                    let shape = make_shape(mirrored);
                    for reads in [&reads_match, &reads_mismatch] {
                        // The union range fetch (the hoist): one fetch
                        // covering all three shifts' ranges.
                        let union_lo = (occ.start as i64 - 2).max(0) as u64;
                        let union_hi = occ.start + shape.span + 1;
                        let seq = fetch(union_lo, union_hi);
                        for shift in [-1i64, 0, 1] {
                            assert_eq!(
                                verify_key_at(&shape, &occ, reads, &seq, union_lo, k, shift),
                                old_verify_key_at(&shape, &occ, reads, &fetch, k, shift),
                                "hoisted verify disagrees with the per-shift fetch \
                                 (start {start}, orientation {orientation}, mirrored \
                                 {mirrored}, shift {shift})",
                            );
                        }
                    }
                }
            }
        }
        // And a true verdict exists (the equivalence is not vacuous):
        // the read spells path[2..16], its first anchor's own k-mer at
        // read position 2, so the occurrence start is 4.
        let occ = CensusOccurrence {
            path: 0,
            start: 4,
            orientation: 0,
            partitions: vec![],
            intervals: vec![],
        };
        let shape = make_shape(false);
        let union_lo = 0u64;
        let seq = fetch(union_lo, occ.start + shape.span + 1);
        assert!(verify_key_at(&shape, &occ, &reads_match, &seq, union_lo, k, 0));
        assert!(!verify_key_at(&shape, &occ, &reads_mismatch, &seq, union_lo, k, 0));
    }

    /// THE COMMITTED ENUMERATION (the pre-DP anchor correspondence,
    /// copied verbatim for the equivalence proof — it descends every
    /// monotone assignment; exponential on repeat-dense rows, which
    /// is exactly why the DP replaces it at chrIV scale).
    fn enumerate_correspondence(
        canonical: &[(i32, u64)],
        orientation: u8,
        row_positions: &dyn Fn(u32) -> Vec<(u64, i32)>,
    ) -> (Vec<Option<u64>>, Vec<usize>) {
        let m = canonical.len();
        let mut candidates: Vec<Vec<u64>> = Vec::with_capacity(m);
        for &(node, _) in canonical {
            let required = required_sign(node, orientation);
            let positions = row_positions(required.unsigned_abs());
            candidates.push(
                positions
                    .iter()
                    .filter(|&&(_, sign)| sign == required.signum())
                    .map(|&(bp, _)| bp)
                    .collect(),
            );
        }
        let mut pinned: Vec<Option<u64>> = vec![None; m];
        let mut pinned_seen: Vec<bool> = vec![false; m];
        let mut assignments = 0u64;
        let mut chain: Vec<(usize, u64)> = Vec::with_capacity(m);
        fn enumerate(
            index: usize,
            candidates: &[Vec<u64>],
            increasing: bool,
            chain: &mut Vec<(usize, u64)>,
            pinned: &mut Vec<Option<u64>>,
            pinned_seen: &mut Vec<bool>,
            assignments: &mut u64,
        ) {
            if index == candidates.len() {
                *assignments += 1;
                for &(j, bp) in chain.iter() {
                    if !pinned_seen[j] {
                        pinned[j] = Some(bp);
                        pinned_seen[j] = true;
                    } else if pinned[j] != Some(bp) {
                        pinned[j] = None;
                    }
                }
                return;
            }
            if candidates[index].is_empty() {
                enumerate(index + 1, candidates, increasing, chain, pinned, pinned_seen, assignments);
                return;
            }
            for &bp in &candidates[index] {
                let ok = match chain.last() {
                    Some(&(_, last)) => {
                        if increasing {
                            bp > last
                        } else {
                            bp < last
                        }
                    }
                    None => true,
                };
                if ok {
                    chain.push((index, bp));
                    enumerate(index + 1, candidates, increasing, chain, pinned, pinned_seen, assignments);
                    chain.pop();
                }
            }
        }
        enumerate(
            0,
            &candidates,
            orientation == 0,
            &mut chain,
            &mut pinned,
            &mut pinned_seen,
            &mut assignments,
        );
        let ambiguous: Vec<usize> = (0..m)
            .filter(|&j| !candidates[j].is_empty() && pinned[j].is_none())
            .collect();
        let _ = assignments;
        (pinned, ambiguous)
    }

    #[test]
    fn anchor_correspondence_dp_equals_the_committed_enumeration() {
        // The exact-DP equivalence: over randomized synthetic walks
        // (a fixed-seed LCG for determinism), covering unique pins,
        // repeat-node ambiguity, absent anchors, zero-chain cases and
        // both orientations, the DP's (pinned, ambiguous) must equal
        // the committed enumeration's on EVERY case.
        let mut seed: u64 = 0x243F6A8885A308D3;
        let mut rng = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..3000u32 {
            let m = (rng() % 6) as usize + 1; // 1..=6 anchors
            let mut canonical: Vec<(i32, u64)> = Vec::with_capacity(m);
            let mut pos = 0u64;
            for j in 0..m {
                pos += rng() % 7;
                let node = 1 + (rng() % 4) as i32; // nodes 1..=4 (collisions likely)
                let sign = if rng() % 2 == 0 { 1 } else { -1 };
                canonical.push((sign * node, pos));
                let _ = j;
            }
            // the row's node positions: each node at 0..=8 positions
            let mut row: std::collections::HashMap<u32, Vec<(u64, i32)>> =
                std::collections::HashMap::new();
            for node in 1u32..=4 {
                let count = (rng() % 9) as usize; // 0..=8 positions
                let mut positions: Vec<(u64, i32)> = Vec::with_capacity(count);
                let mut p = 0u64;
                for _ in 0..count {
                    p += rng() % 11;
                    let sign = if rng() % 2 == 0 { 1 } else { -1 };
                    positions.push((p, sign));
                }
                row.insert(node, positions);
            }
            let positions_of = move |node: u32| row.get(&node).cloned().unwrap_or_default();
            for orientation in [0u8, 1u8] {
                let expected = enumerate_correspondence(&canonical, orientation, &positions_of);
                let actual = anchor_correspondence(&canonical, orientation, &positions_of);
                assert_eq!(
                    expected, actual,
                    "DP differs from the committed enumeration at case {case} \
                     orientation {orientation}: canonical {canonical:?}"
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // The alignment-connected fold relaxation (stage 1): the merge
    // pass and the truth-second correspondence, unit-tested against
    // toy fold spaces that name every rule clause.
    // ---------------------------------------------------------------

    /// A toy fold: nodes/edges sorted+deduped multisets, members by
    /// name, an all-identical column by default.
    fn alignfold_toy(
        nodes: &[(u64, u32)],
        edges: &[(u64, u32)],
        members: &[&str],
    ) -> Fold {
        let mut sorted: Vec<(u64, u32)> = nodes.to_vec();
        sorted.sort_unstable();
        sorted.dedup_by(|a, b| {
            if a.0 == b.0 {
                b.1 += a.1;
                true
            } else {
                false
            }
        });
        let mut esorted: Vec<(u64, u32)> = edges.to_vec();
        esorted.sort_unstable();
        esorted.dedup_by(|a, b| {
            if a.0 == b.0 {
                b.1 += a.1;
                true
            } else {
                false
            }
        });
        Fold {
            seq: vec![b'A', b'C', b'G', b'T'],
            len: 4,
            walk: vec![(0, 1)],
            node_positions: HashMap::new(),
            members: members
                .iter()
                .map(|&name| MemberRow {
                    path_name: name.to_string(),
                    start: 0,
                    end: 4,
                })
                .collect(),
            nodes: sorted,
            edges: esorted,
        }
    }

    #[test]
    fn score_matrix_row_compression_is_value_identical() {
        // Candidate C's exactness proof, at the accessor level: the
        // compressed rows return the dense form's own bits at every
        // probe, the pinned enumeration is the dense filter's own
        // order, the class scan's accumulation is bit-for-bit the
        // committed dense loop's, the merge's predicates agree with
        // the dense bit comparisons, and the streamed JSON is the
        // committed serde bytes — for BOTH row encodings and across
        // encodings (a dense row and a sparse row holding the same
        // values must hash and compare equal).
        let elsewhere = -14.418156f64;
        // (folds x units, the pinned fraction) cases: a no-pin-heavy
        // row (takes Sparse), a pinned-heavy row (takes Dense), and
        // an all-no-pin row.
        let shapes = [(6usize, 40usize, 0.1f64), (6, 40, 0.9), (6, 40, 0.0)];
        for (n_folds, n_units, pinned_fraction) in shapes {
            let mut dense: Vec<Vec<f64>> = Vec::new();
            let mut pinned: Vec<Vec<bool>> = Vec::new();
            let mut mismatches: Vec<Vec<u32>> = Vec::new();
            for fold in 0..n_folds {
                let mut row = vec![elsewhere; n_units];
                let mut pin = vec![false; n_units];
                let mut mism = vec![0u32; n_units];
                // a deterministic pattern: every k-th unit pins
                let stride = ((pinned_fraction * n_units as f64).round() as usize)
                    .max(1);
                for unit in (fold..n_units).step_by(stride.max(1)) {
                    pin[unit] = true;
                    row[unit] = elsewhere + (unit as f64) / (fold as f64 + 1.0);
                    mism[unit] = (unit % 7) as u32;
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(mism);
            }
            let matrix =
                ScoreMatrix::from_dense(dense.clone(), pinned.clone(), mismatches.clone(), elsewhere);
            // every probe equals the dense form's own value
            for fold in 0..n_folds {
                for unit in 0..n_units {
                    assert_eq!(matrix.ll(fold, unit).to_bits(), dense[fold][unit].to_bits());
                    assert_eq!(matrix.pinned(fold, unit), pinned[fold][unit]);
                    assert_eq!(matrix.mismatches(fold, unit), mismatches[fold][unit]);
                }
                // the pinned enumeration is the dense filter's order
                let dense_order: Vec<(usize, u32)> = (0..n_units)
                    .filter(|&u| pinned[fold][u])
                    .map(|u| (u, mismatches[fold][u]))
                    .collect();
                assert_eq!(matrix.pinned_units(fold).collect::<Vec<_>>(), dense_order);
            }
            // the encoding choice: sparse iff 16B/pinned < 13B/unit
            for fold in 0..n_folds {
                let pinned_count = pinned[fold].iter().filter(|&&p| p).count();
                let expected_sparse = (16 * pinned_count as u64) < (13 * n_units as u64);
                let is_sparse = matches!(matrix.rows[fold], ScoreRow::Sparse(_));
                assert_eq!(
                    is_sparse, expected_sparse,
                    "the derived encoding choice, fold {fold}"
                );
            }
            // the class scan: bit-for-bit the committed dense loop
            let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
            for first in 0..n_folds {
                for second in first..n_folds {
                    let mut committed = 0.0f64;
                    for (u, &count) in counts.iter().enumerate() {
                        committed +=
                            count * mix_logsumexp(dense[first][u], dense[second][u]);
                    }
                    assert_eq!(
                        matrix.class_pair_sum(first, second, &counts).to_bits(),
                        committed.to_bits(),
                        "the class pair sum, ({first},{second})"
                    );
                }
            }
            // the merge predicates: equal columns bit-equal and
            // hash-equal (across encodings too)
            let mut clone = dense[0].clone();
            clone[3] = f64::from_bits(elsewhere.to_bits() + 1); // ONE ULP
            let matrix2 = ScoreMatrix::from_dense(
                vec![dense[0].clone(), clone, dense[0].clone()],
                vec![pinned[0].clone(), pinned[0].clone(), pinned[0].clone()],
                vec![mismatches[0].clone(); 3],
                elsewhere,
            );
            assert!(matrix2.ll_bits_equal(0, 2));
            assert!(!matrix2.ll_bits_equal(0, 1));
            assert_eq!(matrix2.ll_bits_hash(0), matrix2.ll_bits_hash(2));
            // the streamed JSON is the committed serde bytes
            let mut streamed = Vec::new();
            matrix.write_ll_matrix_json(&mut streamed).unwrap();
            assert_eq!(streamed, serde_json::to_vec(&dense).unwrap());
        }
    }

    /// The unit-fold fixture: 3 folds x 8 units carrying every class
    /// the plan must tell apart — a constant column at the elsewhere
    /// constant (units 0, 6, 7), a constant column pinned on EVERY fold
    /// with one repeated value (unit 4, c != E), a bit-identical twin
    /// pair (units 1 and 5), a ONE-ULP near-twin that must NOT merge
    /// (unit 2), and a singleton column (unit 3). Returns the dense
    /// rows; the pin flags are the plan-irrelevant "value != E" form.
    fn unit_fold_fixture(elsewhere: f64) -> (Vec<Vec<f64>>, Vec<Vec<bool>>, Vec<Vec<u32>>) {
        let n_folds = 3usize;
        let n_units = 8usize;
        let mut dense = vec![vec![elsewhere; n_units]; n_folds];
        dense[0][1] = -7.5;
        dense[0][5] = -7.5;
        dense[0][2] = f64::from_bits((-7.5f64).to_bits() + 1); // ONE ULP
        dense[2][3] = -9.25;
        for row in dense.iter_mut() {
            row[4] = -6.0;
        }
        let pinned: Vec<Vec<bool>> = dense
            .iter()
            .map(|row| row.iter().map(|&ll| ll != elsewhere).collect())
            .collect();
        let mismatches = vec![vec![0u32; n_units]; n_folds];
        (dense, pinned, mismatches)
    }

    #[test]
    fn unit_fold_plan_is_the_exact_column_equality_partition() {
        // The grouping is the EXACT f64-bit column partition — a
        // one-ulp near-twin stays its own class, a bit-identical twin
        // merges, and a column pinned on every fold with one repeated
        // value is CONSTANT (mix(c,c) = c exactly, so it joins the
        // scalar). Both row encodings must give the SAME plan.
        let elsewhere = -14.0f64;
        let (dense, pinned, mismatches) = unit_fold_fixture(elsewhere);
        let counts: Vec<f64> = (1..=dense[0].len()).map(|u| u as f64).collect();
        // the sparse-encoded form (the fixture's own pin density)
        let mut sparse_matrix =
            ScoreMatrix::from_dense(dense.clone(), pinned.clone(), mismatches.clone(), elsewhere);
        assert!(matches!(sparse_matrix.rows[0], ScoreRow::Sparse(_)));
        assert!(sparse_matrix.fold_units(&counts).is_some());
        let sparse_plan = sparse_matrix.unit_fold().expect("the plan builds");
        // the all-pinned DENSE-encoded form (the same ll values)
        let all_pinned = vec![vec![true; dense[0].len()]; dense.len()];
        let mut dense_matrix =
            ScoreMatrix::from_dense(dense.clone(), all_pinned, mismatches.clone(), elsewhere);
        assert!(matches!(dense_matrix.rows[0], ScoreRow::Dense(_)));
        assert!(dense_matrix.fold_units(&counts).is_some());
        let dense_plan = dense_matrix.unit_fold().expect("the plan builds");
        assert_eq!(sparse_plan.groups, dense_plan.groups);
        assert_eq!(
            sparse_plan.constant_scalar.map(f64::to_bits),
            dense_plan.constant_scalar.map(f64::to_bits)
        );
        // THE EXPECTED CLASSES, hand-derived: {0,6,7} constant at E
        // (mass 16), {4} constant at -6.0 (mass 5) — the scalar is
        // their sum in representative order; {1,5} merge to mass 8 at
        // unit 1; {2} and {3} are singletons.
        assert_eq!(dense_plan.units, 8);
        assert_eq!(dense_plan.distinct_columns, 5);
        assert_eq!(dense_plan.constant_columns, 2);
        assert_eq!(dense_plan.constant_units, 4);
        assert_eq!(dense_plan.merged_units(), 3);
        assert_eq!(dense_plan.merged_columns(), 3);
        assert_eq!(dense_plan.groups(), &[(1u32, 8.0), (2, 3.0), (3, 4.0)]);
        assert_eq!(dense_plan.pair_terms, 4);
        let scalar = dense_plan.constant_scalar.expect("constant columns");
        assert_eq!(scalar.to_bits(), (16.0 * elsewhere + 5.0 * -6.0).to_bits());
    }

    #[test]
    fn unit_fold_keeps_the_constant_mass_in_the_class_ll() {
        // The scalar is NOT dropped: with 20 constant columns carrying
        // most of the mass, the folded class LL stays within a few ulp
        // of the committed walk (a dropped scalar would be off by
        // ~|count_total * E|, i.e. thousands).
        let elsewhere = -14.901063413333642f64;
        let n_folds = 4usize;
        let n_units = 21usize;
        let mut dense = vec![vec![elsewhere; n_units]; n_folds];
        for fold in 0..n_folds {
            dense[fold][20] = elsewhere + 0.5 * (fold as f64 + 1.0);
            dense[fold][7] = elsewhere + 1e-3;
        }
        let pinned: Vec<Vec<bool>> = dense
            .iter()
            .map(|row| row.iter().map(|&ll| ll != elsewhere).collect())
            .collect();
        let mismatches = vec![vec![0u32; n_units]; n_folds];
        let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + (u % 3) as f64).collect();
        let committed =
            ScoreMatrix::from_dense(dense.clone(), pinned.clone(), mismatches.clone(), elsewhere);
        let mut folded =
            ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
        assert!(folded.fold_units(&counts).is_some());
        for first in 0..n_folds {
            for second in first..n_folds {
                let exact = committed.class_pair_sum(first, second, &counts);
                let reduced = folded.class_pair_sum(first, second, &counts);
                let relative = ((reduced - exact) / exact).abs();
                assert!(
                    relative < 1e-13,
                    "the constant mass stays in the value: ({first},{second}) \
                     exact {exact} reduced {reduced}"
                );
            }
        }
    }

    #[test]
    fn unit_fold_is_bit_identical_when_nothing_folds() {
        // The identity case: no constant column and no bit-identical
        // column pair, so the surviving groups ARE the units in order
        // — the folded walk executes the committed op sequence and
        // must reproduce its bits exactly.
        let elsewhere = -13.75f64;
        let n_folds = 2usize;
        let n_units = 5usize;
        let mut dense = vec![vec![elsewhere; n_units]; n_folds];
        for fold in 0..n_folds {
            for unit in 0..n_units {
                dense[fold][unit] = elsewhere + 3.0 * unit as f64 + fold as f64 * 1.5;
            }
        }
        let pinned = vec![vec![true; n_units]; n_folds];
        let mismatches = vec![vec![0u32; n_units]; n_folds];
        let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
        let committed =
            ScoreMatrix::from_dense(dense.clone(), pinned.clone(), mismatches.clone(), elsewhere);
        let mut folded = ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
        let plan = folded.fold_units(&counts).expect("the plan builds");
        assert_eq!(plan.pair_terms, n_units);
        assert!(plan.constant_scalar.is_none());
        for first in 0..n_folds {
            for second in first..n_folds {
                assert_eq!(
                    folded.class_pair_sum(first, second, &counts).to_bits(),
                    committed.class_pair_sum(first, second, &counts).to_bits(),
                    "no unit folds: the walk is the committed one, ({first},{second})"
                );
            }
        }
    }

    #[test]
    fn unit_fold_declines_with_one_fold() {
        // One fold leaves every column trivially constant, so the plan
        // declines and the committed walk stays — measured nowhere,
        // never a silent rewrite.
        let elsewhere = -12.5f64;
        let dense = vec![vec![elsewhere, elsewhere + 1.0, elsewhere + 2.0]];
        let pinned = vec![vec![true, true, true]];
        let mismatches = vec![vec![0u32; 3]];
        let counts = vec![1.0f64, 2.0, 3.0];
        let committed =
            ScoreMatrix::from_dense(dense.clone(), pinned.clone(), mismatches.clone(), elsewhere);
        let mut declined = ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
        assert!(declined.fold_units(&counts).is_none());
        assert_eq!(
            declined.class_pair_sum(0, 0, &counts).to_bits(),
            committed.class_pair_sum(0, 0, &counts).to_bits()
        );
    }

    #[test]
    fn class_scan_sparse_skip_union_walk() {
        // THE OWNER'S RULING 3's own exactness proof, at the scan
        // level: the sparse-skip form walks exactly the UNION of
        // the two rows' pinned units (the skipped count is the
        // brute-force complement), and its sum is BIT-IDENTICAL to
        // the committed accumulation whenever the union covers
        // every unit (the constant-mass term is exactly zero); when
        // constant entries are skipped, the regrouped accumulation
        // is allowed its named small ulp delta (never more than a
        // few ulps at this toy scale) — the delta is the honest
        // cost of the fast mode, measured not hidden.
        let elsewhere = -14.418156f64;
        // ulp_distance's own contract: equal bits are 0, adjacent
        // floats are 1.
        assert_eq!(ulp_distance(1.0, 1.0), 0);
        assert_eq!(
            ulp_distance(f64::from_bits(1.0f64.to_bits() + 1), 1.0),
            1
        );
        // (folds x units, the pin STRIDE; 0 = no pins) shapes: a
        // no-pin-heavy row (Sparse encoding), a mid row (Sparse), a
        // pinned-heavy row (Dense encoding — exercises the built-list
        // path), and an all-no-pin row.
        let shapes = [(5usize, 40usize, 8usize), (5, 40, 2), (5, 40, 1), (5, 40, 0)];
        for (n_folds, n_units, pin_stride) in shapes {
            let mut dense: Vec<Vec<f64>> = Vec::new();
            let mut pinned: Vec<Vec<bool>> = Vec::new();
            let mut mismatches: Vec<Vec<u32>> = Vec::new();
            for fold in 0..n_folds {
                let mut row = vec![elsewhere; n_units];
                let mut pin = vec![false; n_units];
                let mut mism = vec![0u32; n_units];
                if pin_stride > 0 {
                    for unit in (fold..n_units).step_by(pin_stride) {
                        pin[unit] = true;
                        row[unit] = elsewhere + (unit as f64) / (fold as f64 + 1.0);
                        mism[unit] = (unit % 7) as u32;
                    }
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(mism);
            }
            let matrix = ScoreMatrix::from_dense(
                dense.clone(),
                pinned.clone(),
                mismatches.clone(),
                elsewhere,
            );
            let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
            let fast = ClassScanFastIndex::build(&matrix, &counts);
            for first in 0..n_folds {
                for second in first..n_folds {
                    // the skipped count is the brute-force complement
                    let brute_union = (0..n_units)
                        .filter(|&u| pinned[first][u] || pinned[second][u])
                        .count();
                    assert_eq!(
                        fast.constant_units_skipped(first, second),
                        (n_units - brute_union) as u64,
                        "the union walk's skipped count, ({first},{second})"
                    );
                    let exact = matrix.class_pair_sum(first, second, &counts);
                    let fast_sum = fast.class_pair_sum(first, second, &counts);
                    if brute_union == n_units {
                        assert_eq!(
                            fast_sum.to_bits(),
                            exact.to_bits(),
                            "full-union pair bit-identical, ({first},{second})"
                        );
                    } else {
                        assert!(
                            ulp_distance(exact, fast_sum) <= 64,
                            "the named delta is small at toy scale: \
                             ({first},{second}) {} ulps",
                            ulp_distance(exact, fast_sum)
                        );
                        assert!(
                            (fast_sum - exact).abs() <= 1e-12 * exact.abs().max(1.0),
                            "the named delta is a rounding-scale effect"
                        );
                    }
                }
            }
            // the all-no-pin pair: the fast sum is exactly the one
            // constant term elsewhere * total_count.
            if pin_stride == 0 {
                let total_count: f64 = counts.iter().sum();
                for first in 0..n_folds {
                    assert_eq!(
                        fast.class_pair_sum(first, first, &counts).to_bits(),
                        (elsewhere * total_count).to_bits(),
                        "the all-no-pin pair's single-term form"
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "the marathon worst-locus shape receipt (a ~3GB build; run explicitly)"]
    fn marathon_worst_locus_shape_memory_receipt() {
        // The profiled worst-locus shape (SS10.6.3's roster, T03803:
        // 532 folds x 4,438,575 units, 2.1% pinned pairs) built as
        // the compressed matrix, with the RSS delta measured from
        // /proc/self/status — the memory bound receipt. The dense
        // form's 13 B/pair (30.7GB) is exact arithmetic, never built.
        fn vm_rss_kb() -> u64 {
            std::fs::read_to_string("/proc/self/status")
                .unwrap()
                .lines()
                .find(|line| line.starts_with("VmRSS:"))
                .and_then(|line| {
                    line.split_whitespace().nth(1).and_then(|n| n.parse().ok())
                })
                .unwrap_or(0)
        }
        let (n_folds, n_units, pinned_stride) = (532usize, 4_438_575usize, 48usize);
        let elsewhere = -14.418156f64;
        let before = vm_rss_kb();
        let mut rows = Vec::with_capacity(n_folds);
        for fold in 0..n_folds {
            // a sparse row at the roster's own pinned share: every
            // 48th unit pins (~2.1%)
            let mut units = Vec::new();
            let mut lls = Vec::new();
            let mut mism = Vec::new();
            let mut unit = fold % pinned_stride;
            while unit < n_units {
                units.push(unit as u32);
                lls.push(elsewhere + (unit % 977) as f64 * 1e-6);
                mism.push((unit % 7) as u32);
                unit += pinned_stride;
            }
            rows.push(ScoreRow::Sparse(ScoreRowSparse {
                units,
                lls,
                mismatches: mism,
            }));
        }
        let matrix = ScoreMatrix {
            n_units,
            elsewhere,
            rows,
            unit_fold: None,
        };
        let after = vm_rss_kb();
        let live = matrix.live_bytes();
        let pairs = n_folds as u64 * n_units as u64;
        let dense_bytes = pairs * 13;
        eprintln!(
            "marathon shape T03803: {} folds x {} units = {} pairs; \\
             compressed live {} MB ({:.2} B/pair); the dense form's \\
             bound {} MB (13 B/pair); RSS delta {} MB; collapse {:.1}x",
            n_folds,
            n_units,
            pairs,
            live / (1024 * 1024),
            live as f64 / pairs as f64,
            dense_bytes / (1024 * 1024),
            (after.saturating_sub(before)) / 1024,
            dense_bytes as f64 / live as f64,
        );
        // the bound: the compressed form is a small fraction of the
        // dense 13 B/pair at the roster's own no-pin share, and the
        // measured RSS delta is the live bytes' own magnitude class
        // (allocator overhead included).
        assert!(live * 20 < dense_bytes, "the compression bound");
        assert!(
            (after.saturating_sub(before)) * 1024 < live * 3,
            "the RSS delta is the live bytes' class"
        );
        // the pinned-majority control (T03641's class, 84% pinned):
        // the DERIVED choice keeps the committed dense size EXACTLY
        // (never larger) — 16 B/pinned >= 13 B/unit chooses dense.
        let pinned_majority = 0.85 * n_units as f64;
        assert!((16.0 * pinned_majority) >= (13.0 * n_units as f64));
    }

    #[test]
    #[ignore = "the marathon core scan-wall receipt at the profiled worst-locus shape (~1-2 min on a 64-thread pool; run explicitly)"]
    fn marathon_worst_locus_shape_scan_wall_receipt() {
        // THE MARATHON RE-RECEIPT'S MEASURABLE PART (SS10.6.3's
        // profile: the marathon substrates are reaped, so the full
        // marathon wall is NOT re-receiptable — this receipts what the
        // profiled worst-locus shape can still exercise, and names
        // what it cannot): the scoring core's own scan walls at the
        // T03803 shape (532 folds x 4,438,575 units, the roster's own
        // ~2.1% pinned share, sparse rows — candidate C's compressed
        // form), measured on this test's rayon pool:
        //   (a) the merge's equal-column hash pass (ll_bits_hash over
        //       every row — the bucketing key of the committed
        //       merge),
        //   (b) the exhaustive class scan (every fold pair's
        //       class_pair_sum — the committed accumulation verbatim),
        // both as ns per (fold, unit) pair walked. NOT measured here,
        // because the reaped substrates carried them: the pair loop's
        // PINNED alignment class (records + reads), the bind phase,
        // the records/census phases — no marathon wall number is
        // claimed from this receipt, only the scan-core rates.
        let (n_folds, n_units, pinned_stride) = (532usize, 4_438_575usize, 48usize);
        let elsewhere = -14.418156f64;
        let mut rows = Vec::with_capacity(n_folds);
        for fold in 0..n_folds {
            let mut units = Vec::new();
            let mut lls = Vec::new();
            let mut mism = Vec::new();
            let mut unit = fold % pinned_stride;
            while unit < n_units {
                units.push(unit as u32);
                lls.push(elsewhere + (unit % 977) as f64 * 1e-6);
                mism.push((unit % 7) as u32);
                unit += pinned_stride;
            }
            rows.push(ScoreRow::Sparse(ScoreRowSparse {
                units,
                lls,
                mismatches: mism,
            }));
        }
        let matrix = ScoreMatrix {
            n_units,
            elsewhere,
            rows,
            unit_fold: None,
        };
        // (a) the merge's hash pass.
        let hash_started = std::time::Instant::now();
        let mut hash_sum: u64 = 0;
        for fold in 0..n_folds {
            hash_sum ^= matrix.ll_bits_hash(fold);
        }
        let hash_seconds = hash_started.elapsed().as_secs_f64();
        // (b) the exhaustive class scan, the committed shape: every
        // fold pair, unit counts all 1.0 (the synthetic shape carries
        // no real read mass; the count is one FMA per unit and does
        // not change the rate).
        let unit_counts: Vec<f64> = vec![1.0; n_units];
        let pair_count = n_folds * (n_folds + 1) / 2;
        let scan_started = std::time::Instant::now();
        let class_lls: Vec<f64> = (0..pair_count)
            .into_par_iter()
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                matrix.class_pair_sum(first, second, &unit_counts)
            })
            .collect();
        let scan_seconds = scan_started.elapsed().as_secs_f64();
        let walked = n_folds as u64 * n_units as u64;
        let scanned_pairs = pair_count as u64 * n_units as u64;
        let best = class_lls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "marathon shape T03803 scan-wall receipt: merge hash pass \
             {} rows in {:.1}s = {:.2} ns/(fold,unit); class scan \
             {} pairs x {} units = {:.3}e12 unit-pairs in {:.1}s \
             ({} threads) = {:.2} ns/unit-pair; best class ll {:.3}",
            n_folds,
            hash_seconds,
            hash_seconds * 1e9 / walked as f64,
            pair_count,
            n_units,
            scanned_pairs as f64 / 1e12,
            scan_seconds,
            rayon::current_num_threads(),
            scan_seconds * 1e9 / scanned_pairs as f64,
            best,
        );
        assert!(best.is_finite(), "the class scan produced a finite best");
        let _ = hash_sum;
    }

    #[test]
    fn the_alignment_connected_merge_folds_indistinguishable_copies() {
        // The axis fold (0) and two near-identical copies (1, 2) share
        // nodes (alignment-connected), carry bit-identical likelihood
        // columns, and their variant pockets carry ZERO observed mass
        // — they fold. Fold 3 differs from 1 only in a pocket node
        // (5) that carries observed mass — its merge is REFUSED (the
        // safety property: a pocket with even one observed read
        // refuses the fold). Fold 4 has a distinguishable column —
        // refused. Fold 5 shares NO node with the bubble — outside
        // the axis's connected alignment component, refused.
        let folds = vec![
            // 0: the axis correspondence.
            alignfold_toy(&[(1, 1), (2, 1)], &[(100, 1)], &["axis"]),
            // 1: a near-identical copy — pocket node 3 (massless).
            alignfold_toy(&[(1, 1), (2, 1), (3, 1)], &[(100, 1), (101, 1)], &["twinA"]),
            // 2: another — pocket node 4 (massless).
            alignfold_toy(&[(1, 1), (2, 1), (4, 1)], &[(100, 1), (102, 1)], &["twinB"]),
            // 3: differs from 1 by pocket node 5 (OBSERVED mass 7.5).
            alignfold_toy(&[(1, 1), (2, 1), (5, 1)], &[(100, 1), (103, 1)], &["massed"]),
            // 4: shares nodes, unobserved pockets, but a DIFFERENT
            // column (a read-level distinction).
            alignfold_toy(&[(1, 1), (2, 1), (6, 1)], &[(100, 1), (104, 1)], &["othercolumn"]),
            // 5: NO shared node with the bubble (its own alignment
            // component) — outside the axis correspondence's bubble.
            alignfold_toy(&[(30, 1), (31, 1)], &[(105, 1)], &["disconnected"]),
        ];
        let matrix = vec![
            vec![1.5, -2.0],
            vec![1.5, -2.0],
            vec![1.5, -2.0],
            vec![1.5, -2.0],
            vec![1.5, 0.25], // fold 4: a different column
            vec![1.5, -2.0],
        ];
        let pinned = vec![vec![true, false]; 6];
        let mismatches = vec![vec![0u32, 0u32]; 6];
        let territory = vec![(0u64, 4u64); 6];
        let observed_node: HashMap<u64, f64> =
            [(5u64, 7.5f64)].into_iter().collect();
        let observed_edge: HashMap<u64, f64> = HashMap::new();
        let matrix = ScoreMatrix::from_dense(matrix, pinned, mismatches, -14.0);
        let (new_folds, new_matrix, new_territory, new_axis, groups, counters) =
            alignment_connected_merge(
                folds, matrix, territory, Some(0), Some(0),
                &observed_node, &observed_edge,
            );
        // Exactly one merge group: {0, 1, 2} (the axis bubble's
        // indistinguishable copies).
        assert_eq!(groups.len(), 1, "one merge group");
        assert_eq!(groups[0].members, vec![0, 1, 2]);
        assert_eq!(groups[0].representative, 0);
        assert_eq!(groups[0].pocket_node_mass, 0.0);
        assert_eq!(groups[0].pocket_edge_mass, 0.0);
        // 0, 1, 2 coalesce; 3, 4, 5 stand: 4 folds.
        assert_eq!(new_folds.len(), 4);
        // The merged fold carries every member's rows and the
        // REPRESENTATIVE's multiset (one hypothesis's material —
        // members are proven equivalents, not additional copies).
        assert_eq!(new_folds[0].members.len(), 3);
        assert_eq!(new_folds[0].nodes, vec![(1, 1), (2, 1)]);
        assert_eq!(new_folds[0].edges, vec![(100, 1)]);
        // The merged fold's column is the representative's.
        assert_eq!(new_matrix.ll(0, 0), 1.5);
        assert_eq!(new_matrix.ll(0, 1), -2.0);
        assert_eq!(new_territory[0], (0, 4));
        assert_eq!(new_axis, Some(0));
        // The refusal anatomy: fold 3 refused on pocket mass (with the
        // axis too — its pocket 5 is observed); fold 4 on the column;
        // fold 5 on the bubble (4 pairs). The equal-column bucket is
        // {0,1,2,3,5}: 10 candidate pairs.
        assert_eq!(counters.equal_column_pairs, 10);
        assert_eq!(counters.refused_pocket_mass, 3); // (0,3),(1,3),(2,3)
        assert_eq!(counters.refused_outside_bubble, 4); // (x,5) pairs
        assert_eq!(counters.refused_unshared, 0);
        assert_eq!(counters.merged_pairs, 3); // (0,1),(0,2),(1,2)
        // Solo folds keep their own data.
        assert_eq!(new_folds[1].members[0].path_name, "massed");
    }

    #[test]
    fn the_alignment_connected_merge_refuses_bit_differences() {
        // Two alignment-adjacent folds with massless pockets whose
        // columns differ in the LAST BIT: not bit-identical — a
        // candidate can distinguish them (the exactness discipline
        // is bit-level) — the fold does NOT happen.
        let folds = vec![
            alignfold_toy(&[(1, 1)], &[], &["axis"]),
            alignfold_toy(&[(1, 1), (2, 1)], &[], &["twin"]),
        ];
        let matrix = vec![
            vec![1.0f64],
            vec![1.0f64 + f64::EPSILON], // one ULP apart
        ];
        let matrix =
            ScoreMatrix::from_dense(matrix, vec![vec![true]; 2], vec![vec![0u32]; 2], -14.0);
        let (new_folds, _m, _t, _a, groups, counters) =
            alignment_connected_merge(
                folds, matrix, vec![(0, 4); 2], Some(0), Some(0),
                &HashMap::new(), &HashMap::new(),
            );
        assert!(groups.is_empty(), "a ULP distinction refuses the fold");
        assert_eq!(new_folds.len(), 2);
        assert_eq!(counters.equal_column_pairs, 0);
        assert_eq!(counters.merged_pairs, 0);
    }

    #[test]
    fn the_alignment_connected_merge_without_axis_fold_is_a_noop() {
        // No axis correspondence: no bubble, no merges — the exact
        // fold space passes through unchanged.
        let folds = vec![
            alignfold_toy(&[(1, 1)], &[], &["a"]),
            alignfold_toy(&[(1, 1)], &[], &["b"]),
        ];
        let matrix = vec![vec![1.0], vec![1.0]];
        let matrix =
            ScoreMatrix::from_dense(matrix, vec![vec![true]; 2], vec![vec![0u32]; 2], -14.0);
        let (new_folds, _m, _t, _a, groups, _c) = alignment_connected_merge(
            folds, matrix, vec![(0, 4); 2], None, None,
            &HashMap::new(), &HashMap::new(),
        );
        assert!(groups.is_empty());
        assert_eq!(new_folds.len(), 2);
    }

    #[test]
    fn the_connected_truth_second_takes_the_max_shared_correspondence() {
        // The segdup block's attribution: several folds carry the
        // second haplotype's rows; the alignment connection (the
        // maximum shared contained-node mass with the axis fold)
        // selects the correspondence; ties and disconnected
        // candidates refuse (no fabrication).
        let axis = alignfold_toy(&[(1, 1), (2, 1), (3, 1), (4, 1)], &[], &["axis"]);
        let near = alignfold_toy(&[(1, 1), (2, 1), (3, 1), (9, 1)], &[], &["HG002#2"]);
        let mid = alignfold_toy(&[(1, 1), (2, 1)], &[], &["HG002#2"]);
        let far = alignfold_toy(&[(1, 1)], &[], &["HG002#2"]);
        let unconnected = alignfold_toy(&[(40, 1)], &[], &["HG002#2"]);
        let folds = vec![axis, mid, near, far, unconnected];
        let (chosen, receipt) =
            connected_truth_second(&folds, Some(0), "HG002#2");
        assert_eq!(chosen, Some(2), "the maximum shared-node mass wins");
        let receipt = receipt.expect("the correspondence receipt");
        let candidates = receipt["candidates"].as_array().unwrap();
        assert_eq!(candidates.len(), 4, "every carrier named");
        assert_eq!(candidates[0]["fold"], 2);
        assert_eq!(candidates[0]["shared_node_mass"], 3);

        // A TIE at the maximum: the ambiguous correspondence refuses.
        let tied = alignfold_toy(&[(1, 1), (2, 1), (3, 1), (10, 1)], &[], &["HG002#2"]);
        let folds = vec![
            alignfold_toy(&[(1, 1), (2, 1), (3, 1), (4, 1)], &[], &["axis"]),
            alignfold_toy(&[(1, 1), (2, 1), (3, 1), (9, 1)], &[], &["HG002#2"]),
            tied,
        ];
        let (chosen, _receipt) = connected_truth_second(&folds, Some(0), "HG002#2");
        assert_eq!(chosen, None, "a tie refuses the attribution");

        // NO shared node with the axis fold: not connected to the
        // bubble — refused even as the sole carrier.
        let folds = vec![
            alignfold_toy(&[(1, 1)], &[], &["axis"]),
            alignfold_toy(&[(40, 1)], &[], &["HG002#2"]),
        ];
        let (chosen, _receipt) = connected_truth_second(&folds, Some(0), "HG002#2");
        assert_eq!(chosen, None, "a disconnected carrier refuses");

        // No axis fold: no correspondence at all.
        let (chosen, receipt) = connected_truth_second(&folds, None, "HG002#2");
        assert_eq!(chosen, None);
        assert!(receipt.is_none());
    }

    // ------------------------------------------- the official sequence-QV
    // convention (the owner's ruling, the likegt max-qv form):
    // the PAF derivation, the per-hap best selection, the mean QV,
    // and — when the allwave binary is present — the convention end
    // to end on a known diplotype.

    fn paf_line(query: &str, target: &str, qs: u64, qe: u64, matches: u64) -> String {
        format!(
            "{query}\t{qe}\t{qs}\t{qe}\t+\t{target}\t{qe}\t0\t{qe}\t{matches}\t{qe}\t60\n"
        )
    }

    #[test]
    fn allwave_paf_derivation_matches_likegt_max_qv() {
        // identity = matches / query_span (the aligned query span,
        // NOT the alignment length); QV = -10 log10(1 - identity)
        // capped at 50, perfect = QV 50.0 exactly.
        let paf = paf_line("truth0", "called0", 0, 1000, 995);
        let lines = AllwaveQv::parse_paf(&paf).unwrap();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].2, 995);
        assert_eq!(lines[0].3, 1000);
        assert!((AllwaveQv::qv_of_identity(0.995) - 23.010_299_956).abs() < 1e-9);
        assert_eq!(AllwaveQv::qv_of_identity(1.0), 50.0);
        assert_eq!(AllwaveQv::qv_of_identity(0.999999), 50.0, "the QV-50 cap");
        assert!((AllwaveQv::qv_of_identity(0.5) - 3.010_299_956_639_78).abs() < 1e-9);
        assert_eq!(AllwaveQv::qv_of_identity(0.0), 0.0);
        // a malformed line is an error, never a silent skip
        assert!(AllwaveQv::parse_paf("truth0\tshort\n").is_err());
    }

    #[test]
    fn allwave_per_hap_best_pools_and_selects_like_max_qv() {
        // truth0 -> called0: two alignment lines pooling to
        // (matches 990, span 1000) = identity 0.99; truth0 ->
        // called1: one line identity 0.995. The BEST alignment for
        // the haplotype is called1 (the likegt most-similar rule).
        let paf_text = format!(
            "{}{}{}",
            paf_line("truth0", "called0", 0, 500, 495),
            paf_line("truth0", "called0", 500, 1000, 495),
            paf_line("truth0", "called1", 0, 1000, 995),
        );
        let paf = AllwaveQv::parse_paf(&paf_text).unwrap();
        let (identity, matches, span, target) =
            AllwaveQv::best_for_hap(&paf, "truth0", ["called0", "called1"]);
        assert!((identity - 0.995).abs() < 1e-12);
        assert_eq!((matches, span), (995, 1000));
        assert_eq!(target.as_deref(), Some("called1"));
        // no alignment at all: identity 0, QV 0 (the honest floor —
        // zero matched truth material, never a silent skip).
        let (identity, _m, _s, target) =
            AllwaveQv::best_for_hap(&paf, "truth9", ["called0", "called1"]);
        assert_eq!(identity, 0.0);
        assert!(target.is_none());
        // an exact TIE between the called homologs keeps the FIRST
        // (the stated deterministic order).
        let tied = format!("{}{}", paf_line("t", "called0", 0, 100, 99), paf_line("t", "called1", 0, 100, 99));
        let paf = AllwaveQv::parse_paf(&tied).unwrap();
        let (_, _, _, target) = AllwaveQv::best_for_hap(&paf, "t", ["called0", "called1"]);
        assert_eq!(target.as_deref(), Some("called0"));
        // the diplotype QV is the MEAN of the per-hap QVs (likegt's
        // avg_qv form), and a perfect call under the official
        // convention is QV 50 (the cap), not the biWFA receipts' 60.
        let perfect = format!(
            "{}{}",
            paf_line("truth0", "called0", 0, 100, 100),
            paf_line("truth1", "called1", 0, 100, 99),
        );
        let paf = AllwaveQv::parse_paf(&perfect).unwrap();
        let (id0, ..) = AllwaveQv::best_for_hap(&paf, "truth0", ["called0", "called1"]);
        let (id1, ..) = AllwaveQv::best_for_hap(&paf, "truth1", ["called0", "called1"]);
        let mean = (AllwaveQv::qv_of_identity(id0) + AllwaveQv::qv_of_identity(id1)) / 2.0;
        assert!((mean - 35.0).abs() < 1e-9, "(50 + 20) / 2, QV of 99% identity is 20");
    }

    #[test]
    fn allwave_official_convention_end_to_end_on_a_known_diplotype() {
        // The convention live through the allwave binary: four
        // ~2000bp sequences with KNOWN mismatch counts; skipped
        // only when the allwave binary is absent from PATH (the
        // hermetic derivation is covered by the tests above).
        let binary = std::env::var_os("IMPG_TRUTH_QV_ALLWAVE_BIN")
            .unwrap_or_else(|| "allwave".into());
        let probe = std::process::Command::new(&binary)
            .arg("--version")
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .output();
        if probe.is_err() {
            eprintln!("allwave absent: the end-to-end official-convention test is skipped");
            return;
        }
        let mut base = vec![b'A'; 2000];
        for (i, b) in base.iter_mut().enumerate() {
            *b = b"ACGT"[i % 4];
        }
        // mutate a DISJOINT position set: every position truly changes
        // base letter (the +1 map A->C, C->G, G->T, T->A), so the
        // mismatch counts are exact.
        let mutated = |offset: usize, n: usize| {
            let mut seq = base.clone();
            for i in 0..n {
                let position = offset + 37 * i;
                seq[position] = b"ACGT"[(seq[position] - b'A') as usize % 3 + 1];
            }
            seq
        };
        let truth0 = base.clone();
        let truth1 = mutated(500, 20);      // 20 mismatches vs truth0
        let called0 = mutated(0, 5);        // 5 vs truth0
        // called1 = truth1 + 20 more mutations: 20 mismatches vs
        // truth1, 40 vs truth0 — so truth1's best is called1 and
        // truth0's best is called0.
        let mut called1 = truth1.clone();
        for i in 0..20 {
            let position = 1200 + 37 * i;
            called1[position] = b"ACGT"[(called1[position] - b'A') as usize % 3 + 1];
        }
        let allwave = AllwaveQv { binary, threads: 4 };
        let record = allwave
            .locus_qv(4242, &called0, &called1, &truth0, &truth1)
            .unwrap();
        // truth0 (the unmutated base) takes called0 (5 mismatches,
        // identity 1995/2000); truth1 takes called1 (20 mismatches,
        // identity 1980/2000).
        assert!((record["identity_truth_hap0"].as_f64().unwrap() - 1995.0 / 2000.0).abs() < 1e-9);
        assert!((record["identity_truth_hap1"].as_f64().unwrap() - 1980.0 / 2000.0).abs() < 1e-9);
        assert_eq!(record["best_called_hap_truth_hap0"], "called0");
        assert_eq!(record["best_called_hap_truth_hap1"], "called1");
        let qv0 = -10.0 * (5.0f64 / 2000.0).log10();
        let qv1 = -10.0 * (20.0f64 / 2000.0).log10();
        assert!(
            (record["qv"].as_f64().unwrap() - (qv0 + qv1) / 2.0).abs() < 1e-6,
            "the diplotype QV is the mean of the per-hap QVs"
        );
        assert_eq!(record["perfect"], false);
        // a PERFECT call: both truth haps identical to their called
        // homolog — identity 1.0, QV 50.0 exactly (the official cap).
        let perfect = allwave
            .locus_qv(4243, &base, &base, &base.clone(), &base.clone())
            .unwrap();
        assert_eq!(perfect["identity_truth_hap0"], 1.0);
        assert_eq!(perfect["identity_truth_hap1"], 1.0);
        assert_eq!(perfect["qv"], 50.0);
        assert_eq!(perfect["perfect"], true);
    }

    #[test]
    fn the_sequence_emission_spells_ascii_and_fails_closed() {
        // The product emission's material: the fold's spelled
        // sequence as a JSON-safe string, fail-closed on non-ASCII.
        let spelled = b"ACGTNacgtn".as_slice();
        assert_eq!(
            fold_sequence_string(spelled).unwrap(),
            "ACGTNacgtn".to_string()
        );
        assert_eq!(fold_sequence_string(&[]).unwrap(), "");
        let mut poison = b"ACGT".to_vec();
        poison.push(0xC3);
        poison.push(0xA9); // UTF-8 continuation bytes are non-ASCII
        assert!(fold_sequence_string(&poison).is_err());
    }

    #[test]
    fn the_flank_overlap_parses_both_forms() {
        // The bp form and the percent form; the loud rejections.
        assert_eq!("10000".parse::<FlankOverlap>(), Ok(FlankOverlap::Bp(10000)));
        assert_eq!("10%".parse::<FlankOverlap>(), Ok(FlankOverlap::Percent(10.0)));
        assert_eq!("100%".parse::<FlankOverlap>(), Ok(FlankOverlap::Percent(100.0)));
        assert!("0".parse::<FlankOverlap>().is_err());
        assert!("0%".parse::<FlankOverlap>().is_err());
        assert!("101%".parse::<FlankOverlap>().is_err());
        assert!("abc".parse::<FlankOverlap>().is_err());
        assert!("".parse::<FlankOverlap>().is_err());
        assert_eq!(FlankOverlap::Bp(10000).setting(), "10000");
        assert_eq!(FlankOverlap::Percent(10.0).setting(), "10%");
    }

    #[test]
    fn the_flank_rule_keeps_only_rows_reaching_into_both_flanks() {
        // THE PRESENCE SEMANTICS (ruled at the measured fork): the
        // window [1000, 1100) (100bp); N = 10bp - a member row is kept
        // iff its hull OVERLAPS [1000, 1010) AND [1090, 1100).
        let flank = FlankOverlap::Bp(10);
        // full coverage: kept
        assert!(row_covers_flanks(1000, 1100, 1000, 1100, flank));
        assert!(row_covers_flanks(900, 1200, 1000, 1100, flank));
        // rows reaching into both flanks (edge truncation is not
        // fragmentation - the hull-edge geometry): kept
        assert!(row_covers_flanks(1005, 1100, 1000, 1100, flank));
        assert!(row_covers_flanks(1000, 1095, 1000, 1100, flank));
        assert!(row_covers_flanks(1000, 1091, 1000, 1100, flank));
        // missing the first flank entirely (starts after it ends):
        // removed
        assert!(!row_covers_flanks(1010, 1100, 1000, 1100, flank));
        assert!(!row_covers_flanks(1095, 1200, 1000, 1100, flank));
        // missing the last flank entirely (ends before it starts):
        // removed
        assert!(!row_covers_flanks(1000, 1090, 1000, 1100, flank));
        // an interior row missing both flanks: removed
        assert!(!row_covers_flanks(1020, 1080, 1000, 1100, flank));
        // a disjoint row: removed
        assert!(!row_covers_flanks(2000, 2100, 1000, 1100, flank));
    }

    #[test]
    fn the_flank_percent_form_is_ceil_of_the_window_length() {
        // 10% of a 10,036bp window is 1,004bp (ceil); 10% of 100bp is
        // exactly 10bp.
        assert_eq!(FlankOverlap::Percent(10.0).flank_bp(10_036), 1_004);
        assert_eq!(FlankOverlap::Percent(10.0).flank_bp(100), 10);
        assert_eq!(FlankOverlap::Percent(50.0).flank_bp(101), 51);
        // at 100% the flank intervals span the whole window: the rule
        // is presence in the window.
        assert!(row_covers_flanks(1, 100, 0, 100, FlankOverlap::Percent(100.0)));
        assert!(!row_covers_flanks(100, 200, 0, 100, FlankOverlap::Percent(100.0)));
    }

    #[test]
    fn the_flank_rule_degenerates_to_presence_below_n() {
        // THE DEGENERATE CASE (the owner's rule): a window shorter
        // than N keeps every PRESENT row - even an interior fragment.
        let flank = FlankOverlap::Bp(1000);
        assert!(row_covers_flanks(20, 80, 0, 100, flank));
        assert!(row_covers_flanks(0, 10, 0, 100, flank));
        assert!(!row_covers_flanks(200, 300, 0, 100, flank));
        // a window exactly N long: both flank intervals span the
        // whole window - presence keeps every overlapping row (not
        // the degenerate case - "shorter than N" - the same result
        // either way).
        assert!(row_covers_flanks(0, 1000, 0, 1000, flank));
        assert!(row_covers_flanks(1, 999, 0, 1000, flank));
        assert!(!row_covers_flanks(1000, 2000, 0, 1000, flank));
    }

    #[test]
    fn the_exclude_paths_glob_matches_fnmatch_style() {
        // The masking operator's matcher: `*` any run, `?` one char,
        // literal otherwise; no escapes.
        assert!(path_matches_glob("HG002#1#chr1:1-100", "HG002#*"));
        assert!(path_matches_glob("HG002#1#chr1:1-100", "*#chr1*"));
        assert!(path_matches_glob("HG002#1#chr1:1-100", "HG002#?#chr1:1-100"));
        assert!(path_matches_glob("HG002#1#chr1:1-100", "HG002#1#chr1:1-100"));
        assert!(path_matches_glob("anything", "*"));
        assert!(path_matches_glob("", "*"));
        assert!(!path_matches_glob("HG002#1#chr1:1-100", "HG003*"));
        // `?` consumes exactly one char: it matches the final `0`
        // here, but cannot span two.
        assert!(path_matches_glob("HG002#1#chr1:1-100", "HG002#1#chr1:1-10?"));
        assert!(!path_matches_glob("HG002#1#chr1:1-100", "HG002#1#chr1:1-?"));
        assert!(!path_matches_glob("HG002#1#chr1:1-100", "HG002#?"));
        // a pattern ending in `*` must be able to consume the whole
        // name (the `(1..=len).rev()` arm) and also match nothing.
        assert!(path_matches_glob("abc", "abc*"));
        assert!(path_matches_glob("abc", "a*bc*"));
        assert!(!path_matches_glob("abc", "b*"));
    }

    #[test]
    fn the_exclusion_never_drops_the_axis_row() {
        // The coordinate frame stays in the candidate universe even
        // when the glob matches it; every other matching row leaves.
        let component = "HG002#1#chr1:0-1000";
        let globs = vec!["HG002#*".to_string()];
        assert!(!row_is_excluded(component, component, &globs));
        assert!(row_is_excluded("HG002#2#chr1:0-1000", component, &globs));
        assert!(!row_is_excluded("CHM13#0#chr1:0-1000", component, &globs));
        // no globs: nothing is excluded (the default, the committed
        // candidate universe verbatim).
        assert!(!row_is_excluded("HG002#2#chr1:0-1000", component, &[]));
    }

    #[test]
    fn the_region_resolves_to_the_covered_windows() {
        // The component's own frame: HG002#1#chr1:4632437-4883593,
        // windows of 10kbp.
        let component = "HG002#1#chr1:4632437-4883593";
        let rows: Vec<(u64, u64)> = (0..25).map(|i| (i * 10_000, (i + 1) * 10_000)).collect();
        // a region covering exactly window 1, then windows 1..2
        // (half-open bounds never bleed into the neighbors)
        let loci = resolve_region_loci(component, &rows, "chr1:4642437-4652437").unwrap();
        assert_eq!(loci, vec![1]);
        let loci = resolve_region_loci(component, &rows, "chr1:4642437-4662437").unwrap();
        assert_eq!(loci, vec![1, 2]);
        // the whole frame: every window
        let all = resolve_region_loci(
            component,
            &rows,
            "chr1:4632437-4883593",
        )
        .unwrap();
        assert_eq!(all, (0..25u32).collect::<Vec<_>>());
        // a region touching only a boundary point covers nothing
        // (half-open overlap)
        assert!(resolve_region_loci(component, &rows, "chr1:4640000-4640000").is_err());
        // empty region, outside the frame, wrong chromosome, and
        // malformed forms all fail closed
        assert!(resolve_region_loci(component, &rows, "chr1:4640000-4640000x").is_err());
        assert!(resolve_region_loci(component, &rows, "chr1:1-100").is_err());
        assert!(resolve_region_loci(component, &rows, "chr2:1-100").is_err());
        assert!(resolve_region_loci(component, &rows, "chr1").is_err());
        // a component name without a coordinate frame fails closed
        assert!(resolve_region_loci("HG002#1", &rows, "chr1:1-100").is_err());
        // a region past the component end fails closed
        assert!(resolve_region_loci(
            component,
            &rows,
            "chr1:4883593-4883594"
        )
        .is_err());
    }

    #[test]
    fn the_component_axis_rows_load_from_the_locality_maps() {
        // The axis rows: the component's member rows over the
        // locality maps, in axis order (start), from real map files.
        let dir = std::env::temp_dir().join("impg-test-region-axis-rows");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for (partition, start) in [(3u32, 20_000u64), (1, 0u64), (2, 10_000u64)] {
            let json = format!(
                "{{\"partition\": {partition}, \"members\": [                  {{\"path_name\": \"HG002#1#chr1:0-30000\", \"start\": {start}, \"end\": {}}},                  {{\"path_name\": \"other#1#chr1:0-30000\", \"start\": 0, \"end\": 100}}]}}",
                start + 10_000
            );
            std::fs::write(
                dir.join(format!("window{partition}.gfa.map.json")),
                json,
            )
            .unwrap();
        }
        let rows = component_axis_rows(&dir, "HG002#1#chr1:0-30000").unwrap();
        assert_eq!(rows, vec![(0, 10_000), (10_000, 20_000), (20_000, 30_000)]);
        // a component absent from the maps fails closed
        assert!(component_axis_rows(&dir, "absent#1#chr1:0-30000").is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_coalesced_called_set_is_the_distinct_sorted_called_folds() {
        // The owner's ambiguity ruling: a tied-class call emits
        // every called fold's sequence; an untied call's set is
        // exactly the representative pair.
        assert_eq!(coalesced_called_set(&[(3, 1)]), vec![1, 3]);
        // Two tied classes sharing hap folds coalesce to the
        // distinct union, sorted.
        assert_eq!(
            coalesced_called_set(&[(0, 5), (0, 2), (5, 2)]),
            vec![0, 2, 5]
        );
        // A fold called on both haps of one class appears once.
        assert_eq!(coalesced_called_set(&[(7, 7)]), vec![7]);
        assert!(coalesced_called_set(&[]).is_empty());
    }

    #[test]
    fn lane_cache_span_delta_extend_serves_direct_fetch_bytes() {
        // THE DELTA-EXTEND CONTRACT, unit-proved: every served crop is
        // byte-equal to the direct Sources::fetch of the crop, through
        // a union-growth request order that forces the left-extend,
        // right-extend, and both-side extension paths — and the
        // counters receipt the DELTA fetch volume (the missing
        // sub-ranges only), never the regrown span's.
        let len = 40_000u64;
        let seq: Vec<u8> = (0..len)
            .map(|i| b"ACGT"[(i.wrapping_mul(2_654_435_761u64) % 4) as usize])
            .collect();
        let dir = std::env::temp_dir().join(format!("impg-lanecache-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("lane.fa");
        {
            use std::io::Write;
            let mut out = std::fs::File::create(&fa).unwrap();
            writeln!(out, ">lane0").unwrap();
            for chunk in seq.chunks(60) {
                out.write_all(chunk).unwrap();
                out.write_all(b"\n").unwrap();
            }
        }
        let sources = routes::Sources::open(
            &[fa.display().to_string()],
            vec![("lane0".into(), len)],
        )
        .unwrap();
        // The served-crop property itself, first: a crop of a fetch
        // answer equals the direct fetch of the crop (the property the
        // covered branch and the stitch both ride).
        let whole = sources.fetch(0, 5_000, 30_000).unwrap();
        assert_eq!(whole.len(), 25_000);
        assert_eq!(&whole[..3_000], &sources.fetch(0, 5_000, 8_000).unwrap()[..]);
        assert_eq!(&whole[20_000..], &sources.fetch(0, 25_000, 30_000).unwrap()[..]);
        let cache = LaneCache {
            lane_lengths: vec![len],
            sources: std::sync::Arc::new(sources),
            whole_lane: false,
            seqs: vec![std::sync::Mutex::new(None)],
            loads: AtomicU64::new(0),
            bytes: AtomicU64::new(0),
            nanos: AtomicU64::new(0),
            lane_trace: Default::default(),
        };
        // (initial span) right-extend, both-side extend, covered,
        // right-extend again — every crop >1,024bp (outside the small-
        // crop fetch memo; the memo path returns the same bytes either
        // way).
        for (start, end) in [
            (1_000u64, 5_000u64),
            (20_000, 25_000),
            (500, 30_000),
            (1_200, 4_800),
            (30_000, 35_000),
        ] {
            let served = cache.fetch(0, start, end).unwrap();
            assert_eq!(
                served,
                cache.sources.fetch(0, start, end).unwrap(),
                "served crop {start}..{end} differs from the direct fetch"
            );
        }
        // The delta-volume receipt: loads = 1 initial + 1 right-extend
        // ([5_000, 25_000)) + 2 both-side ([500, 1_000) + [25_000,
        // 30_000)) + 1 right-extend ([30_000, 35_000)) = 5; bytes = the
        // missing sub-ranges only = 4_000 + 20_000 + 500 + 5_000 +
        // 5_000 = 34_500 (the union-regrowth form would have fetched
        // 4_000 + 24_000 + 29_500 + 0 + 34_500 = 92_000).
        assert_eq!(cache.loads.load(Ordering::Relaxed), 5);
        assert_eq!(cache.bytes.load(Ordering::Relaxed), 34_500);
        // The whole-lane policy stays verbatim: one load, the whole
        // lane, every crop a slice of it.
        let whole_cache = LaneCache {
            lane_lengths: vec![len],
            sources: std::sync::Arc::new(
                routes::Sources::open(
                    &[fa.display().to_string()],
                    vec![("lane0".into(), len)],
                )
                .unwrap(),
            ),
            whole_lane: true,
            seqs: vec![std::sync::Mutex::new(None)],
            loads: AtomicU64::new(0),
            bytes: AtomicU64::new(0),
            nanos: AtomicU64::new(0),
            lane_trace: Default::default(),
        };
        assert_eq!(
            whole_cache.fetch(0, 1_234, 4_321).unwrap(),
            whole_cache.sources.fetch(0, 1_234, 4_321).unwrap()
        );
        assert_eq!(whole_cache.loads.load(Ordering::Relaxed), 1);
        assert_eq!(whole_cache.bytes.load(Ordering::Relaxed), len);
        let _ = std::fs::remove_file(&fa);
        let _ = std::fs::remove_file(dir.join("lane.fa.fai"));
        let _ = std::fs::remove_dir(&dir);
    }

    #[test]
    fn far_bin_prune_is_pure_window_geometry() {
        // THE OWNER'S FAR-BIN RULE, the pure-function proof: a
        // record's own routing keeps one window, or two ADJACENT
        // windows when the anchors straddle a boundary (the routing's
        // own interval-overlap rule); placements scattering into
        // NON-ADJACENT windows are the "different bins, not basically
        // adjacent" class — pruned from the census view ENTIRELY
        // (there is no choice-free home bin). The test is the
        // touched-window diameter against the census's own pair-bins
        // adjacency constant 1: no thresholds, nothing tuned.
        let line = |partitions_per_occ: &[&[u32]]| CensusRecordLine {
            record: 0,
            multiplicity: 1,
            anchors: 1,
            occurrences: partitions_per_occ
                .iter()
                .map(|&partitions| CensusOccurrence {
                    path: 0,
                    start: 0,
                    orientation: 0,
                    partitions: partitions.to_vec(),
                    intervals: Vec::new(),
                })
                .collect(),
        };
        // One window — the read's own bin: kept.
        assert!(!census_record_far_bin_scattered(&line(&[&[7]])));
        // No windows at all (votes nowhere by its own routing): not
        // the far-bin class.
        assert!(!census_record_far_bin_scattered(&line(&[&[]])));
        // The boundary-straddler: two ADJACENT bins, kept at BOTH
        // (the owner's own boundary case, pure interval overlap).
        assert!(!census_record_far_bin_scattered(&line(&[&[12, 13]])));
        // Straddling placements at separate bins of one boundary:
        // still the adjacent class.
        assert!(!census_record_far_bin_scattered(&line(&[&[12], &[13]])));
        assert!(!census_record_far_bin_scattered(&line(&[&[13, 14], &[13]])));
        // THE MEASURED FAR-BIN CLASS: the [12,16] scatter — L1's
        // 43,686-multiplicity single-anchor pattern's own shape —
        // pruned (diameter 4).
        assert!(census_record_far_bin_scattered(&line(&[&[12], &[16]])));
        // The [4,11] class (L1's second family) and the [2,13] class.
        assert!(census_record_far_bin_scattered(&line(&[&[4], &[11]])));
        assert!(census_record_far_bin_scattered(&line(&[&[2], &[13]])));
        // A mixed shape — an adjacent pair PLUS a far outlier — is
        // still the scatter class as a whole (the pattern's bins
        // are not "basically adjacent"); no per-bin surgery, no
        // choice of a home component.
        assert!(census_record_far_bin_scattered(&line(&[&[12, 13], &[16]])));
        // A three-window chain (diameter 2) is NOT basically
        // adjacent as a set: the far pair 4 and 6 rules it out.
        assert!(census_record_far_bin_scattered(&line(&[&[4], &[5], &[6]])));
        // Windows are per-occurrence in the census: the scatter is
        // visible only across occurrences, never within one
        // occurrence's own interval-overlap partitions.
        assert!(census_record_far_bin_scattered(&line(
            &[&[12, 13], &[15, 16]]
        )));
    }

    #[test]
    fn far_bin_prune_frame_gate_keeps_the_local_frame_byte_frozen() {
        // THE PRODUCT DEFAULT'S OWN GATE (the §10.6.19 default flip):
        // the prune is ON by default at the GLOBAL frame, and the
        // committed LOCAL frame is byte-frozen — a default-on prune
        // reaching the local frame is REFUSED LOUDLY (naming the
        // --far-bin-prune=false opt-out), never silently applied;
        // the explicit opt-out at the local frame is the committed
        // local form.
        assert!(far_bin_prune_frame_gate(true, true).is_ok());
        assert!(far_bin_prune_frame_gate(true, false).is_ok());
        assert!(far_bin_prune_frame_gate(false, false).is_ok());
        let refused = far_bin_prune_frame_gate(false, true).unwrap_err();
        assert_eq!(refused.kind(), std::io::ErrorKind::InvalidInput);
        assert!(refused.to_string().contains("--far-bin-prune=false"));
    }


    #[test]
    fn mix0_additive_bound_holds_for_nonnegative_excesses() {
        // THE ADMISSIBLE BOUND'S OWN LEMMA (the frontier-1 design's
        // core, SS10.6.21): for alpha, beta >= 0,
        //   mix0(alpha, beta) = ln(0.5e^alpha + 0.5e^beta) <= alpha + beta
        // — because WLOG alpha >= beta gives mix0 = alpha +
        // ln(0.5 + 0.5e^{beta - alpha}) <= alpha <= alpha + beta. This
        // is the per-unit form of the per-fold excess-mass bound
        // S(f,g) <= E*C + A_f + A_g: with every ll entry either
        // exactly E (the no-pin constant) or E + alpha (alpha > 0),
        // the class pair sum's excess over the E-floor is bounded by
        // the sum of the two folds' own excess masses. The float-level
        // check carries a DERIVED slack — the accumulation/rounding
        // slack of the two sides' own magnitudes (a few ulps), never
        // a tuned tolerance: the real-arithmetic lemma is exact.
        let e = -14.418156f64;
        // A deterministic adversarial grid + LCG sweep of excesses,
        // including the degenerate corners (0, 0), (0, big), the
        // near-equal plateau (the segdup-copy class), and the
        // far-apart class (the winner-vs-rival class).
        let mut cases: Vec<(f64, f64)> = vec![(0.0, 0.0), (0.0, 40.0), (40.0, 0.0)];
        let mut alpha = 0.0f64;
        while alpha <= 40.0 {
            let mut beta = 0.0f64;
            while beta <= 40.0 {
                cases.push((alpha, beta));
                cases.push((beta, alpha));
                beta += 0.37;
            }
            cases.push((alpha, alpha.max(1e-9)));
            alpha += 0.61;
        }
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for _ in 0..4096 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let a = (state >> 32) as f64 / u32::MAX as f64 * 44.0;
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let b = (state >> 32) as f64 / u32::MAX as f64 * 44.0;
            cases.push((a, b));
        }
        for (alpha, beta) in cases {
            let lhs = mix_logsumexp(e + alpha, e + beta) - e;
            let rhs = alpha + beta;
            // The real lemma is lhs <= rhs; each side's float
            // evaluation carries at most a few ulps of its own
            // magnitude (the mix's ln/exp rounding and the sum's).
            let slack = 8.0f64
                * f64::EPSILON
                * (lhs.abs() + rhs.abs() + e.abs()).max(1.0);
            assert!(
                lhs <= rhs + slack,
                "mix0({alpha},{beta}) = {lhs} exceeds alpha + beta = {rhs}"
            );
            // The max form is tighter and also holds (the same
            // proof): mix0 <= max(alpha, beta) + ln(1) + slack.
            let tighter = alpha.max(beta) + slack;
            assert!(
                lhs <= tighter,
                "mix0({alpha},{beta}) = {lhs} exceeds max form {tighter}"
            );
        }
    }

    #[test]
    fn exp_underflow_point_derived_by_evaluation() {
        // THE DERIVED SKIP THRESHOLD'S CONSTANT (the frontier-1
        // design): a class pair whose exact sum sits at or below
        // best + t — t the largest f64 with exp(t) == +0.0 —
        // contributes EXACTLY +0.0 to every QUAL posterior term the
        // exact run computes (exp(ll - best) underflows to the exact
        // +0.0, and adding +0.0 to the non-negative running sums is
        // float identity), so a skipped pair's QUAL arithmetic is
        // bit-identical by construction. The constant is DERIVED BY
        // EVALUATION from the f64 semantics (the subnormal underflow
        // point of exp, -1075*ln(2) in real arithmetic), never tuned:
        // walk down by ulps from the real-arithmetic point until the
        // rounded exponential is exactly zero.
        // DERIVED BY EVALUATION: the runtime's own exp is monotone,
        // so the boundary between the nonzero and the underflowed
        // region is found by a bit-axis binary search between a
        // known-nonzero point (e^-744 ~ 2^-1073.4, far above the
        // boundary) and a known-zero point (-1075*ln(2), the
        // real-arithmetic boundary itself).
        let mut lo = -744.0f64;
        let mut hi = -1075.0f64 * 2.0f64.ln();
        assert!(lo.exp() > 0.0, "the search's low anchor must not underflow");
        assert_eq!(hi.exp(), 0.0, "the search's high anchor must underflow");
        while hi.to_bits() - lo.to_bits() > 1 {
            // the midpoint on the bit axis WITHOUT the u64 sum
            // overflow (both anchors are negative doubles: their bit
            // patterns sum past u64::MAX).
            let mid = f64::from_bits(
                lo.to_bits() + (hi.to_bits() - lo.to_bits()) / 2,
            );
            if mid == lo || mid == hi {
                break;
            }
            if mid.exp() == 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let t = hi;
        // (a) the constant's own contract: exp(t) is exactly +0.0
        // and one ulp up the exponential is already nonzero.
        assert_eq!(t.exp(), 0.0);
        let up = f64::from_bits(t.to_bits() - 1);
        assert!(up.exp() > 0.0, "one ulp above the derived point must not underflow");
        // (b) the constant sits at the real-arithmetic underflow
        // point -1075*ln(2) (within the runtime exp's own boundary
        // error — a derived-scale sanity bound, not a tolerance).
        let real = -1075.0f64 * 2.0f64.ln();
        assert!(
            (t - real).abs() <= 1e-9,
            "derived point {t} vs the real-arithmetic point {real}"
        );
        // (c) monotone safety: everything at or below the point
        // underflows exactly (spot-check a ladder DOWN the bit axis
        // — for negative doubles, more-negative is the LARGER bit
        // pattern).
        let mut probe = t;
        for _ in 0..4 {
            probe = f64::from_bits(probe.to_bits() + 1);
            assert_eq!(probe.exp(), 0.0);
        }
        // (d) the QUAL identity the skip rule uses: adding the
        // underflowed +0.0 to a non-negative running sum is
        // bit-identity (including at +0.0 itself).
        let running = 3.7e-300f64;
        assert_eq!((running + 0.0f64).to_bits(), running.to_bits());
        assert_eq!((0.0f64 + 0.0f64).to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn bit_equal_rows_share_one_class_sum() {
        // THE EXACT-PRESERVING MEMOIZATION'S OWN PROOF (the only
        // frontier-1 lever that can pass the committed byte-identity
        // bar, SS10.6.21): the class pair sum is a function of the
        // rows' VALUES in unit order alone — two folds whose ll
        // columns are bit-equal (whichever encoding each row took)
        // produce BIT-IDENTICAL sums against every third fold, so
        // one computed sum per unordered (row-class, row-class) pair
        // serves every merged-fold pair without changing an output
        // bit. The construction crosses the ENCODINGS: one twin keeps
        // pinned=true at units whose value is exactly E (the sparse
        // form stores those entries), the other holds pinned=false
        // there (its sparse form drops them — the values are the
        // constant either way); the columns are bit-equal, the
        // encodings differ.
        let elsewhere = -14.418156f64;
        let n_units = 40usize;
        let mut dense: Vec<Vec<f64>> = Vec::new();
        let mut pinned: Vec<Vec<bool>> = Vec::new();
        let mut mismatches: Vec<Vec<u32>> = Vec::new();
        // Fold 0: the pinned-heavy row (Dense encoding).
        let mut row0 = vec![elsewhere; n_units];
        let mut pin0 = vec![false; n_units];
        for unit in (0..n_units).step_by(2) {
            row0[unit] = elsewhere + unit as f64 / 8.0;
            pin0[unit] = true;
        }
        dense.push(row0.clone());
        pinned.push(pin0);
        mismatches.push(vec![0; n_units]);
        // Fold 1: fold 0's VALUES with pinned=false at the constant
        // units (Sparse encoding, values unchanged).
        dense.push(row0.clone());
        pinned.push(vec![false; n_units]);
        mismatches.push(vec![0; n_units]);
        // Fold 2: an all-no-pin row (Sparse, empty).
        dense.push(vec![elsewhere; n_units]);
        pinned.push(vec![false; n_units]);
        mismatches.push(vec![0; n_units]);
        // Fold 3: bit-equal to fold 0's values but pinned=true at a
        // CONSTANT unit too (the cross-encoding case: the sparse
        // form stores the E-valued entry, the column is unchanged).
        let mut pin3 = vec![false; n_units];
        pin3[1] = true;
        dense.push(row0.clone());
        pinned.push(pin3);
        mismatches.push(vec![0; n_units]);
        let matrix = ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
        assert!(
            matrix.ll_bits_equal(0, 1),
            "the twins' ll columns are bit-equal"
        );
        assert!(
            matrix.ll_bits_equal(0, 3),
            "the cross-encoding twin's column is bit-equal"
        );
        let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
        for other in 0..matrix.n_folds() {
            assert_eq!(
                matrix.class_pair_sum(0, other, &counts).to_bits(),
                matrix.class_pair_sum(1, other, &counts).to_bits(),
                "twin rows share one sum against fold {other}"
            );
            assert_eq!(
                matrix.class_pair_sum(0, other, &counts).to_bits(),
                matrix.class_pair_sum(3, other, &counts).to_bits(),
                "the cross-encoding twin shares one sum against fold {other}"
            );
        }
        // The hom pair too (a twin against itself is trivially equal;
        // the twin's hom pair is the same accumulation).
        assert_eq!(
            matrix.class_pair_sum(0, 0, &counts).to_bits(),
            matrix.class_pair_sum(1, 1, &counts).to_bits(),
            "the twins' hom pairs share one sum"
        );
    }

    #[test]
    fn scan_memo_canonicalizes_only_bit_equal_columns() {
        // THE BIT-EQUAL-ROW MEMOIZATION'S OWN KEY DISCIPLINE (the
        // fourth expressibility reduction): the canonical row class
        // is EXACTLY the ll column (the merge's own predicate) — a
        // pinned entry whose value is exactly E is INVISIBLE to the
        // column (it reads E, like an absent entry), so rows that
        // differ only in such entries share one class and one sum;
        // a one-ulp difference at any non-E entry is a DIFFERENT
        // column — a different class, never merged (the exact
        // predicate decides every bucket candidate, so a hash
        // collision can never merge either). The representative is
        // the class's smallest fold index — deterministic.
        let elsewhere = -12.5f64;
        let n_units = 10usize;
        let mut dense: Vec<Vec<f64>> = Vec::new();
        let mut pinned: Vec<Vec<bool>> = Vec::new();
        let mut mismatches: Vec<Vec<u32>> = Vec::new();
        // fold 0: two non-E pinned entries.
        let mut row0 = vec![elsewhere; n_units];
        let mut pin0 = vec![false; n_units];
        row0[2] = elsewhere + 3.25;
        pin0[2] = true;
        row0[7] = elsewhere + 0.5;
        pin0[7] = true;
        // the working copies the later twins read (the originals move
        // into the matrix vectors below)
        let (row1, pin1) = (row0.clone(), pin0.clone());
        let (row2_src, pin2_src) = (row0.clone(), pin0.clone());
        let (row3_src, pin3_src) = (row0.clone(), pin0.clone());
        dense.push(row0);
        pinned.push(pin0);
        mismatches.push(vec![0; n_units]);
        // fold 1: an exact clone.
        dense.push(row1);
        pinned.push(pin1);
        mismatches.push(vec![0; n_units]);
        // fold 2: the clone with an EXTRA pinned entry at exactly E
        // (the cross-encoding case: the entry is invisible to the
        // column).
        let mut row2 = row2_src;
        let mut pin2 = pin2_src;
        row2[4] = elsewhere;
        pin2[4] = true;
        dense.push(row2);
        pinned.push(pin2);
        mismatches.push(vec![0; n_units]);
        // fold 3: a NEAR-TWIN — one ulp off at a non-E pinned entry.
        // A different column: its own class, never merged.
        let mut row3 = row3_src;
        let pin3 = pin3_src;
        row3[2] = f64::from_bits(row3[2].to_bits() + 1);
        dense.push(row3);
        pinned.push(pin3);
        mismatches.push(vec![0; n_units]);
        // fold 4: the all-E row, unpinned.
        dense.push(vec![elsewhere; n_units]);
        pinned.push(vec![false; n_units]);
        mismatches.push(vec![0; n_units]);
        // fold 5: the all-E row PINNED EVERYWHERE at exactly E — the
        // same column (every pinned entry is invisible).
        dense.push(vec![elsewhere; n_units]);
        pinned.push(vec![true; n_units]);
        mismatches.push(vec![0; n_units]);
        // fold 6: the all-E row again, unpinned (a third member of
        // the all-E class).
        dense.push(vec![elsewhere; n_units]);
        pinned.push(vec![false; n_units]);
        mismatches.push(vec![0; n_units]);
        let matrix = ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
        // the canonical key: the pinned-at-E entries are invisible —
        // folds 0/1/2 share one key, folds 4/5/6 share one, fold 3
        // is its own.
        assert_eq!(
            matrix.column_key_hash(0),
            matrix.column_key_hash(1),
            "exact clones share the canonical key"
        );
        assert_eq!(
            matrix.column_key_hash(0),
            matrix.column_key_hash(2),
            "an extra pinned-at-E entry is invisible to the key"
        );
        assert_ne!(
            matrix.column_key_hash(0),
            matrix.column_key_hash(3),
            "a one-ulp difference is a different key"
        );
        assert_eq!(
            matrix.column_key_hash(4),
            matrix.column_key_hash(5),
            "pinned-at-E everywhere is the all-E column"
        );
        let memo = ClassScanMemo::build(&matrix);
        assert_eq!(memo.canon, vec![0, 0, 0, 3, 4, 4, 4]);
        assert_eq!(memo.distinct_rows, 3);
        assert_eq!(memo.canonical_pairs, 6);
        // the memo's flat spelling: bit-equal pairs share one
        // canonical flat, whatever their own indices.
        assert_eq!(memo.canonical_flat(0, 3), memo.canonical_flat(2, 3));
        assert_eq!(memo.canonical_flat(1, 5), memo.canonical_flat(0, 4));
        assert_eq!(memo.canonical_flat(5, 6), memo.canonical_flat(4, 4));
        // and the shared VALUE, bit for bit, against a third fold
        // and as the hom pair.
        let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
        assert_eq!(
            matrix.class_pair_sum(0, 3, &counts).to_bits(),
            matrix.class_pair_sum(2, 3, &counts).to_bits(),
            "the cross-encoding twin shares the sum against fold 3"
        );
        assert_eq!(
            matrix.class_pair_sum(4, 4, &counts).to_bits(),
            matrix.class_pair_sum(5, 6, &counts).to_bits(),
            "the all-E class shares one hom sum"
        );
        // (The near-twin fold 3 carries a DIFFERENT column — proven
        // above by canon and the key — but its pair SUMS may still
        // collide in f64 after accumulation: the class discipline is
        // the memo's, never a sum-inequality claim.)
    }

    #[test]
    fn memoized_scan_is_bit_identical_to_direct_on_randomized_rows() {
        // THE MEMO'S OWN RANDOMIZED ORACLE (the task's own gate:
        // memoized == direct on randomized bit-row fixtures): for
        // every case, the memo's canonical classes are EXACTLY the
        // brute-force bit-equality classes (the exact predicate,
        // never the hash), and the memoized two-phase scan — the
        // scan site's own shape: phase A computes each distinct
        // canonical pair sum once, phase B fills every pair's own
        // flat slot — reproduces the direct per-pair accumulation
        // BIT FOR BIT at every pair, clones and near-twins included.
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for case in 0..60u32 {
            let elsewhere = -10.0 - (next() % 200) as f64 / 7.0;
            let n_units = 12 + (next() as usize) % 40;
            let base_folds = 3 + (next() as usize) % 8;
            let extra_folds = (next() as usize) % 6;
            let n_folds = base_folds + extra_folds;
            let mut dense: Vec<Vec<f64>> = Vec::new();
            let mut pinned: Vec<Vec<bool>> = Vec::new();
            let mut mismatches: Vec<Vec<u32>> = Vec::new();
            for fold in 0..base_folds {
                let mut row = vec![elsewhere; n_units];
                let mut pin = vec![false; n_units];
                let mut mism = vec![0u32; n_units];
                for _ in 0..1 + (next() as usize) % 9 {
                    let unit = (next() as usize) % n_units;
                    if !pin[unit] {
                        pin[unit] = true;
                        mism[unit] = next() % 5;
                        // every so often a pinned entry carries
                        // EXACTLY E — invisible to the column
                        row[unit] = if next() % 7 == 0 {
                            elsewhere
                        } else {
                            elsewhere
                                + (next() % 97) as f64 / 13.0
                                + (fold as f64 + 0.5) / 4.0
                        };
                    }
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(mism);
            }
            for _ in 0..extra_folds {
                let source = (next() as usize) % base_folds;
                let mut row = dense[source].clone();
                let mut pin = pinned[source].clone();
                match next() % 3 {
                    // the exact clone (a bit-equal row)
                    0 => {}
                    // the cross-encoding twin: an extra pinned entry
                    // at exactly E — invisible, still bit-equal
                    1 => {
                        for unit in 0..n_units {
                            if !pin[unit] {
                                pin[unit] = true;
                                row[unit] = elsewhere;
                                break;
                            }
                        }
                    }
                    // the near-twin NEGATIVE CONTROL: one ulp off at
                    // a pinned non-E entry — a different column
                    _ => {
                        let mut moved = false;
                        for unit in 0..n_units {
                            if pin[unit] && row[unit].to_bits() != elsewhere.to_bits() {
                                row[unit] = f64::from_bits(row[unit].to_bits() + 1);
                                moved = true;
                                break;
                            }
                        }
                        if !moved {
                            // no non-E pinned entry to perturb: pin
                            // one fresh unit one ulp above E instead
                            for unit in 0..n_units {
                                if !pin[unit] {
                                    pin[unit] = true;
                                    row[unit] =
                                        f64::from_bits(elsewhere.to_bits() + 1);
                                    break;
                                }
                            }
                        }
                    }
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(vec![0u32; n_units]);
            }
            let matrix =
                ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
            let counts: Vec<f64> = (0..n_units)
                .map(|unit| 1.0 + (next() % 9) as f64 * (unit as f64 + 1.0))
                .collect();
            let pair_count = n_folds * (n_folds + 1) / 2;
            // the canonical classes are EXACTLY the brute-force
            // bit-equality classes (the exact predicate decides)
            let memo = ClassScanMemo::build(&matrix);
            for a in 0..n_folds {
                for b in 0..n_folds {
                    assert_eq!(
                        matrix.ll_bits_equal(a, b),
                        memo.canon[a] == memo.canon[b],
                        "case {case}: the canon classes must be exactly \
                         the bit-equality classes at ({a}, {b})"
                    );
                }
            }
            for fold in 0..n_folds {
                assert!(
                    memo.canon[fold] <= fold,
                    "case {case}: the representative is the class's min"
                );
            }
            // the direct reference
            let direct: Vec<f64> = (0..pair_count)
                .map(|flat| {
                    let (first, second) = unrank_pair(flat);
                    matrix.class_pair_sum(first, second, &counts)
                })
                .collect();
            // the memoized two-phase scan (the scan site's own shape)
            let mut values = vec![f64::NAN; pair_count];
            let distinct: Vec<usize> = (0..pair_count)
                .filter(|&flat| {
                    let (first, second) = unrank_pair(flat);
                    memo.canonical_flat(first, second) == flat
                })
                .collect();
            assert_eq!(
                distinct.len() as u64,
                memo.canonical_pairs,
                "case {case}: the canonical pair count"
            );
            for flat in &distinct {
                let (first, second) = unrank_pair(*flat);
                values[*flat] = matrix.class_pair_sum(first, second, &counts);
            }
            for flat in 0..pair_count {
                let (first, second) = unrank_pair(flat);
                let memoized = values[memo.canonical_flat(first, second)];
                assert_eq!(
                    memoized.to_bits(),
                    direct[flat].to_bits(),
                    "case {case}: the memoized pair {flat} must be the \
                     direct sum bit for bit"
                );
            }
        }
    }


    #[test]
    fn receipt_f64_text_is_the_non_roundtrip_serialization_class() {
        // THE RECEIPT-SERIALIZATION CLASS, NAMED AND PROOFED (the
        // frontier-1 probe's verify flags resolved to it, SS10.6.21):
        // this repo's serde_json carries NO float_roundtrip feature,
        // so a serialized f64 can emit JSON text that parses back to
        // a DIFFERENT f64 — one ulp off on the monotone bit axis.
        // The committed receipts' float fields (the instrument
        // receipt's class_log_likelihoods, the ingredients' ll_matrix)
        // are that emission: an offline re-accumulation of the
        // emitted text can sit 1-2 ulps from the emitted class sum at
        // pairs whose accumulation crosses affected values, while the
        // emitting BINARY reproduces its own receipt bit-identically
        // (the rerun receipt: same binary, same inputs, same bytes).
        // The class therefore lives in the EMISSION, never in the
        // scan's float path and never in a consumer's parse: the
        // byte-identity gates compare BYTES (unaffected); only
        // offline value-level recomputations see the ulps. The
        // counterexample below is pinned against the repo's locked
        // serde_json; if a future dependency flip enables
        // float_roundtrip, this test names exactly what changed.
        let value = 1.1919406916781266e-278f64;
        let text = serde_json::to_string(&value).unwrap();
        let parsed: f64 = serde_json::from_str(&text).unwrap();
        assert_eq!(text, "1.1919406916781266e-278");
        assert_ne!(
            parsed.to_bits(),
            value.to_bits(),
            "the emitted text must parse back one ulp off (the class)"
        );
        assert_eq!(ulp_distance(parsed, value), 1);
        // The class is COMMON, not a corner: a deterministic sweep
        // finds it across the exponent range (receipted ~31.5% of a
        // random sweep at the probe; here a small pinned sample).
        let mut bad = 0usize;
        let mut state = 0x243F_6A88_85A3_08D3u64;
        for _ in 0..100_000u64 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = f64::from_bits(state & !(0xFu64 << 60));
            if !v.is_finite() {
                continue;
            }
            let t = serde_json::to_string(&v).unwrap();
            let back: f64 = serde_json::from_str(&t).unwrap();
            if back.to_bits() != v.to_bits() {
                bad += 1;
            }
        }
        assert!(
            bad > 10_000,
            "the non-roundtrip class is the common case at this              dependency configuration (saw {bad} of 100k)"
        );
    }

    #[test]
    fn excess_mass_bound_dominates_every_class_pair_sum() {
        // THE PER-FOLD EXCESS-MASS BOUND'S OWN PROOF, at the scan
        // level: A_f >= 0 (a pinned entry's marginalization rounds to
        // >= E, never below — the row-compression exactness
        // argument), and E*C + A_f + A_g >= the exact class pair sum
        // for EVERY pair (the mix0 <= alpha + beta lemma the
        // committed design-proof test proves) — the admissibility
        // contract the bounded scan's every skip rests on.
        let elsewhere = -14.418156f64;
        let shapes = [
            (5usize, 40usize, 8usize),
            (5, 40, 2),
            (6, 40, 1),
            (4, 40, 0),
            (7, 33, 3),
        ];
        for (n_folds, n_units, pin_stride) in shapes {
            let mut dense: Vec<Vec<f64>> = Vec::new();
            let mut pinned: Vec<Vec<bool>> = Vec::new();
            let mut mismatches: Vec<Vec<u32>> = Vec::new();
            for fold in 0..n_folds {
                let mut row = vec![elsewhere; n_units];
                let mut pin = vec![false; n_units];
                let mut mism = vec![0u32; n_units];
                if pin_stride > 0 {
                    for unit in (fold..n_units).step_by(pin_stride) {
                        pin[unit] = true;
                        row[unit] = elsewhere + (unit as f64) / (fold as f64 + 1.0);
                        mism[unit] = (unit % 7) as u32;
                    }
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(mism);
            }
            let matrix =
                ScoreMatrix::from_dense(dense, pinned.clone(), mismatches, elsewhere);
            let counts: Vec<f64> = (0..n_units).map(|u| 1.0 + u as f64).collect();
            let index = ClassScanBoundIndex::build(&matrix, &counts);
            for &excess in &index.excess {
                assert!(
                    excess >= 0.0,
                    "A_f must be nonnegative (the pinned entries round to >= E)"
                );
            }
            for first in 0..n_folds {
                for second in first..n_folds {
                    let bound = index.pair_bound(first, second);
                    let exact = matrix.class_pair_sum(first, second, &counts);
                    assert!(
                        bound >= exact,
                        "bound {bound} below exact {exact} at ({first},{second})"
                    );
                    // THE SLACK-CORRECTED DOMINANCE (the shape class
                    // the committed five shapes do NOT include — the
                    // bound's f64 evaluation can round below the
                    // exact unit-order sum, SS10.6.22's receipted
                    // 12-ulp shortfall): exact <= bound + the
                    // derived forward-error term, at every shape.
                    let slack = index.pair_slack(first, second);
                    assert!(
                        exact <= bound + slack,
                        "exact {exact} above bound {bound} + derived slack {slack} \
                         at ({first},{second})"
                    );
                    assert!(
                        slack.is_finite() && slack >= 0.0,
                        "the derived slack must be a finite nonnegative term"
                    );
                }
            }
        }
    }

    #[test]
    fn bounded_class_scan_preserves_the_exhaustive_answers() {
        // THE SANCTIONED BOUNDED CLASS SCAN'S OWN ANSWER-PRESERVATION
        // PROOF, at the scan level, over shapes that include a
        // DEEP landscape (pairs thousands of nats below best — the
        // pruning regime) and the flat controls (nothing prunable):
        // (a) stage one's best is the exhaustive best bit for bit;
        // (b) the called set (the computed ties, flat order) is the
        // exhaustive called set; (c) every computed exact sum is the
        // exhaustive scan's own value bit for bit (the committed
        // unit-order accumulation); (d) every SKIPPED pair's exact
        // sum sits below best + the derived exp-underflow point, so
        // its every QUAL posterior term is exactly +0.0 — the bit
        // identity the answers rest on; (e) the winner is never
        // skipped (it is in the computed set — the bound's own
        // arithmetic); (f) stage one is a pure function of the
        // matrix (two runs bit-identical — the determinism gate).
        let elsewhere = -14.418156f64;
        // (fold-giant, sparse, flat, and a DEEP shape: fold 0 pinned
        // everywhere with a 3-nat excess and heavy counts — the
        // (0,0) pair sits thousands of nats above the all-elsewhere
        // pairs, past the 745-nat window, so stage two genuinely
        // prunes.)
        let shapes = [(5usize, 40usize, 8usize), (5, 40, 2), (4, 40, 0)];
        for (n_folds, n_units, pin_stride) in shapes {
            let mut dense: Vec<Vec<f64>> = Vec::new();
            let mut pinned: Vec<Vec<bool>> = Vec::new();
            let mut mismatches: Vec<Vec<u32>> = Vec::new();
            for fold in 0..n_folds {
                let mut row = vec![elsewhere; n_units];
                let mut pin = vec![false; n_units];
                let mut mism = vec![0u32; n_units];
                if pin_stride > 0 {
                    for unit in (fold..n_units).step_by(pin_stride) {
                        pin[unit] = true;
                        row[unit] = elsewhere + (unit as f64) / (fold as f64 + 1.0);
                        mism[unit] = (unit % 7) as u32;
                    }
                }
                dense.push(row);
                pinned.push(pin);
                mismatches.push(mism);
            }
            let matrix =
                ScoreMatrix::from_dense(dense, pinned, mismatches, elsewhere);
            let counts: Vec<f64> = (0..n_units).map(|u| 25.0 + 5.0 * u as f64).collect();
            let pair_count = n_folds * (n_folds + 1) / 2;
            // the exhaustive reference
            let exhaustive: Vec<f64> = (0..pair_count)
                .map(|flat| {
                    let (first, second) = unrank_pair(flat);
                    matrix.class_pair_sum(first, second, &counts)
                })
                .collect();
            let exhaustive_best = exhaustive.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let exhaustive_called: Vec<(usize, usize)> = (0..pair_count)
                .filter(|&flat| exhaustive[flat] == exhaustive_best)
                .map(unrank_pair)
                .collect();
            // stage one, twice (the determinism receipt), once with
            // the bit-equal-row memo and once without — the memo's
            // values are the direct form's own bits either way
            let (bounds_a, slacks_a, computed_a, best_a, _) =
                bounded_stage_one(&matrix, &counts, None);
            let (_bounds_b, slacks_b, computed_b, best_b, _) =
                bounded_stage_one(&matrix, &counts, None);
            let memo = ClassScanMemo::build(&matrix);
            let (_bounds_m, _slacks_m, computed_m, best_m, _) =
                bounded_stage_one(&matrix, &counts, Some(&memo));
            assert_eq!(best_a.to_bits(), best_b.to_bits());
            assert_eq!(slacks_a, slacks_b, "the derived slack is deterministic");
            assert_eq!(best_a.to_bits(), exhaustive_best.to_bits(), "stage one's best");
            // THE MEMOIZED STAGE ONE's own gate: the memo's kept
            // values are the direct form's own bits (the fourth
            // reduction's construction proof, exercised here at the
            // bounded path's stage one)
            assert_eq!(
                best_m.to_bits(),
                best_a.to_bits(),
                "the memoized stage one's best"
            );
            assert!(
                computed_a.len() == computed_m.len()
                    && computed_a
                        .iter()
                        .zip(&computed_m)
                        .all(|(a, b)| match (a, b) {
                            (Some(x), Some(y)) => x.to_bits() == y.to_bits(),
                            (None, None) => true,
                            _ => false,
                        }),
                "the memoized stage one's kept values are the direct form's own"
            );
            assert!(
                computed_a.len() == computed_b.len()
                    && computed_a
                        .iter()
                        .zip(&computed_b)
                        .all(|(a, b)| match (a, b) {
                            (Some(x), Some(y)) => x.to_bits() == y.to_bits(),
                            (None, None) => true,
                            _ => false,
                        }),
                "stage one is a pure function of the matrix"
            );
            let called: Vec<(usize, usize)> = computed_a
                .iter()
                .enumerate()
                .filter_map(|(flat, value)| {
                    let sum = (*value)?;
                    (sum == best_a).then(|| unrank_pair(flat))
                })
                .collect();
            assert_eq!(called, exhaustive_called, "the bounded called set");
            // stage two: the derived keep threshold (no truth pair at
            // this locus — the production call's own form)
            let keep_threshold = bounded_keep_threshold(best_a, None);
            assert_eq!(
                keep_threshold.to_bits(),
                (best_a + exp_underflow_point()).to_bits()
            );
            let mut known = computed_a.clone();
            for flat in 0..pair_count {
                if known[flat].is_none() && bounds_a[flat] + slacks_a[flat] >= keep_threshold {
                    let (first, second) = unrank_pair(flat);
                    known[flat] = Some(matrix.class_pair_sum(first, second, &counts));
                }
            }
            let skipped_any = known.iter().any(|value| value.is_none());
            for flat in 0..pair_count {
                match known[flat] {
                    Some(sum) => {
                        // every computed value is the exhaustive one
                        assert_eq!(
                            sum.to_bits(),
                            exhaustive[flat].to_bits(),
                            "kept pair {flat} bit-identical"
                        );
                    }
                    None => {
                        // every skipped pair is below the QUAL window:
                        // its exact sum underflows to +0.0 at best
                        assert!(
                            exhaustive[flat] < keep_threshold,
                            "skipped pair {flat} inside the QUAL window"
                        );
                        assert_eq!((exhaustive[flat] - best_a).exp(), 0.0);
                        assert!(
                            bounds_a[flat] + slacks_a[flat] < keep_threshold,
                            "the skip certificate's own slack-corrected threshold"
                        );
                    }
                }
            }
            // the deep shape actually prunes (the regime is reached)
            if pin_stride == 8 {
                assert!(
                    skipped_any,
                    "the deep shape must prune (the 745-nat window is exceeded)"
                );
            }
            // the winner is never skipped: the exhaustive argmax is
            // in the computed set
            let (winner_first, winner_second) = exhaustive_called[0];
            let winner_flat = winner_second * (winner_second + 1) / 2 + winner_first;
            assert!(computed_a[winner_flat].is_some(), "the winner is kept");
            // THE FINAL-VALUE TIE-SET REBUILD (the T03226 locus-16
            // fix's own construction): best and the called set as
            // pure functions of the stage-two-patched class_lls —
            // the exhaustive ones, at every shape above.
            let final_lls: Vec<f64> = known
                .iter()
                .map(|value| value.unwrap_or(f64::NEG_INFINITY))
                .collect();
            let (rebuilt_best, rebuilt_called) = called_set_from_class_lls(&final_lls);
            assert_eq!(
                rebuilt_best.to_bits(),
                exhaustive_best.to_bits(),
                "the rebuilt best"
            );
            assert_eq!(rebuilt_called, exhaustive_called, "the rebuilt called set");
        }
    }

    #[test]
    fn bounded_tie_set_rebuilt_from_the_final_class_lls() {
        // THE T03226 LOCUS-16 BUG CLASS, SYNTHESIZED AND FIXED: the
        // bound's f64 evaluation can sit strictly BELOW the exact
        // unit-order sum (the receipted 12-ulp shortfall — the bound
        // inadmissible-by-rounding at tie shapes), so stage one's
        // strict stop could fire legitimately BETWEEN bit-identical
        // ties and the emitted called set — built from stage one's
        // computed prefix only — silently lost them. The fix under
        // test, both halves: (a) THE DERIVED SLACK — every skip
        // comparison consumes bound + the bound's own
        // forward-error term, restoring exact <= bound + slack so
        // the stop can no longer fire between ties; (b) THE
        // FINAL-VALUE REBUILD — the emitted called set is a pure
        // function of the final class_lls, never of the bound's
        // rounding. SYNTHESIZING THE SHAPE CLASS: two BIT-IDENTICAL
        // folds (a guaranteed three-way tie: pairs (0,0), (0,1),
        // (1,1) share one accumulation bit for bit) whose pinned
        // entries sit a FEW ULPs ABOVE E (the committed
        // row-compression class — a pinned entry's marginalization
        // rounds to >= E — with the excess mass tiny), so the bound
        // fl(fl(E*C + A_f) + A_g) carries the big E*C-against-A
        // cancellation while the exact sums accumulate unit-order:
        // the bound can round strictly BELOW the exact tie —
        // precisely the receipted mechanism. The search is
        // deterministic (a fixed grid, no RNG) and asserts its own
        // premise found.
        let ulps_up = |value: f64, k: usize| {
            let mut bits = value.to_bits();
            for _ in 0..k {
                bits += 1; // value < 0: +1 bit pattern is one ulp toward +inf
            }
            f64::from_bits(bits)
        };
        let shape_rows = |n_units: usize, elsewhere: f64| -> Vec<Vec<f64>> {
            vec![
                (0..n_units)
                    .map(|unit| {
                        if unit % 5 < 3 {
                            ulps_up(elsewhere, 1 + unit % 3)
                        } else {
                            elsewhere
                        }
                    })
                    .collect(),
                (0..n_units)
                    .map(|unit| {
                        if unit % 5 < 3 {
                            ulps_up(elsewhere, 1 + unit % 3)
                        } else {
                            elsewhere
                        }
                    })
                    .collect(),
            ]
        };
        let shape_matrix = |n_units: usize, elsewhere: f64| -> ScoreMatrix {
            ScoreMatrix::from_dense(
                shape_rows(n_units, elsewhere),
                vec![
                    (0..n_units).map(|unit| unit % 5 < 3).collect(),
                    (0..n_units).map(|unit| unit % 5 < 3).collect(),
                ],
                vec![vec![0u32; n_units], vec![0u32; n_units]],
                elsewhere,
            )
        };
        let mut found: Option<(usize, f64, f64, f64)> = None;
        for n_units in [61usize, 131, 293] {
            let counts: Vec<f64> =
                (0..n_units).map(|unit| 1.0 + (unit % 9) as f64).collect();
            for step in 1..=400i32 {
                let elsewhere = -(2.0 + step as f64 * 0.0137713137);
                let matrix = shape_matrix(n_units, elsewhere);
                let exact0 = matrix.class_pair_sum(0, 0, &counts);
                let index = ClassScanBoundIndex::build(&matrix, &counts);
                let bound = index.pair_bound(0, 0);
                // the rounding class: the exact tie sits strictly
                // ABOVE the bound (the bound is
                // inadmissible-by-rounding at this shape), and all
                // three pairs share the bound (the tie shape)
                if exact0 > bound
                    && index.pair_bound(0, 1) == bound
                    && index.pair_bound(1, 1) == bound
                {
                    found = Some((n_units, elsewhere, bound, exact0));
                    break;
                }
            }
            if found.is_some() {
                break;
            }
        }
        let Some((n_units, elsewhere, bound, exact0)) = found else {
            panic!(
                "the search found no inadmissible-by-rounding tie shape — \
                 the test's premise shape class is gone"
            );
        };
        let counts: Vec<f64> =
            (0..n_units).map(|unit| 1.0 + (unit % 9) as f64).collect();
        let matrix = shape_matrix(n_units, elsewhere);
        let pair_count = 3usize;
        // all three pairs tie bit-identical at the accumulated sum
        let exhaustive: Vec<f64> = (0..pair_count)
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                matrix.class_pair_sum(first, second, &counts)
            })
            .collect();
        assert_eq!(exhaustive[0].to_bits(), exact0.to_bits());
        assert_eq!(exhaustive[1].to_bits(), exact0.to_bits());
        assert_eq!(exhaustive[2].to_bits(), exact0.to_bits());
        let (exhaustive_best, exhaustive_called) = called_set_from_class_lls(&exhaustive);
        assert_eq!(exhaustive_best.to_bits(), exact0.to_bits());
        assert_eq!(exhaustive_called, vec![(0, 0), (0, 1), (1, 1)]);
        // THE RECEIPTED MECHANISM, at this shape: the bound sits
        // strictly below the exact best (the old stop rule fired on
        // exactly this comparison after computing (0,0) alone)
        assert!(
            bound < exact0,
            "the premise: bound {bound} must sit below the exact tie {exact0}"
        );
        // (a) THE DERIVED SLACK: the bound-plus-slack dominance
        // holds at exactly the shape where the bare bound fails
        let index = ClassScanBoundIndex::build(&matrix, &counts);
        for flat in 0..pair_count {
            let (first, second) = unrank_pair(flat);
            let slack = index.pair_slack(first, second);
            assert!(slack.is_finite() && slack < 1e-6);
            assert!(
                exhaustive[flat] <= index.pair_bound(first, second) + slack,
                "the slack-corrected dominance must hold at pair {flat}"
            );
        }
        // and the slack-corrected stop no longer fires between the
        // ties: stage one computes the full tie set
        let (bounds, slacks, computed, best, _) =
            bounded_stage_one(&matrix, &counts, None);
        assert_eq!(best.to_bits(), exact0.to_bits(), "stage one's best");
        for flat in 0..pair_count {
            assert!(
                computed[flat].is_some(),
                "the slack-corrected stop must not skip tie pair {flat}"
            );
            assert_eq!(
                computed[flat].unwrap().to_bits(),
                exact0.to_bits(),
                "stage one's tie value must be the exhaustive one"
            );
        }
        // (b) THE FINAL-VALUE REBUILD, at the receipted old state:
        // the OLD code's stage one stopped after (0,0) — the two
        // other ties uncomputed, their equal bounds strictly below
        // the exact best — so the old class_lls carried the prefix
        // only, and the old construction (from the prefix) emitted
        // the short set: the bug, receipted
        let old_class_lls = vec![exact0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        let (_, old_called) = called_set_from_class_lls(&old_class_lls);
        assert_eq!(old_called, vec![(0, 0)], "the old prefix-only emission");
        // stage two's keep window (the production form: no truth
        // pair) selects the two skipped ties of the old state —
        // their bounds + slack sit far above the derived threshold
        let keep_threshold = bounded_keep_threshold(best, None);
        let mut patched = old_class_lls;
        for flat in 0..pair_count {
            if patched[flat] == f64::NEG_INFINITY
                && bounds[flat] + slacks[flat] >= keep_threshold
            {
                let (first, second) = unrank_pair(flat);
                patched[flat] = matrix.class_pair_sum(first, second, &counts);
            }
        }
        assert_eq!(patched[1].to_bits(), exact0.to_bits());
        assert_eq!(patched[2].to_bits(), exact0.to_bits());
        // and the rebuild from the FINAL values emits the full tie
        let (rebuilt_best, rebuilt_called) = called_set_from_class_lls(&patched);
        assert_eq!(rebuilt_best.to_bits(), exact0.to_bits());
        assert_eq!(rebuilt_called, exhaustive_called);
    }

    #[test]
    fn bounded_keep_threshold_carries_the_truth_window() {
        // THE TRUTH WINDOW'S OWN SOUNDNESS: where the truth pair sits
        // BELOW the QUAL window (best + t), the derived keep
        // threshold drops to the truth pair's own exact ll, so every
        // pair that can be strictly above or tied with the truth is
        // kept (a skipped pair's bound < truth_ll forces its exact <
        // truth_ll — neither above nor tied) and the receipt's truth
        // fields read the exhaustive values; where the truth pair
        // sits INSIDE the window the threshold is the QUAL window's
        // own (the receipted shape — the log gaps are under the
        // 745-nat window).
        let best = -1234.5f64;
        let t = exp_underflow_point();
        assert_eq!(
            bounded_keep_threshold(best, None).to_bits(),
            (best + t).to_bits()
        );
        assert_eq!(
            bounded_keep_threshold(best, Some(best)).to_bits(),
            (best + t).to_bits(),
            "the truth pair at best adds nothing"
        );
        let deep_truth = best + 2.0 * t; // ~1490 nats below best
        assert_eq!(
            bounded_keep_threshold(best, Some(deep_truth)).to_bits(),
            deep_truth.to_bits(),
            "a truth pair below the window drops the threshold to itself"
        );
        assert!(deep_truth < best + t);
    }

    // ---------------------------------------------------------------------
    // (THE PANEL-SIDE PER-LOCUS PRECOMPUTE, SS10.6.23 slice 1 — the
    // unit bars the design names: the asset round trip is lossless;
    // the serve == live-built fold layouts on randomized fixtures
    // INCLUDING the masked forms (the masking commute lemma); the
    // representative-independence clause; the sequence-FNV drift
    // detector is a loud error, never a served structure.)
    // ---------------------------------------------------------------------

    /// A tiny deterministic LCG (no tuned constants, reproducible).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn folds_equal(a: &[Fold], b: &[Fold]) {
        assert_eq!(
            a.len(),
            b.len(),
            "fold count differs: live {} vs served {}",
            a.len(),
            b.len()
        );
        for (index, (live, served)) in a.iter().zip(b.iter()).enumerate() {
            assert_eq!(
                live.seq, served.seq,
                "fold {index}: the sequence differs"
            );
            assert_eq!(live.len, served.len, "fold {index}: the length differs");
            assert_eq!(live.walk, served.walk, "fold {index}: the walk differs");
            assert_eq!(
                live.node_positions, served.node_positions,
                "fold {index}: the node positions differ"
            );
            assert_eq!(live.nodes, served.nodes, "fold {index}: the nodes differ");
            assert_eq!(live.edges, served.edges, "fold {index}: the edges differ");
            assert_eq!(
                live.members, served.members,
                "fold {index}: the members differ"
            );
        }
    }

    /// The randomized-fixture world: a synthetic row set over a few
    /// paths, with (a) rows sharing their (sequence, relative walk)
    /// content but not their extents (the coalesced-copy class, the
    /// representative-independence clause), (b) unique rows, (c) rows
    /// whose member entries repeat. The fetch closure keys on
    /// (path name, start, end) exactly as the live fetch does.
    struct Fixture {
        names: Vec<String>,
        path_of_name: HashMap<String, usize>,
        rows: Vec<RowFold>,
        dedup: HashMap<(usize, u64, u64), usize>,
        members: Vec<MemberRow>,
        axis_start: u64,
        axis_end: u64,
        component: String,
    }

    fn precompute_fixture(seed: u64) -> Fixture {
        let mut rng = Lcg(seed);
        let names: Vec<String> = (0..5)
            .map(|i| format!("DONOR#{i}#chr1"))
            .chain(std::iter::once("STRANGER#0#chr1".to_string()))
            .collect();
        let mut path_of_name = HashMap::new();
        for (index, name) in names.iter().enumerate() {
            path_of_name.insert(name.clone(), index);
        }
        let component = "AXIS#0#chr1".to_string();
        let axis_path = names.len();
        path_of_name.insert(component.clone(), axis_path);
        let stranger_path = names.len() - 1;
        let mut rows: Vec<RowFold> = Vec::new();
        let mut members: Vec<MemberRow> = Vec::new();
        // The CONTENT-CLONE class: identical sequence bytes and
        // identical contained relative walks, at different extents and
        // paths — one fold whichever row represents it.
        let clone_len = 80u64;
        let clone_seq: Vec<u8> = (0..clone_len)
            .map(|i| b"ACGT"[(i * 7 + (i / 13)) as usize % 4])
            .collect();
        let clone_walk: Vec<(u64, i32)> = (0..6)
            .map(|i| (20 + i * 9, if i % 3 == 2 { -(10 + i as i32) } else { 10 + i as i32 }))
            .collect();
        for clone in 0..3 {
            let start = 500 + clone * 137;
            let end = start + clone_len;
            let walk = clone_walk
                .iter()
                .map(|&(bp, node)| (bp + start, node))
                .collect::<Vec<_>>();
            let path = (clone as usize) % stranger_path;
            rows.push(RowFold {
                path,
                start,
                end,
                seq: clone_seq.clone(),
                walk: walk.clone(),
                frames: walk.iter().map(|&(_, node)| (node < 0) as u8).collect(),
                stored_walk: walk.clone(),
                members: Vec::new(),
            });
            members.push(MemberRow {
                path_name: names[path].clone(),
                start,
                end,
            });
        }
        // Unique rows with randomized walks.
        for row_index in 0..12 {
            let start = 40 + row_index as u64 * 90;
            let len = 60 + rng.below(60);
            let end = start + len;
            let mut walk: Vec<(u64, i32)> = Vec::new();
            let mut bp = start + rng.below(6);
            while bp + K <= end {
                let node = 1 + rng.below(30) as i32;
                walk.push((bp, if rng.below(2) == 0 { node } else { -node }));
                bp += 1 + rng.below(9);
            }
            walk.sort_unstable_by_key(|&(bp, _)| bp);
            let path = rng.below(stranger_path as u64) as usize;
            rows.push(RowFold {
                path,
                start,
                end,
                seq: (0..len)
                    .map(|i| b"ACGT"[((i as u64 * 11 + row_index as u64) % 4) as usize])
                    .collect(),
                walk: walk.clone(),
                frames: walk.iter().map(|&(_, node)| (node < 0) as u8).collect(),
                stored_walk: walk.clone(),
                members: Vec::new(),
            });
            members.push(MemberRow {
                path_name: names[path].clone(),
                start,
                end,
            });
            // A repeated member entry over the SAME row identity (the
            // deduplication class: one row, two member entries).
            if rng.below(3) == 0 {
                members.push(MemberRow {
                    path_name: names[path].clone(),
                    start,
                    end,
                });
            }
        }
        // The axis row itself (never excluded, the coordinate frame).
        let axis_start = 100u64;
        let axis_end = 900u64;
        let axis_walk: Vec<(u64, i32)> = (0..10)
            .map(|i| (axis_start + i * 80, 100 + i as i32))
            .collect();
        rows.push(RowFold {
            path: axis_path,
            start: axis_start,
            end: axis_end,
            seq: (0..axis_end - axis_start)
                .map(|i| b"ACGT"[(i % 4) as usize])
                .collect(),
            walk: axis_walk.clone(),
            frames: axis_walk.iter().map(|_| 0u8).collect(),
            stored_walk: axis_walk,
            members: Vec::new(),
        });
        members.push(MemberRow {
            path_name: component.clone(),
            start: axis_start,
            end: axis_end,
        });
        // The stranger member (the masking operator's target — a REAL
        // panel path masked at call time, so its row identity exists).
        let stranger_start = 5230u64;
        let stranger_end = stranger_start + 64;
        rows.push(RowFold {
            path: stranger_path,
            start: stranger_start,
            end: stranger_end,
            seq: vec![b'T'; (stranger_end - stranger_start) as usize],
            walk: vec![(stranger_start + 9, -77)],
            frames: vec![1],
            stored_walk: vec![(stranger_start + 9, -77)],
            members: Vec::new(),
        });
        members.push(MemberRow {
            path_name: names[stranger_path].clone(),
            start: stranger_start,
            end: stranger_end,
        });
        // The row table: dedup by (path, start, end), the members
        // attached in member order (the live insertion discipline).
        let mut dedup: HashMap<(usize, u64, u64), usize> = HashMap::new();
        for (index, row) in rows.iter().enumerate() {
            dedup.insert((row.path, row.start, row.end), index);
        }
        for (member_index, member) in members.iter().enumerate() {
            let path_idx = path_of_name[&member.path_name];
            let row_index = dedup[&(path_idx, member.start, member.end)];
            rows[row_index].members.push((0, member_index));
        }
        Fixture {
            names,
            path_of_name,
            rows,
            dedup,
            members,
            axis_start,
            axis_end,
            component,
        }
    }



    fn fixture_fetch(fixture: &Fixture) -> impl Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync {
        let rows = fixture.rows.clone();
        let path_of_name = fixture.path_of_name.clone();
        move |name: &str, start: u64, end: u64| {
            let path_idx = path_of_name[name];
            for row in &rows {
                if row.path == path_idx && row.start == start && row.end == end {
                    return Ok(row.seq.clone());
                }
            }
            Err(invalid("the fixture fetch misses a row"))
        }
    }

    fn build_serve_state(
        fixture: &Fixture,
    ) -> (PrecomputeServe, Vec<RowFold>, Vec<u64>) {
        // The builder's own path, in miniature: the row context asset
        // (with the FNV content keys), the unmasked live grouping, the
        // asset round trip through the REAL encoder/decoder.
        let rows: Vec<PrecomputeRowAsset> = fixture
            .rows
            .iter()
            .map(|row| PrecomputeRowAsset {
                path_name: if let Some(name) = fixture.names.get(row.path) {
                    name.clone()
                } else {
                    fixture.component.clone()
                },
                start: row.start,
                end: row.end,
                seq_fnv: fnv1a64(&row.seq),
                stored_walk: row.stored_walk.clone(),
                walk: row.walk.clone(),
                frames: row.frames.clone(),
                members: row.members.clone(),
            })
            .collect();
        let encoded = encode_rows(&rows);
        let decoded = decode_rows(&encoded).expect("the asset round trip decodes");
        // The served rows, built as run() builds them (sequences
        // left in the sources).
        let served_rows: Vec<RowFold> = decoded
            .iter()
            .map(|row| RowFold {
                path: fixture.path_of_name[&row.path_name],
                start: row.start,
                end: row.end,
                seq: Vec::new(),
                walk: row.walk.clone(),
                frames: row.frames.clone(),
                stored_walk: row.stored_walk.clone(),
                members: row.members.clone(),
            })
            .collect();
        let (_, fold_rep, member_fold, _, _) = build_window_folds_live(
            &fixture.rows,
            &fixture.dedup,
            &fixture.members,
            &fixture.component,
            fixture.axis_start,
            fixture.axis_end,
            None,
            &[],
            &fixture.path_of_name,
            &|_| u64::MAX / 4,
            K,
        )
        .expect("the unmasked live grouping builds");
        let windows = vec![PrecomputeWindowAsset {
            partition: 0,
            axis_start: fixture.axis_start,
            axis_end: fixture.axis_end,
            fold_rep,
            member_fold,
        }];
        let serve = PrecomputeServe {
            windows,
            seq_fnv: decoded.iter().map(|row| row.seq_fnv).collect(),
            seqs: Mutex::new(HashMap::new()),
            bind: None,
        };
        (serve, served_rows, decoded.iter().map(|row| row.seq_fnv).collect())
    }

    #[test]
    fn panel_precompute_asset_round_trip_is_lossless() {
        for seed in [1u64, 7, 42, 20261023] {
            let fixture = precompute_fixture(seed);
            let rows: Vec<PrecomputeRowAsset> = fixture
                .rows
                .iter()
                .map(|row| PrecomputeRowAsset {
                    path_name: fixture.names[row.path.min(fixture.names.len() - 1)].clone(),
                    start: row.start,
                    end: row.end,
                    seq_fnv: fnv1a64(&row.seq),
                    stored_walk: row.stored_walk.clone(),
                    walk: row.walk.clone(),
                    frames: row.frames.clone(),
                    members: row.members.clone(),
                })
                .collect();
            let decoded = decode_rows(&encode_rows(&rows)).expect("decode");
            assert_eq!(decoded.len(), rows.len());
            for (a, b) in rows.iter().zip(decoded.iter()) {
                assert_eq!(a.path_name, b.path_name);
                assert_eq!(a.start, b.start);
                assert_eq!(a.end, b.end);
                assert_eq!(a.seq_fnv, b.seq_fnv);
                assert_eq!(a.stored_walk, b.stored_walk);
                assert_eq!(a.walk, b.walk);
                assert_eq!(a.frames, b.frames);
                assert_eq!(a.members, b.members);
            }
        }
    }

    #[test]
    fn panel_precompute_serve_equals_live_built_on_randomized_fixtures() {
        // The masking commute lemma, exercised on every fixture and
        // every filter combination: unmasked, the flank filter, the
        // masking operator, and both together.
        for seed in [1u64, 7, 42, 20261023] {
            let fixture = precompute_fixture(seed);
            let (serve, served_rows, _) = build_serve_state(&fixture);
            let fetch = fixture_fetch(&fixture);
            let mut path_names = fixture.names.clone();
            path_names.push(fixture.component.clone());
            for (flank, exclude) in [
                (None::<FlankOverlap>, Vec::<String>::new()),
                (
                    Some(FlankOverlap::Bp(120)),
                    Vec::<String>::new(),
                ),
                (None::<FlankOverlap>, vec!["STRANGER#*".to_string()]),
                (
                    Some(FlankOverlap::Bp(120)),
                    vec!["STRANGER#*".to_string()],
                ),
            ] {
                let (live_folds, _rep, _member_fold, excluded, flank_excluded) =
                    build_window_folds_live(
                        &fixture.rows,
                        &fixture.dedup,
                        &fixture.members,
                        &fixture.component,
                        fixture.axis_start,
                        fixture.axis_end,
                        flank,
                        &exclude,
                        &fixture.path_of_name,
                        &|_| u64::MAX / 4,
                        K,
                    )
                    .expect("the live masked build");
                let (served_folds, s_excluded, s_flank_excluded) = serve
                    .serve_window_folds(
                        &served_rows,
                        0,
                        &fixture.members,
                        &fixture.component,
                        fixture.axis_start,
                        fixture.axis_end,
                        flank,
                        &exclude,
                        K,
                        &fetch,
                        &path_names,
                    )
                    .expect("the serve");
                assert_eq!(excluded, s_excluded, "the excluded counts differ");
                assert_eq!(
                    flank_excluded, s_flank_excluded,
                    "the flank-excluded counts differ"
                );
                folds_equal(&live_folds, &served_folds);
                // The representative-independence clause, live on the
                // fixture: the content-clone rows (equal sequence,
                // equal relative walk, different extents) group into
                // ONE fold, and whichever row represents it the
                // served content is identical to the live one —
                // asserted on the unmasked form (the masked forms may
                // filter the class away, which is the filters' own
                // semantics, proved identical above).
                if flank.is_none() && exclude.is_empty() {
                    let coalesced = live_folds
                        .iter()
                        .any(|fold| fold.seq == fixture.rows[0].seq && fold.members.len() >= 3);
                    assert!(
                        coalesced,
                        "the fixture must exercise the coalesced-copy class \
                         (seed {seed})"
                    );
                }
            }
        }
    }

    #[test]
    fn panel_precompute_serve_requires_the_axis_fold_under_territory_masks() {
        // The masking filters never exclude the axis row (it is the
        // coordinate frame), so the territory machinery's axis-fold
        // discovery holds under the serve exactly as under the live
        // build — asserted here over the fixtures' masked forms.
        for seed in [1u64, 7, 42] {
            let fixture = precompute_fixture(seed);
            let (serve, served_rows, _) = build_serve_state(&fixture);
            let fetch = fixture_fetch(&fixture);
            let mut path_names = fixture.names.clone();
            path_names.push(fixture.component.clone());
            let (folds, _, _) = serve
                .serve_window_folds(
                    &served_rows,
                    0,
                    &fixture.members,
                    &fixture.component,
                    fixture.axis_start,
                    fixture.axis_end,
                    Some(FlankOverlap::Bp(120)),
                    &["STRANGER#*".to_string()],
                    K,
                    &fetch,
                    &path_names,
                )
                .expect("the masked serve");
            assert!(folds.iter().any(|fold| {
                fold.members.iter().any(|m| {
                    m.path_name == fixture.component
                        && m.start == fixture.axis_start
                        && m.end == fixture.axis_end
                })
            }));
        }
    }

    #[test]
    fn panel_precompute_sequence_fnv_drift_is_a_loud_error() {
        // The drift detector: a source whose bytes no longer match the
        // asset's recorded FNV-64 content key is a HARD ERROR at serve
        // time — never a silently served structure.
        let fixture = precompute_fixture(7);
        let (serve, served_rows, mut seq_fnv) = build_serve_state(&fixture);
        // Corrupt one representative row's content key.
        let rep = serve.windows[0].fold_rep[0];
        seq_fnv[rep] ^= 1;
        let serve = PrecomputeServe {
            windows: serve.windows,
            seq_fnv,
            seqs: Mutex::new(HashMap::new()),
            bind: None,
        };
        let fetch = fixture_fetch(&fixture);
        let mut path_names = fixture.names.clone();
        path_names.push(fixture.component.clone());
        let result = serve.serve_window_folds(
            &served_rows,
            0,
            &fixture.members,
            &fixture.component,
            fixture.axis_start,
            fixture.axis_end,
            None,
            &[],
            K,
            &fetch,
            &path_names,
        );
        assert!(result.is_err(), "a drifted sequence content key must be loud");
        assert!(
            result
                .err()
                .unwrap()
                .to_string()
                .contains("mismatches the fetched bytes"),
            "the error must name the drift class"
        );
    }

    // ---------------------------------------------------------------------
    // (THE BIND SCAFFOLDING, SS10.6.25 slice 2 — the unit bars the
    // design names: the asset round trip is lossless; the serve ==
    // live derivation at EVERY position of randomized in-process
    // panels (the locality lemma — served band, uncovered fallback,
    // path edges included); the loader serves the slice-1 /1 format
    // without the section and validates the /2 section's counters.)
    // ---------------------------------------------------------------------

    #[test]
    fn bind_scaffold_asset_round_trip_is_lossless() {
        for seed in [1u64, 7, 42, 20261024] {
            let mut rng = Lcg(seed);
            let paths: Vec<BindScaffoldPathAsset> = (0..4)
                .map(|path| {
                    let ranges: Vec<(u64, u64)> = (0..1 + rng.below(3))
                        .map(|_| {
                            let lo = rng.below(2000);
                            (lo, lo + 1 + rng.below(500))
                        })
                        .collect();
                    let steps: Vec<(u64, i32)> = (0..rng.below(200))
                        .map(|_| {
                            let bp = rng.below(2000);
                            let node = rng.below(1 << 20) as i64;
                            (
                                bp,
                                ((node >> 1) as i64 ^ -(node & 1)) as i32,
                            )
                        })
                        .collect();
                    BindScaffoldPathAsset {
                        path_name: format!("PATH#{path}#chrT:0-1"),
                        ranges,
                        steps,
                    }
                })
                .collect();
            let decoded =
                decode_bind_scaffolding(&encode_bind_scaffolding(&paths)).expect("decode");
            assert_eq!(decoded.len(), paths.len());
            for (a, b) in paths.iter().zip(decoded.iter()) {
                assert_eq!(a.path_name, b.path_name);
                assert_eq!(a.ranges, b.ranges);
                assert_eq!(a.steps, b.steps);
            }
        }
    }

    #[test]
    fn bind_scaffold_serve_equals_live_at_every_position() {
        // The syng C library's global state serializes against every
        // other in-process test that touches the FFI (the shared
        // SYNG_TEST_LOCK — the in-process panels below are real
        // indexes built through the same C library).
        let _guard = crate::syng::SYNG_TEST_LOCK
            .lock()
            .unwrap_or_else(|e| e.into_inner());
        let k = 63u64; // the committed k63/s8 convention's syncmer length
        for seed in [3u64, 19, 91] {
            let mut rng = Lcg(seed);
            // The randomized-panel world: three paths of varied length
            // (path-edge bands exercised), deterministic ACGT with an
            // N-run on one path (the invalid-base skip class) and a
            // shared segment between two paths (the multi-match
            // dictionary class — the reverse map's node can come from
            // another path's dictionary entry).
            let lengths = [1400usize, 1900, 2300];
            let mut sequences: Vec<Vec<u8>> = lengths
                .iter()
                .map(|&len| {
                    (0..len)
                        .map(|_| b"ACGT"[rng.below(4) as usize])
                        .collect()
                })
                .collect();
            for i in 600..640 {
                sequences[1][i] = b'N';
            }
            let shared_segment: Vec<u8> = sequences[0][100..400].to_vec();
            for (i, &base) in shared_segment.iter().enumerate() {
                sequences[2][1500 + i] = base;
            }
            let names: Vec<String> =
                (0..3).map(|p| format!("PATH#{p}#chrT:0-{}", lengths[p])).collect();
            let panel = crate::syng::SyngIndex::build(
                crate::syng::SyncmerParams::default(),
                names.iter().cloned().zip(sequences.iter().cloned()),
            );
            let path_of_name: HashMap<String, usize> = names
                .iter()
                .enumerate()
                .map(|(path, name)| (name.clone(), path))
                .collect();
            // The in-memory sources, fetched with the production
            // fetch's own semantics (the hi clamp at the path length;
            // lo >= hi is empty).
            let fetch = |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
                let path = path_of_name[name];
                let len = lengths[path] as u64;
                let hi = hi.min(len);
                if lo >= hi {
                    return Ok(Vec::new());
                }
                Ok(sequences[path][lo as usize..hi as usize].to_vec())
            };
            // Randomized member-row extents per path: overlapping,
            // gapped (a gap wider than the window's reach must exist —
            // the uncovered fallback class), an extent at the path
            // start and one at the path end (the clamped-context
            // clause).
            let mut member_rows: Vec<(usize, u64, u64)> = Vec::new();
            for (path, &len) in lengths.iter().enumerate() {
                member_rows.push((path, 0, 100 + rng.below(150)));
                member_rows.push((path, len as u64 - 200, len as u64));
                let mut cursor = 300u64;
                while cursor + 100 < len as u64 - 400 {
                    let span = 100 + rng.below(400);
                    member_rows.push((path, cursor, cursor + span));
                    // The gap: sometimes within the window's reach
                    // (merged), sometimes beyond it (uncovered).
                    cursor += span
                        + if rng.below(2) == 0 {
                            rng.below(40)
                        } else {
                            100 + rng.below(200)
                        };
                }
            }
            let assets =
                build_bind_scaffold(&panel, &fetch, k, member_rows.iter().copied())
                    .expect("the scaffold builds");
            let serve = BindScaffoldServe::from_asset(assets.clone(), &path_of_name, &panel)
                .expect("the serve state builds");
            let mut served_positions = 0u64;
            let mut fallback_positions = 0u64;
            for (path, &len) in lengths.iter().enumerate() {
                let path_len = len as u64;
                for position in 0..path_len {
                    // THE LIVE DERIVATION (the committed bind's own
                    // spelling, at the record's own window).
                    let live = canonical_steps_over_range(
                        &panel,
                        &fetch,
                        k,
                        path,
                        position.saturating_sub(150),
                        position + 150,
                    )
                    .expect("the live derivation")
                    .into_iter()
                    .find(|&(bp, _)| bp == position);
                    // The band classification must be exactly the
                    // locality lemma's predicate over the builder's
                    // own ranges.
                    let covered = assets
                        .iter()
                        .find(|asset| asset.path_name == names[path])
                        .map(|asset| {
                            asset.ranges.iter().any(|&(lo, hi)| {
                                bind_range_covers(
                                    lo,
                                    hi,
                                    position,
                                    path_len,
                                    panel.params.w as u64 - 1,
                                    panel.syncmer_length_bp() as u64,
                                )
                            })
                        })
                        .unwrap_or(false);
                    match serve.step_at(path, position, path_len) {
                        Some(recorded) => {
                            assert!(covered, "a position outside the served band was served");
                            assert_eq!(
                                recorded, live,
                                "path {path} position {position}: the served step differs \
                                 from the live derivation"
                            );
                            served_positions += 1;
                        }
                        None => {
                            assert!(
                                !covered,
                                "path {path} position {position}: a covered position \
                                 fell back (the serve must answer it)"
                            );
                            fallback_positions += 1;
                        }
                    }
                }
            }
            // The fixture must exercise BOTH classes at every seed
            // (the served band and the uncovered fallback).
            assert!(
                served_positions > 0 && fallback_positions > 0,
                "seed {seed}: the fixture must exercise both the served band \
                 ({served_positions}) and the fallback ({fallback_positions})"
            );
            assert_eq!(
                serve.served.load(Ordering::Relaxed),
                served_positions,
                "the served counter must count the lookups"
            );
            assert_eq!(
                serve.fallbacks.load(Ordering::Relaxed),
                fallback_positions,
                "the fallback counter must count the uncovered positions"
            );
        }
    }

    #[test]
    fn bind_scaffold_loader_serves_slice1_format_and_validates_counters() {
        // The slice-1 /1 manifest (no bind counters, no section file)
        // loads with `bind: None` — the row context + fold layouts
        // serve, the bind derives live; the /2 format requires the
        // section and validates its counters against the manifest.
        let dir = std::env::temp_dir().join(format!(
            "impg-bind-scaffold-loader-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let identity = crate::genome_inference::PanelIdentity {
            checksum_algorithm: "fnv1a64".to_string(),
            sidecars: Vec::new(),
        };
        let keys = PanelPrecomputeKeys {
            panel: identity.clone(),
            k: 63,
            component: "PATH#1#chrT:0-1000".to_string(),
            frame: "local".to_string(),
            record_origin: 0,
            stored_frame: false,
            maps_digest: 1,
            routes_digest: 2,
        };
        let write_asset = |format: &str, bind: Option<&Vec<BindScaffoldPathAsset>>| {
            std::fs::create_dir_all(&dir).expect("temp dir");
            let (bind_paths, bind_ranges, bind_steps) = match bind {
                Some(paths) => (
                    paths.len() as u64,
                    paths.iter().map(|p| p.ranges.len() as u64).sum(),
                    paths.iter().map(|p| p.steps.len() as u64).sum(),
                ),
                None => (0, 0, 0),
            };
            let manifest = PanelPrecomputeManifest {
                format: format.to_string(),
                panel: identity.clone(),
                k: 63,
                component: keys.component.clone(),
                frame: "local".to_string(),
                record_origin: 0,
                stored_frame: false,
                maps_digest: 1,
                routes_digest: 2,
                windows: 0,
                rows: 0,
                folds: 0,
                bind_paths,
                bind_ranges,
                bind_steps,
            };
            crate::genome_inference::atomic_write(
                &dir.join("rows.bin"),
                &encode_rows(&[]),
            )
            .expect("rows");
            crate::genome_inference::atomic_write(
                &dir.join("windows.bin"),
                &encode_windows(&[]),
            )
            .expect("windows");
            if let Some(paths) = bind {
                crate::genome_inference::atomic_write(
                    &dir.join("bind-scaffolding.bin"),
                    &encode_bind_scaffolding(paths),
                )
                .expect("bind");
            }
            crate::genome_inference::write_json(&dir.join("manifest.json"), &manifest)
                .expect("manifest");
        };
        // The slice-1 /1 form: no section, `bind: None`.
        write_asset("impg-panel-precompute-loci/1", None);
        let asset = load_panel_precompute(&dir, &keys).expect("the /1 asset loads");
        assert!(asset.bind.is_none(), "the /1 asset carries no bind section");
        // The /2 form with the section.
        let paths = vec![BindScaffoldPathAsset {
            path_name: "PATH#1#chrT:0-1000".to_string(),
            ranges: vec![(0, 600)],
            steps: vec![(10, 7), (40, -9)],
        }];
        write_asset("impg-panel-precompute-loci/2", Some(&paths));
        let asset = load_panel_precompute(&dir, &keys).expect("the /2 asset loads");
        let bind = asset.bind.expect("the /2 asset carries the section");
        assert_eq!(bind.len(), 1);
        assert_eq!(bind[0].ranges, paths[0].ranges);
        assert_eq!(bind[0].steps, paths[0].steps);
        // A corrupted counter is a loud error, never a served section:
        // the manifest still carries the good counters while the
        // section file was replaced underneath it.
        let mut bad = paths.clone();
        bad[0].steps.clear();
        crate::genome_inference::atomic_write(
            &dir.join("bind-scaffolding.bin"),
            &encode_bind_scaffolding(&bad),
        )
        .expect("bind");
        let result = load_panel_precompute(&dir, &keys);
        assert!(result.is_err(), "a mismatched counter must be loud");
        assert!(
            result
                .err()
                .unwrap()
                .to_string()
                .contains("bind-scaffolding counters disagree"),
            "the error must name the counter class"
        );
        // An unknown format string is refused outright.
        write_asset("impg-panel-precompute-loci/9", None);
        let result = load_panel_precompute(&dir, &keys);
        assert!(result.is_err(), "an unknown format must be loud");
        assert!(
            result.err().unwrap().to_string().contains("unknown format"),
            "the error must name the format class"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// THE PACKED PREDICATE SET: set semantics, ascending iteration, and
    /// the cardinality counter (the receipt's entry count) that the
    /// packed form must keep identical to the list's length.
    #[test]
    fn cut_set_packs_members_and_iterates_them_ascending() {
        // insertion order must not matter, duplicates must not double
        let mut set = CutSet::default();
        for position in [130usize, 0, 65, 3, 64, 130, 129] {
            set.insert(position);
        }
        assert_eq!(set.len(), 6, "duplicates must not count");
        assert!(!set.is_empty());
        assert_eq!(set.iter().collect::<Vec<_>>(), vec![0, 3, 64, 65, 129, 130]);
        // empty stays empty, and the counter stays zero
        let empty = CutSet::default();
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
        assert_eq!(empty.iter().count(), 0);
        // from_ascending agrees with the incremental build
        let built = CutSet::from_ascending([0usize, 3, 64, 65, 129, 130]);
        assert_eq!(built, set);
        // the packed form is 1 bit per position, not 8 bytes per member
        let wide = CutSet::from_ascending(0..1000usize);
        assert_eq!(wide.len(), 1000);
        assert_eq!(wide.iter().collect::<Vec<_>>(), (0..1000).collect::<Vec<_>>());
        // THE FOOTPRINT IS THE POINT: ~1 bit per spectrum position, not
        // 8 bytes per member (the unpacked form's constant), so a dense
        // predicate set is ~64x smaller and the degenerate monster
        // (|called| sets over a pair_count-long spectrum) stops being
        // an allocation that cannot exist.
        let packed_bytes = wide.words.len() * std::mem::size_of::<u64>();
        let unpacked_bytes = wide.len() * std::mem::size_of::<usize>();
        assert!(
            packed_bytes * 8 <= unpacked_bytes,
            "packed {packed_bytes} B must be at least 8x under unpacked {unpacked_bytes} B"
        );
    }

    /// THE SPLICED RECEIPT WRITER IS THE VALUE-TREE FORM, BYTE FOR
    /// BYTE. The certificate streams from its typed lists instead of
    /// being built as a `serde_json::Value` first (the 32-bytes-per-entry
    /// copy the receipts-phase audit found); the emission order is
    /// serde_json's ordered map, so the two forms must agree exactly -
    /// otherwise the receipts stop being comparable to the committed
    /// ones, which is the whole point of the doctrine.
    #[test]
    fn the_spliced_report_writer_matches_the_value_tree_form_byte_for_byte() {
        // ASCENDING, as the producer builds them (it walks the spectrum
        // in position order); insertion-order independence is CutSet's
        // own test.
        let lists: Vec<(usize, Vec<usize>)> =
            vec![(0, vec![1, 3, 4]), (7, Vec::new()), (2, vec![0, 9])];
        // The sentinel keeps the key in the map (its alphabetical slot),
        // and the reference builds the member the committed way.
        let report = json!({
            "alpha": 1,
            "zeta": "x",
            "mid": [1.5, 2],
            "qual_tied_cut_indices": serde_json::Value::Null,
        });
        let mut reference_value = report.clone();
        reference_value["qual_tied_cut_indices"] = json!(lists
            .iter()
            .map(|(index, list)| json!([index, list]))
            .collect::<Vec<_>>());
        let mut reference = Vec::new();
        serde_json::to_writer(&mut reference, &reference_value).unwrap();

        let packed: Vec<(usize, CutSet)> = lists
            .iter()
            .map(|(index, list)| (*index, CutSet::from_ascending(list.iter().copied())))
            .collect();
        let mut spliced = Vec::new();
        write_report_with_cut_indices(
            &mut spliced,
            &report,
            "qual_tied_cut_indices",
            &packed,
        )
        .unwrap();

        assert_eq!(
            String::from_utf8(spliced).unwrap(),
            String::from_utf8(reference).unwrap(),
            "the streamed certificate must emit the value tree's own bytes"
        );
    }
}
