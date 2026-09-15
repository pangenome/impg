//! The same standalone diagnostic, with portable public replay and analytic controls.
#[path = "../examples/panel_route_convex_diagnostic/mod.rs"]
mod diagnostic;
use diagnostic::math::*;
use impg::{
    genome_inference::{self as genome, joint, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::json;
use std::fs;
#[test]
fn analytic_fractional_gap_derivative_delta_and_pricing_scope() {
    let c = vec![1.0];
    let columns = vec![vec![0.0], vec![2.0]];
    let scores = columns.iter().map(|s| loss(&c, s)).collect::<Vec<_>>();
    let r = optimize(&c, &columns, &scores, true).unwrap();
    let last = r.trace.last().unwrap();
    close(last.signal[0], 0.9).unwrap();
    close(last.objective, loss(&c, &[0.9])).unwrap();
    assert!(r.discrete_minimum - last.objective > 0.3);
    assert!(!r.rigorous_real_certificate);
    for row in &r.trace {
        assert!(row.pricing.numerical_tangent_lower_bound.unwrap() <= loss(&c, &[0.9]) + 1e-12);
    }
    let restricted = optimize(&c, &[vec![0.0]], &[0.0], false).unwrap();
    assert!(restricted.trace[0].objective > loss(&c, &[0.9]));
    assert!(restricted.correlated_optimum_column_indices.is_empty());
    for row in &restricted.trace {
        assert!(
            row.pricing.numerical_tangent_lower_bound.is_none() && row.pricing.raw_gap.is_none()
        );
    }
    let coupled_columns = vec![vec![0.0, 0.0], vec![0.9, 0.9]];
    let coupled_c = vec![1.0, 1.0];
    let coupled_scores = coupled_columns
        .iter()
        .map(|s| loss(&coupled_c, s))
        .collect::<Vec<_>>();
    let coupled = optimize(&coupled_c, &coupled_columns, &coupled_scores, true).unwrap();
    assert!(coupled_scores[1] < coupled_scores[0]); // analytic two-move domain excludes single moves
    let overshoot = vec![100.0];
    assert!(dot(&gradient(&c, &[0.0]), &overshoot) < 0.0);
    assert!(delta(&c, &[0.0], &overshoot) > 0.0);
    close(
        delta(&c, &[0.0], &overshoot),
        loss(&c, &overshoot) - loss(&c, &[0.0]),
    )
    .unwrap();
    if let Some(path) = std::env::var_os("IMPG_TEST_CONVEX_ANALYTIC_OUTPUT") {
        let path = std::path::PathBuf::from(path);
        fs::create_dir(&path).unwrap();
        diagnostic::write(&path.join("analytic.json"),&json!({
            "scope":"analytic vectors, not generated biological profiles",
            "fractional":r,"restricted_incomplete_pricing":restricted,
            "coupled_analytic":coupled,"analytic_minimum":loss(&c,&[0.9]),
            "overshoot":{"c":c,"from":[0.0],"to":overshoot,"derivative":dot(&gradient(&c,&[0.0]),&overshoot),"delta":delta(&c,&[0.0],&overshoot)}
        })).unwrap();
    }
    let c = vec![0.0, 3.0, 1.0];
    let s = vec![2.0, 0.0, 0.5];
    let t = vec![1.0, 0.0, 1.5];
    assert_eq!(gradient(&c, &s)[0], 1.0);
    close(delta(&c, &s, &t), loss(&c, &t) - loss(&c, &s)).unwrap();
    assert_eq!(loss(&[3.0], &[0.0]), 0.0); // unsupported positive is background-relative zero, not masked
    let h = 1e-6;
    close(
        gradient(&[1.0], &[0.5])[0],
        (loss(&[1.0], &[0.5 + h]) - loss(&[1.0], &[0.5 - h])) / (2.0 * h),
    )
    .unwrap();
}
#[test]
fn diagnostic_caps_and_invalid_values_fail_closed() {
    assert!(optimize(
        &[1.0],
        &vec![vec![0.0]; MAX_COLUMNS + 1],
        &vec![0.0; MAX_COLUMNS + 1],
        true
    )
    .is_err());
    assert!(optimize(&[1.0], &[vec![-1.0]], &[0.0], true).is_err());
    assert!(optimize(&[f64::NAN], &[vec![0.0]], &[0.0], true).is_err());
}
#[test]
fn public_l150_counts_reverse_duplicate_description_and_physical_copy() {
    let temp = tempfile::tempdir().unwrap();
    let retained =
        std::env::var_os("IMPG_TEST_CONVEX_PORTABLE_OUTPUT").map(std::path::PathBuf::from);
    if let Some(root) = &retained {
        fs::create_dir(root).expect("fresh portable public fixture directory");
    }
    let root = retained.as_deref().unwrap_or_else(|| temp.path());
    let mut seed = 113u64;
    let dna: Vec<u8> = (0..220)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize]
        })
        .collect();
    let sequences = vec![
        ("A#0#one".to_string(), dna.clone()),
        ("B#0#two".to_string(), dna.clone()),
    ];
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    let prefix = root.join("p.syng");
    panel.save(prefix.to_str().unwrap()).unwrap();
    let identity = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let fasta = root.join("sources.fa");
    fs::write(
        &fasta,
        sequences
            .iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let alternatives:Vec<_>=[(0,"+"),(0,"-"),(0,"+"),(1,"+")].iter().enumerate().map(|(i,(source,strand))|json!({"id":format!("a{i}"),"topology":"linear","left_endpoint":"asserted-molecule-terminus","right_endpoint":"asserted-molecule-terminus","pieces":[{"instance":format!("i{source}"),"start":0,"end":220,"strand":strand}],"adjacencies":[]})).collect();
    let layout=serde_json::from_value(json!({"version":1,"model":joint::LAYOUT_MODEL,"panel":identity,"sources":sequences.iter().enumerate().map(|(i,(name,s))|json!({"id":i,"name":name,"length":s.len()})).collect::<Vec<_>>(),"instances":(0..2).map(|i|json!({"id":format!("i{i}"),"source":i,"start":0,"end":220})).collect::<Vec<_>>(),"slots":[{"id":"molecule","alternatives":alternatives}]})).unwrap();
    let compiled = joint::compile(
        &panel,
        identity.clone(),
        layout,
        &[fasta.to_str().unwrap().into()],
        vec![150],
        None,
    )
    .unwrap();
    let reads = root.join("reads.fa");
    fs::write(
        &reads,
        dna.windows(150)
            .enumerate()
            .map(|(i, s)| format!(">r{i}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let sample = sample::build(&panel, identity, &[reads]).unwrap();
    let out = root.join("diagnostic");
    fs::create_dir(&out).unwrap();
    let summary = diagnostic::analyze(&compiled, &sample, None, &out).unwrap();
    assert_eq!(summary["feasible_columns"], 4);
    assert_eq!(
        summary["exact_count_equivalent_column_groups"][0],
        json!([0, 1, 2, 3])
    );
    assert_eq!(summary["correlated_optimum_columns"], json!([0, 1, 2, 3]));
    let rows: serde_json::Value = genome::read_json(&out.join("hypotheses.json")).unwrap();
    assert_eq!(rows[0]["physical_routes"], rows[2]["physical_routes"]); // duplicate description, not a new copy
    assert_ne!(rows[0]["physical_routes"], rows[3]["physical_routes"]); // genuinely different source capacity
    assert_eq!(
        rows[0]["counts_by_feature_length"],
        rows[3]["counts_by_feature_length"]
    );
    let features: serde_json::Value = genome::read_json(&out.join("features.json")).unwrap();
    for (i, t) in features["tokens"].as_array().unwrap().iter().enumerate() {
        let tokens: Vec<u64> = serde_json::from_value(t.clone()).unwrap();
        assert_eq!(
            rows[0]["counts_by_feature_length"][i][0].as_u64().unwrap(),
            sample.counts.count(&tokens).unwrap()
        );
    }
    let mut pair_layout = compiled.layout.clone();
    let mut first = pair_layout.slots[0].clone();
    first.alternatives = vec![first.alternatives[0].clone()];
    let mut second = pair_layout.slots[0].clone();
    second.id = "second-molecule".into();
    second.alternatives = vec![second.alternatives[3].clone()];
    pair_layout.slots = vec![first, second];
    let pair = joint::compile(
        &panel,
        compiled.layout.panel.clone(),
        pair_layout.clone(),
        &[fasta.to_str().unwrap().into()],
        vec![150],
        None,
    )
    .unwrap();
    let pair_out = root.join("two-physical-copies");
    fs::create_dir(&pair_out).unwrap();
    diagnostic::analyze(&pair, &sample, None, &pair_out).unwrap();
    let pair_rows: serde_json::Value =
        genome::read_json(&pair_out.join("hypotheses.json")).unwrap();
    for (one, two) in rows[0]["counts_by_feature_length"]
        .as_array()
        .unwrap()
        .iter()
        .zip(pair_rows[0]["counts_by_feature_length"].as_array().unwrap())
    {
        assert_eq!(two[0].as_u64().unwrap(), 2 * one[0].as_u64().unwrap());
    }
    // Explicit distinct instance IDs can assert extra copies at identical source
    // coordinates. They are NOT equivalent to route canonical-span capacity.
    pair_layout.instances[1].source = 0;
    let incompatible = joint::compile(
        &panel,
        compiled.layout.panel.clone(),
        pair_layout,
        &[fasta.to_str().unwrap().into()],
        vec![150],
        None,
    )
    .unwrap();
    let error = diagnostic::analyze(&incompatible, &sample, None, &pair_out).unwrap_err();
    assert!(error.to_string().contains("feasibility differ"));
    diagnostic::write(
        &root.join("expected-capacity-mismatch.json"),
        &json!({"expected_rejection":error.to_string(),"bounds_withheld":true}),
    )
    .unwrap();
    // Selecting both physical source paths doubles integer multiplicity, while
    // selecting forward/RC views of the same path violates canonical capacity.
    use diagnostic::canonical_feasible as resources_feasible;
    use impg::genome_inference::panel_routes::Segment;
    let one = Segment {
        source: 0,
        start: 0,
        end: 220,
        reverse: false,
    };
    let mut reverse = one.clone();
    reverse.reverse = true;
    assert!(!resources_feasible([one.clone(), reverse]));
    let mut second = one.clone();
    second.source = 1;
    assert!(resources_feasible([one, second]));
}
#[test]
#[ignore = "requires explicit fresh synthetic finite/coupled fixture roots and fresh output directory"]
fn configured_finite_and_coupled_public_pricing() {
    let base = std::path::PathBuf::from(
        std::env::var_os("IMPG_TEST_CONVEX_OUTPUT").expect("IMPG_TEST_CONVEX_OUTPUT required"),
    );
    fs::create_dir(&base).unwrap();
    for (name, var, products, columns) in [
        ("finite", "IMPG_TEST_CONVEX_FINITE", 513, 513),
        ("coupled", "IMPG_TEST_CONVEX_COUPLED", 4096, 64),
    ] {
        let root = std::path::PathBuf::from(
            std::env::var_os(var).expect("explicit fixture configuration required"),
        );
        let out = base.join(name);
        let summary = diagnostic::run_fixture(&root, &out).unwrap();
        assert_eq!(summary["products_examined"], products);
        assert_eq!(summary["feasible_columns"], columns);
        assert_eq!(summary["public_panel_evaluations"], columns);
        if name == "finite" {
            assert_eq!(summary["generated_zero"], true);
            assert_eq!(summary["unsupported_positive"], true);
        } else {
            assert_eq!(
                summary["correlated_optimum_columns"]
                    .as_array()
                    .unwrap()
                    .len(),
                64
            );
            let rows: serde_json::Value = genome::read_json(&out.join("hypotheses.json")).unwrap();
            let rows = rows.as_array().unwrap();
            assert!(rows.windows(2).all(|w| w[0]["counts_by_feature_length"]
                == w[1]["counts_by_feature_length"]
                && w[0]["physical_routes"] != w[1]["physical_routes"]));
            // Any changed feasible pair requires BOTH choices to change: the
            // independent single-move marginals contain canonical conflicts.
            for a in rows {
                for b in rows {
                    if a != b {
                        assert_ne!(a["choices"][0], b["choices"][0]);
                        assert_ne!(a["choices"][1], b["choices"][1]);
                    }
                }
            }
        }
    }
}
