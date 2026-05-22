use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::io;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let Some(gfa_path) = args.next() else {
        eprintln!("usage: cargo run --example povu_decomp_report -- <graph.gfa> [reference-name]");
        std::process::exit(2);
    };
    let refs: Vec<String> = args.collect();
    let gfa = fs::read_to_string(&gfa_path)?;
    let graph = povu::NativeGfa::parse(&gfa)
        .map_err(|err| io::Error::other(format!("POVU parse failed: {err}")))?;
    let decomposition = graph
        .decompose_flubbles(&refs)
        .map_err(|err| io::Error::other(format!("POVU decomposition failed: {err}")))?;

    let mut by_level: BTreeMap<usize, usize> = BTreeMap::new();
    let mut leaf_count = 0usize;
    for site in &decomposition.sites {
        *by_level.entry(site.level).or_default() += 1;
        if site.is_leaf {
            leaf_count += 1;
        }
    }

    println!("#gfa\t{gfa_path}");
    println!("#reference\t{}", decomposition.reference_path);
    println!(
        "#reference_path_index\t{}",
        decomposition.reference_path_index
    );
    println!("#segments\t{}", graph.segments.len());
    println!("#links\t{}", graph.links.len());
    println!("#paths\t{}", graph.paths.len());
    println!("#sites\t{}", decomposition.sites.len());
    println!("#leaf_sites\t{leaf_count}");
    for (level, count) in by_level {
        println!("#level\t{level}\t{count}");
    }
    println!(
        "id\tparent_id\tlevel\tis_leaf\tref_start_step\tref_end_step\tref_span_steps\tstart\tend"
    );

    let mut sites: Vec<_> = decomposition.sites.iter().collect();
    sites.sort_by_key(|site| {
        (
            std::cmp::Reverse(
                site.reference_end_step
                    .saturating_sub(site.reference_start_step),
            ),
            site.reference_start_step,
            site.reference_end_step,
        )
    });
    for site in sites {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            site.id,
            site.parent_id.as_deref().unwrap_or("."),
            site.level,
            site.is_leaf,
            site.reference_start_step,
            site.reference_end_step,
            site.reference_end_step
                .saturating_sub(site.reference_start_step),
            site.start.token(),
            site.end.token(),
        );
    }
    Ok(())
}
