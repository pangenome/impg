use super::*;
pub(super) fn policy_identity() -> BTreeMap<String, String> {
    [
        ("version", "fair-sample-ranked-source-pair-order-v3"),
        ("policy", include_str!("../../src/commands/genome_infer/panel_route_search_policy/mod.rs")),
        ("machine", include_str!("../../src/commands/genome_infer/panel_route_search_policy/machine.rs")),
        ("adapter", include_str!("../../src/commands/genome_infer/panel_route_search_policy/adapter.rs")),
        ("checkpoint", include_str!("../../src/commands/genome_infer/panel_route_search_policy/checkpoint.rs")),
        ("cli", include_str!("../../src/commands/genome_infer/panel_route_search_policy/../../genome_infer.rs")),
        ("continuation", include_str!("../../src/commands/genome_infer/panel_route_search_policy/continuation.rs")),
        ("transition", include_str!("../../src/commands/genome_infer/panel_route_search_policy/transition.rs")),
        ("source_order", include_str!("../../src/commands/genome_infer/panel_route_search_policy/source_order.rs")),
        ("v1_schema", include_str!("../../src/commands/genome_infer/panel_route_search_policy/v1_schema.rs")),
        ("validation", include_str!("../../src/commands/genome_infer/panel_route_search_policy/validation.rs")),
        (
            "validation_containers",
            include_str!("../../src/commands/genome_infer/panel_route_search_policy/validation_containers.rs"),
        ),
        ("supported_v1_policy", include_str!("../../src/commands/genome_infer/panel_route_search_policy/supported_v1/mod.rs")),
        (
            "supported_v1_machine",
            include_str!("../../src/commands/genome_infer/panel_route_search_policy/supported_v1/machine.rs"),
        ),
        (
            "supported_v1_adapter",
            include_str!("../../src/commands/genome_infer/panel_route_search_policy/supported_v1/adapter.rs"),
        ),
        (
            "supported_v1_checkpoint",
            include_str!("../../src/commands/genome_infer/panel_route_search_policy/supported_v1/checkpoint.rs"),
        ),
        (
            "supported_v1_cli",
            include_str!("../../src/commands/genome_infer/panel_route_search_policy/supported_v1/genome_infer.rs"),
        ),
    ]
    .into_iter()
    .map(|(name, source)| {
        (
            name.into(),
            if name == "version" {
                source.into()
            } else {
                format!(
                    "fnv1a64-source-{:016x}",
                    genome::checksum(source.as_bytes())
                )
            },
        )
    })
    .collect()
}
