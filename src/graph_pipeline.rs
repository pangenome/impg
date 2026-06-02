//! Typed parser for graph-construction pipeline specs.
//!
//! The CLI grammar is intentionally small:
//!
//! ```text
//! stage[,key=value...][:stage[,key=value...]...
//! ```
//!
//! The parser validates syntax and normalizes stage/parameter names, but does
//! not decide whether a stage is executable. Runtime dispatch remains owned by
//! the command layer.

use std::collections::HashSet;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GraphPipelineSpec {
    pub stages: Vec<GraphPipelineStage>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GraphPipelineStage {
    pub name: String,
    pub params: Vec<GraphPipelineParam>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GraphPipelineParam {
    pub key: String,
    pub value: String,
}

impl GraphPipelineSpec {
    pub fn parse(raw: &str) -> Result<Self, String> {
        let raw = raw.trim();
        if raw.is_empty() {
            return Err("empty graph pipeline spec".to_string());
        }

        let mut stages = Vec::new();
        for (stage_idx, raw_stage) in raw.split(':').enumerate() {
            let raw_stage = raw_stage.trim();
            if raw_stage.is_empty() {
                return Err(format!("empty stage at position {}", stage_idx + 1));
            }
            stages.push(GraphPipelineStage::parse(raw_stage, stage_idx)?);
        }

        Ok(Self { stages })
    }

    pub fn to_spec(&self) -> String {
        self.stages
            .iter()
            .map(GraphPipelineStage::to_spec)
            .collect::<Vec<_>>()
            .join(":")
    }

    pub fn stages_from(&self, start: usize) -> Self {
        Self {
            stages: self.stages[start..].to_vec(),
        }
    }
}

impl GraphPipelineStage {
    fn parse(raw: &str, stage_idx: usize) -> Result<Self, String> {
        let mut pieces = raw.split(',').map(str::trim);
        let name_raw = pieces.next().unwrap_or_default();
        let name = normalize_name(name_raw);
        if name.is_empty() {
            return Err(format!("empty stage name at position {}", stage_idx + 1));
        }

        let mut params = Vec::new();
        let mut seen = HashSet::new();
        for piece in pieces {
            if piece.is_empty() {
                return Err(format!("stage '{}' has an empty parameter", name));
            }
            let (key_raw, value_raw) = piece.split_once('=').ok_or_else(|| {
                format!("stage '{}' parameter '{}' must be key=value", name, piece)
            })?;
            let key = normalize_name(key_raw);
            let value = value_raw.trim().to_string();
            if key.is_empty() {
                return Err(format!("stage '{}' has an empty parameter key", name));
            }
            if value.is_empty() {
                return Err(format!(
                    "stage '{}' parameter '{}' has empty value",
                    name, key
                ));
            }
            if !seen.insert(key.clone()) {
                return Err(format!("stage '{}' repeats parameter '{}'", name, key));
            }
            params.push(GraphPipelineParam { key, value });
        }

        Ok(Self { name, params })
    }

    pub fn param(&self, key: &str) -> Option<&str> {
        let key = normalize_name(key);
        self.params
            .iter()
            .find(|param| param.key == key)
            .map(|param| param.value.as_str())
    }

    pub fn to_spec(&self) -> String {
        let mut out = self.name.clone();
        for param in &self.params {
            out.push(',');
            out.push_str(&param.key);
            out.push('=');
            out.push_str(&param.value);
        }
        out
    }
}

fn normalize_name(raw: &str) -> String {
    raw.trim().replace('_', "-").to_ascii_lowercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_stage_params() {
        let spec =
            GraphPipelineSpec::parse("syng,k=63,s=8,seed=7:crush,max-span=10k,max-traversals=128")
                .unwrap();
        assert_eq!(spec.stages.len(), 2);
        assert_eq!(spec.stages[0].name, "syng");
        assert_eq!(spec.stages[0].param("k"), Some("63"));
        assert_eq!(spec.stages[1].name, "crush");
        assert_eq!(spec.stages[1].param("max_span"), Some("10k"));
        assert_eq!(
            spec.to_spec(),
            "syng,k=63,s=8,seed=7:crush,max-span=10k,max-traversals=128"
        );
    }

    #[test]
    fn preserves_legacy_pipeline_shapes() {
        let spec = GraphPipelineSpec::parse("wfmash:seqwish:10k").unwrap();
        assert_eq!(
            spec.stages
                .iter()
                .map(|s| s.name.as_str())
                .collect::<Vec<_>>(),
            vec!["wfmash", "seqwish", "10k"]
        );

        let spec = GraphPipelineSpec::parse("sweepga:fastga:pggb,window=20k").unwrap();
        assert_eq!(
            spec.stages
                .iter()
                .map(|s| s.name.as_str())
                .collect::<Vec<_>>(),
            vec!["sweepga", "fastga", "pggb"]
        );
        assert_eq!(spec.stages[2].param("window"), Some("20k"));
    }

    #[test]
    fn rejects_bad_params() {
        assert!(GraphPipelineSpec::parse("").is_err());
        assert!(GraphPipelineSpec::parse("syng:").is_err());
        assert!(GraphPipelineSpec::parse("syng,k").is_err());
        assert!(GraphPipelineSpec::parse("syng,k=63,").is_err());
        assert!(GraphPipelineSpec::parse("syng,k=63,k=31").is_err());
    }
}
