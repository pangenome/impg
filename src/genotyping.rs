//! Shared vocabulary for backend-neutral genotyping and inference.
//!
//! The current implementation is syng-backed, but the genotyping pipeline is
//! intended to support multiple graph backends, evidence projections, and
//! scoring methods. These enums keep user-visible metadata stable while giving
//! future backends a common language.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureSpace {
    SyngSyncmerNode,
    VariationGraphNode,
    LocalHaplotypeSegment,
    HaplotypeSequence,
    MemHit,
}

impl FeatureSpace {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SyngSyncmerNode => "syng-syncmer-node",
            Self::VariationGraphNode => "variation-graph-node",
            Self::LocalHaplotypeSegment => "local-haplotype-segment",
            Self::HaplotypeSequence => "haplotype-sequence",
            Self::MemHit => "mem-hit",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceBackend {
    Pack,
    Projection,
    SyncmerWalk,
    VariationGraphAlignment,
    HaplotypeAlignment,
    MemProjection,
}

impl EvidenceBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pack => "pack",
            Self::Projection => "projection",
            Self::SyncmerWalk => "syncmer-walk",
            Self::VariationGraphAlignment => "variation-graph-alignment",
            Self::HaplotypeAlignment => "haplotype-alignment",
            Self::MemProjection => "mem-projection",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringMethod {
    Cos,
    CountLikelihood,
    AlignmentLikelihood,
    LearnedLocalEmission,
    GbwtImputation,
}

impl ScoringMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Cos => "cos",
            Self::CountLikelihood => "count-likelihood",
            Self::AlignmentLikelihood => "alignment-likelihood",
            Self::LearnedLocalEmission => "learned-local-emission",
            Self::GbwtImputation => "gbwt-imputation",
        }
    }

    pub fn metric_str(self) -> &'static str {
        match self {
            Self::Cos => "cosine",
            Self::CountLikelihood => "count-likelihood",
            Self::AlignmentLikelihood => "alignment-likelihood",
            Self::LearnedLocalEmission => "learned-local-emission",
            Self::GbwtImputation => "gbwt-imputation",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_user_visible_names() {
        assert_eq!(FeatureSpace::SyngSyncmerNode.as_str(), "syng-syncmer-node");
        assert_eq!(EvidenceBackend::Pack.as_str(), "pack");
        assert_eq!(ScoringMethod::Cos.as_str(), "cos");
        assert_eq!(ScoringMethod::Cos.metric_str(), "cosine");
    }
}
