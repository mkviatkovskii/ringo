pub mod canonical;
pub mod reader;
pub mod writer;

use crate::io::smiles::reader::molecule::DirectionalStereo;
use crate::io::smiles::writer::DoubleBondStereoPolicy;

/// Combined stereo tolerance settings for reading and writing SMILES.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TolerateStereo {
    /// Policy for invalid, incomplete, or unpreserved `/` and `\` marks on input.
    pub directional: DirectionalStereo,
    /// Policy for unrepresentable double-bond configurations on output.
    pub double_bond: DoubleBondStereoPolicy,
}

impl TolerateStereo {
    /// Reject all invalid or unrepresentable stereo (default).
    pub const STRICT: Self = Self {
        directional: DirectionalStereo::Reject,
        double_bond: DoubleBondStereoPolicy::Reject,
    };

    /// Tolerate all supported stereo issues on both input and output.
    pub const TOLERATE: Self = Self {
        directional: DirectionalStereo::Tolerate,
        double_bond: DoubleBondStereoPolicy::Tolerate,
    };

    /// Applies the umbrella `tolerate_all` override on top of individual policies.
    #[must_use]
    pub const fn from_flags(
        tolerate_all: bool,
        directional: DirectionalStereo,
        double_bond: DoubleBondStereoPolicy,
    ) -> Self {
        if tolerate_all {
            Self::TOLERATE
        } else {
            Self {
                directional,
                double_bond,
            }
        }
    }
}
