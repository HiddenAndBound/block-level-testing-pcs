use binius_field::{BinaryField128b, BinaryField32b};
use crate::utils::merkle::Hash;

pub mod merkle;
pub mod ntt;
pub mod packed_arithmetic;
pub mod encoding;

#[derive(Debug)]
pub struct Commitment{
    pub commit: Hash,
    pub cols: usize
}

pub struct EvalProof{
    pub folded_poly: Vec<BinaryField128b>,
    pub queried_columns: Vec<Vec<BinaryField32b>>,
    pub merkle_paths: Vec<Vec<Hash>>
}


impl EvalProof{
    pub fn new(
        folded_poly: Vec<BinaryField128b>,
        queried_columns: Vec<Vec<BinaryField32b>>,
        merkle_paths: Vec<Vec<Hash>>
    )->EvalProof{
        EvalProof{
            folded_poly,
            queried_columns,
            merkle_paths,
        }
    }
}