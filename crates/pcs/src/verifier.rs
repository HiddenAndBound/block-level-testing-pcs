
use crate::utils::{encoding::{compute_fourier_bases, encode_extension}, merkle::{hash_field, verify_merkle_path}, ntt::AdditiveNTT, packed_arithmetic::PackedAlgebra32, Commitment, EvalProof};
use binius_field::{BinaryField128b, BinaryField32b, ExtensionField, Field};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};

pub fn verify(
    commit: Commitment,
    eval: BinaryField128b,
    proof:EvalProof,
    point: Vec<BinaryField128b>,
    queries: Vec<usize>,
    ntt: &AdditiveNTT
){

    let scalars = compute_fourier_bases(&point[..point.len() - (commit.cols + 5)].to_vec());
    let encoded_poly = encode_extension(&proof.folded_poly, ntt);

    for i in 0..queries.len(){
        let leaf_hash = hash_field(&proof.queried_columns[i]);

        verify_merkle_path(&commit.commit, leaf_hash, queries[i],&proof.merkle_paths[i]);
        let linear_combination =  unpacked_linear_combination(&scalars, &proof.queried_columns[i]);
        assert_eq!(encoded_poly[queries[i]],linear_combination, "Query {i} failed")
    }

    let scalars = compute_fourier_bases(&point[point.len() - (commit.cols + 5)..].to_vec());

    assert_eq!(eval, linear_combination(&scalars, &proof.folded_poly));

}


pub fn linear_combination<F0:ExtensionField<F1>, F1:Field>(scalars: &Vec<F0>, vals: &Vec<F1>)->F0{

    vals.par_iter().zip(scalars.par_iter()).map(|(val, scalar)|  *scalar* *val).reduce(||F0::ZERO, |acc,g| acc + g)
}

pub fn unpacked_linear_combination(scalars: &Vec<BinaryField128b>, vals: &Vec<BinaryField32b>)->PackedAlgebra32{

    let mut res = [BinaryField128b::ZERO; 32];

    res.par_iter_mut().enumerate().for_each(|(i, val)|
        for k in 0..scalars.len(){
            if (vals[k].val()>>i)&1==1{
                *val += scalars[k]
            }
        }
    );
    PackedAlgebra32::new(res)
}