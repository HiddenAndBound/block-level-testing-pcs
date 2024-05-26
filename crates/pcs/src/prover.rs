use binius_field::{BinaryField128b, BinaryField32b, Field};
use rand::thread_rng;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::utils::{encoding::{compute_fourier_bases, make_coeff_matrix, make_linear_combination, Code}, merkle::{hash_field, merklize, Hash, MerkleTree}, ntt::AdditiveNTT, packed_arithmetic::PackedAlgebra32, Commitment, EvalProof};

pub fn commit(poly:&Vec<BinaryField32b>, ntt:&AdditiveNTT)->(Commitment, MerkleTree, Code){
    let encoded_poly = Code::new(poly, ntt);
    let leaf_hashes:Vec<Hash> = (0..encoded_poly.cols).into_par_iter().map(|column| hash_field(&encoded_poly.col(column))).collect();
    let merkle_tree = merklize(leaf_hashes);
    let commitment = Commitment{
        commit: merkle_tree.get_root(),
        cols: encoded_poly.cols.trailing_zeros() as usize
    };

    (commitment, merkle_tree, encoded_poly)
}


pub fn prove(poly:&Vec<BinaryField32b>, encoded_poly: &Code, merkle_tree: &MerkleTree, point:&Vec<BinaryField128b>, queries:&Vec<usize>)->EvalProof{

    let scalars = compute_fourier_bases(&point[..(point.len() - (encoded_poly.cols.trailing_zeros() as usize + 5))].to_vec());
    let poly_coeff_matrix = make_coeff_matrix(&poly, encoded_poly.cols);

    let linear_combination = make_linear_combination(poly_coeff_matrix, scalars);

    // println!("{:?}", linear_combination);
    let mut queried_columns = vec![vec![BinaryField32b::ZERO; encoded_poly.rows]; queries.len()];
    let mut merkle_paths = vec![Vec::<Hash>::new(); queries.len()];

    for i in 0..queries.len(){
        queried_columns[i] = encoded_poly.col(queries[i]);
        merkle_paths[i] = merkle_tree.get_merkle_path(queries[i])
    }


    EvalProof::new(linear_combination, queried_columns, merkle_paths)
}




#[test]

fn prover_test(){
    
    let mut rng = thread_rng();
    let poly:Vec<BinaryField32b> = (0..1<<6).into_iter().map(|_| BinaryField32b::random(&mut rng)).collect();
    let ntt = AdditiveNTT::new(poly.len().trailing_zeros() as usize);

    let (commitment, merkle_tree, encoded_poly) = commit(&poly, &ntt);

    let point = vec![BinaryField128b::random(&mut rng);6];
    let queries = vec![2, 3];

    prove(&poly, &encoded_poly, &merkle_tree, &point, &queries);

}