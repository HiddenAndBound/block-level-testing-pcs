use std::time::Instant;

use binius_field::{BinaryField128b, BinaryField1b, BinaryField32b, Field};
use rand::{distributions::Uniform, thread_rng, Rng};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{prover::{commit, prove}, utils::{encoding::compute_fourier_bases, ntt::AdditiveNTT, packed_arithmetic::evaluate_unpacked}, verifier::verify};

#[test]
fn commitment_test(){

    for l in 10..40{

        println!("--------------|| length 2^{:?} ||-------------- \n\n", l+5);
        let time = Instant::now();
        let poly:Vec<BinaryField32b> = (0..1<<l).into_par_iter().map(|_| {
            BinaryField32b::random(thread_rng())
        }).collect();
        println!("Time: {:?} \n", time.elapsed());

        println!("Making ntt precomputations");
        let time = Instant::now();
        let ntt = AdditiveNTT::new(l);
        println!("Time: {:?} \n", time.elapsed());

        println!("Committing");

        let time = Instant::now();
        let (commitment, merkle_tree, encoded_poly) = commit(&poly, &ntt);
        println!("Time: {:?} \n", time.elapsed());

        let point: Vec<BinaryField128b> = (0..l+5).into_iter().map(|_|BinaryField1b::random(thread_rng()).into()).collect();

        let base = compute_fourier_bases(&point);
        // println!("{:?}", base);
        let eval = evaluate_unpacked(&poly, &base);
        let mut rng = thread_rng();
        println!("Generating queries");
        let queries = (0..241).into_iter().map(|_| thread_rng().sample(Uniform::new(0, encoded_poly.cols))).collect();
        println!("Generating proof");
        let time = Instant::now();
        let eval_proof = prove(&poly, &encoded_poly, &merkle_tree, &point, &queries);
        println!("Time: {:?}\n", time.elapsed());

        let time = Instant::now();

        println!("Verifying");
        verify(commitment, eval, eval_proof, point, queries, &ntt);
        println!("Time: {:?} \n", time.elapsed());
    }
}