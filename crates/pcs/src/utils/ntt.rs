use std::{ ops::{AddAssign, Mul}, process::Output, time::Instant};

use binius_field::{BinaryField, BinaryField1b, BinaryField32b, ExtensionField, Field, TowerField};
use rand::{random, thread_rng};
use rayon::{iter::{IntoParallelIterator, IntoParallelRefIterator,IntoParallelRefMutIterator, IndexedParallelIterator, ParallelIterator}, prelude::ParallelSliceMut};

pub struct  AdditiveNTT{
    log_transform_size: usize,
    twiddles: Vec<Vec<BinaryField32b>>
}

impl  AdditiveNTT {
    
    pub fn new(log_transform_size: usize)->AdditiveNTT{
        let twiddles = compute_twiddles(log_transform_size);
        AdditiveNTT{
            log_transform_size,
            twiddles
        }
    }

    //Forward ntt using precomputed twiddles, we dont parallelise here but rather will do so accross ntts when we encode our polynomial in the pcs.
    pub fn forward_ntt<F:Copy + Mul<BinaryField32b, Output =F> + AddAssign<<F as Mul<BinaryField32b>>::Output>>(
        &self,
        coeffs: &mut Vec<F>,
        coset:u32
    ){
        //Depth of the recursion, which is the base-2 logarithm of the length of the input.

        let rounds = coeffs.len().trailing_zeros();

        for r in (0..rounds).rev(){

            //Number of subproblems at this depth of the recursion.
            let parts = 1<<r;
            
            let normalising_value = vanishing_map(<BinaryField32b as TowerField>::basis(0, r as usize).unwrap(), r as usize).invert().unwrap();
            let coeset_twiddle = vanishing_map(BinaryField32b::new(coset), r as usize)*normalising_value;

            for p in 0..parts{
                for b in 0..(1<<(rounds-1 - r)){
                    
                    //the butterfly in this ntt indexes (b*2^(r+1) + p) and (b*2^(r+1) + p + 2^r)
                    let left_idx =(b<<(r+1))+ p;
                    let right_idx = left_idx + (1<<r);

                    let right_val = coeffs[right_idx];
                    coeffs[left_idx] += right_val * (self.twiddles[r as usize][b] + coeset_twiddle);

                    let left_val = coeffs[left_idx];
                    coeffs[right_idx] += left_val;
                }
            }


        }
    }

    pub fn inverse_ntt<F:Copy + Mul<BinaryField32b, Output =F> + AddAssign<<F as Mul<BinaryField32b>>::Output>>(
        &self,
        coeffs: &mut Vec<F>,
        coset:u32
    ){
        //Depth of the recursion, which is the base-2 logarithm of the length of the input.

        let rounds = coeffs.len().trailing_zeros();

        for r in 0..rounds{

            //Number of subproblems at this depth of the recursion.
            let parts = 1<<r;

            let normalising_value = vanishing_map(<BinaryField32b as TowerField>::basis(0, r as usize).unwrap(), r as usize).invert().unwrap();
            let coset_twiddle = vanishing_map(BinaryField32b::new(coset), r as usize)*normalising_value;

            for p in 0..parts{
                for b in 0..1<<(rounds-1 - r){
                    
                    //the butterfly in this ntt indexes (b*2^(r+1) + p) and (b*2^(r+1) + p + 2^r)
                    let left_idx = (b<<(r+1)) + p;
                    let right_idx = left_idx + (1<<r);
                    
                    let left_val = coeffs[left_idx];
                    coeffs[right_idx] += left_val;

                    let right_val = coeffs[right_idx];
                    coeffs[left_idx] += right_val * (self.twiddles[r as usize][b] + coset_twiddle);

                   
                }
            }


        }
    }
}

//Recursively generates twiddles.
pub fn compute_twiddles(log_transform_size:usize)->Vec<Vec<BinaryField32b>>{
    let s_evals = precompute_subspace_evals(log_transform_size);
    let s_evals_expanded = s_evals
        .iter()
        .enumerate()
        .map(|(i, s_evals_i)| {
            let mut expanded = Vec::with_capacity(1 << s_evals_i.len());
            expanded.push(BinaryField32b::ZERO);
            for &eval in s_evals_i.iter() {
                for i in 0..expanded.len() {
                    expanded.push(expanded[i] + eval);
                }
            }
            expanded
                }).collect();
    s_evals_expanded
}

fn precompute_subspace_evals(log_domain_size: usize) -> Vec<Vec<BinaryField32b>>{

	let mut s_evals = Vec::with_capacity(log_domain_size);

	// normalization_consts[i] = W_i(2^i)
	let mut normalization_consts = Vec::with_capacity(log_domain_size);
	normalization_consts.push(BinaryField32b::ONE);

	let s0_evals = (1..log_domain_size)
		.map(|i| <BinaryField32b as ExtensionField<BinaryField1b>>::basis(i).expect("basis vector must exist because of FieldTooSmall check above"))
		.collect::<Vec<_>>();
	s_evals.push(s0_evals);

	for _ in 1..log_domain_size {
		let (norm_const_i, s_i_evals) = {
			let norm_prev = *normalization_consts
				.last()
				.expect("normalization_consts is not empty");
			let s_prev_evals = s_evals.last().expect("s_evals is not empty");

			let norm_const_i = subspace_map(s_prev_evals[0], norm_prev);
			let s_i_evals = s_prev_evals
				.iter()
				.skip(1)
				.map(|&s_ij_prev| subspace_map(s_ij_prev, norm_prev))
				.collect::<Vec<_>>();

			(norm_const_i, s_i_evals)
		};

		normalization_consts.push(norm_const_i);
		s_evals.push(s_i_evals);
	}

	for (norm_const_i, s_evals_i) in normalization_consts.iter().zip(s_evals.iter_mut()) {
		let inv_norm_const = norm_const_i
			.invert()
			.expect("normalization constants are nonzero");
		for s_ij in s_evals_i.iter_mut() {
			*s_ij *= inv_norm_const;
		}
	}

	s_evals
}

fn subspace_map<F: Field>(elem: F, constant: F) -> F {
	elem.square() + constant * elem
}

pub fn vanishing_map(val:BinaryField32b, size:usize)->BinaryField32b{

    let mut res = BinaryField32b::ONE;

    if size != 0{
        for i in 0..1<<size{
            res *= val + BinaryField32b::from(i)
        }
    
    }


    res
}

#[test]

fn twiddles_test(){

    for i in 4..10{
    let twiddles = compute_twiddles(7);

    for j in 0..twiddles.len(){
            println!("{:?}  \n", twiddles[j])
    }
}
}

//Evaluates the polynomial naiively for testing purposes
pub fn poly_eval(coeffs:&Vec<BinaryField32b>, point:BinaryField32b)->BinaryField32b{
    let bits = (coeffs.len().trailing_zeros() +1 ) as usize;
    let mut res = BinaryField32b::ZERO;
    let mut basis_vals:Vec<BinaryField32b> = (0..bits).into_iter().map(|k|{
        vanishing_map(point, k) 
        *(vanishing_map(BinaryField32b::from(1<<k), k)).invert().unwrap()
    }).collect();

    basis_vals[0] = point;
    for i in 0..coeffs.len(){
        let mut basis_eval = BinaryField32b::ONE;
        for j in 0..bits{
            if ((i>>j)&1) == 1{
                basis_eval *= basis_vals[j as usize]
            }
        }

        res += basis_eval*coeffs[i];
    }

    res
}

#[test]
fn forward_test(){
    for i in 3..16{
        println!("Generating random polynomial \n");
        let time = Instant::now();

        
        let mut poly: Vec<BinaryField32b> = (0..1<<i).into_par_iter()
        .map(|_| {
            let mut rng = thread_rng();
            BinaryField32b::random(rng)
        }).collect();

        println!("Time taken {:?}", time.elapsed());


        let mut test_evals:Vec<_> = (0..1<<i).into_par_iter().map(|i| poly_eval(&poly, BinaryField32b::from(i))).collect();
        test_evals.par_chunks_mut(2).for_each(|pair| pair.swap(0, 1));
        println!("Computing twiddles");

        let ntt = AdditiveNTT::new(poly.len().trailing_zeros() as usize);

        println!("Computing forward transform");

        ntt.forward_ntt(&mut poly, 0);
        // ntt.inverse_ntt(&mut poly, 0);
        assert_eq!(test_evals, poly)
    }
}
