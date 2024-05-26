use binius_field::{BinaryField128b, BinaryField32b, Field};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::utils::{packed_arithmetic::PackedAlgebra32, ntt::AdditiveNTT};

const RATE:usize =  4;
const PACKING_DEGREE:usize = 5;


pub fn encode(message:&Vec<BinaryField32b>, ntt:&AdditiveNTT)->Vec<BinaryField32b>{

    let mut code = message.clone();

    let mut inverse = message.clone();

    ntt.inverse_ntt(&mut inverse, 0);
    for i in 1..RATE{
        let mut encode = inverse.clone();
        ntt.forward_ntt(&mut encode, (i*message.len()) as u32);
        code.append(&mut encode.clone());
    }

    code
}

pub fn encode_extension(message:&Vec<BinaryField128b>, ntt:&AdditiveNTT)->Vec<PackedAlgebra32>{

    let packed_message = PackedAlgebra32::pack(message.clone());
    let mut code = packed_message.clone();

    let mut inverse = packed_message.clone();

    ntt.inverse_ntt(&mut inverse, 0);
    for i in 1..RATE{
        let mut encode = inverse.clone();
        ntt.forward_ntt(&mut encode, (i*message.len()) as u32);
        code.append(&mut encode.clone());
    }

    code
}


pub fn encode_interleaved(poly: &Vec<Vec<BinaryField32b>>, ntt:&AdditiveNTT, rows:usize, cols:usize)->Vec<Vec<BinaryField32b>>{
    let code:Vec<Vec<BinaryField32b>> = (0..rows).into_par_iter().map(|row| encode(&poly[row], ntt)).collect();
    code
}

#[derive(Debug)]
pub struct Code{
    code:Vec<Vec<BinaryField32b>>,
    pub rows: usize,
    pub cols: usize
}

impl Code{
    pub fn new(
        poly: &Vec<BinaryField32b>,
        ntt:&AdditiveNTT
    )->Code{

        let variables = poly.len().trailing_zeros() as usize;
        let mut log_cols = (1<<PACKING_DEGREE) + PACKING_DEGREE - (RATE.trailing_zeros() as usize);

        if variables<log_cols{
            log_cols = (variables+1)/2
        }
        let log_rows = variables - log_cols;

        let cols = 1<<log_cols;
        let rows = 1<<log_rows;

        let coeff_matrix:Vec<Vec<BinaryField32b>> = make_coeff_matrix(poly, cols);
        let code = encode_interleaved(&coeff_matrix, ntt,rows, cols);
        Code{
            code,
            rows,
            cols
        }
    }

    pub fn make_linear_combination(
        &self,
        scalars:Vec<BinaryField128b>
    )->Vec<BinaryField128b>{

        assert_eq!(self.rows, scalars.len());

        let mut result = vec![BinaryField128b::ZERO; self.cols*32];

        result.par_iter_mut().chunks(32).enumerate()
        .for_each(|(col, mut entry)|{

            for row in 0..self.rows{
                for bit in 0..32{
                    if (self.code[row][col].val()>>bit)&1 == 1{
                    *entry[bit] += scalars[row];
                    }
                }
            }
        }
        );
        result
    }

    pub fn col(
        &self,
        col:usize
    )->Vec<BinaryField32b>{

        assert!(col < self.cols);

        let mut result = vec![BinaryField32b::ZERO; self.rows];

        result.par_iter_mut().enumerate()
        .for_each(|(row, entry)|{
                *entry = self.code[row][col];
        }
        );
        result
    }

}

pub fn make_coeff_matrix(poly: &Vec<BinaryField32b>, cols:usize)->Vec<Vec<BinaryField32b>> {
    poly.chunks(cols).map(|row| row.to_vec()).collect()
}

pub fn make_linear_combination(
    poly: Vec<Vec<BinaryField32b>>,
    scalars:Vec<BinaryField128b>
)->Vec<BinaryField128b>{

    assert_eq!(poly.len(), scalars.len());

    let mut result = vec![BinaryField128b::ZERO; poly[0].len()*32];

    result.par_iter_mut().chunks(32).enumerate()
    .for_each(|(col, mut entry)|{

        for row in 0..poly.len(){
            for bit in 0..32{
                if (poly[row][col].val()>>bit)&1 == 1{
                *entry[bit] += scalars[row];
                }
            }
        }
    }
    );
    result
}

//Computes the Fourier coefficients/Lagrange basis evaluations at a random point.
pub fn compute_fourier_bases(r: &Vec<BinaryField128b>) -> Vec<BinaryField128b> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc_eq = [BinaryField128b::ONE - r[0], r[0]].to_vec();
    //Iterate over the length of the r vector
    for k in 1..r.len() {
        let temp = fc_eq;
        //initialize fc_eq of double size with zero

        fc_eq = vec![BinaryField128b::ZERO; temp.len() * 2];
        for iter in 0..temp.len(){
            fc_eq[2*iter+1] = temp[iter]*r[k];
            fc_eq[2*iter] = temp[iter] - fc_eq[2*iter+1];

        }
}
    fc_eq
}