use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};
use paste::paste;
use binius_field::{BinaryField128b, BinaryField2b,BinaryField16b, BinaryField8b, BinaryField4b, BinaryField32b, Field};
use rand::thread_rng;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use concat_arrays::concat_arrays;

//Implementation of the packed algebra required as part of the block level testing scheme, essentially requiring us to perform operations with vectors of 32 F_128 elements as if they were elements in the F_32 extension field.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PackedAlgebra32(pub [BinaryField128b;32]);


impl PackedAlgebra32 {
    pub fn new(vec: [BinaryField128b;32])->PackedAlgebra32{
        PackedAlgebra32(vec)
    }

    pub fn pack(vec: Vec<BinaryField128b>)->Vec<PackedAlgebra32>{

        vec.par_iter().copied().chunks(32)
        .map(|chunk| PackedAlgebra32(<[BinaryField128b;32]>::try_from(chunk).unwrap())).collect()
    }

    pub fn unpack(vec:Vec<PackedAlgebra32>)->Vec<BinaryField128b>{
        vec.par_iter().map(|pack| pack.0.to_vec()).flatten().collect()
    }

}

impl Mul<BinaryField32b> for PackedAlgebra32 {

    type Output = Self;

    fn mul(self, rhs:BinaryField32b) -> Self::Output {

        PackedAlgebra32(PackedMul32(&self.0, rhs))
    }

}


impl Add for PackedAlgebra32 {

    type Output = Self;

    fn add(self, rhs:PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl Sub for PackedAlgebra32 {

    type Output = Self;

    fn sub(self, rhs:PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl Add<&Self> for PackedAlgebra32 {

    type Output = Self;

    fn add(self, rhs:&PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl Sub<&Self> for PackedAlgebra32 {

    type Output = Self;

    fn sub(self, rhs:&PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl Add<&mut Self> for PackedAlgebra32 {

    type Output = Self;

    fn add(self, rhs:&mut PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl Sub<&mut Self> for PackedAlgebra32 {

    type Output = Self;

    fn sub(self, rhs:&mut PackedAlgebra32) -> Self::Output {
        PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }

}

impl AddAssign for PackedAlgebra32 {
    fn add_assign(&mut self, rhs: Self) {
        *self = PackedAlgebra32(packed_tower_add_32(&self.0, &rhs.0))
    }
}

impl SubAssign for PackedAlgebra32 {
    fn sub_assign(&mut self, rhs: Self) {
        *self += rhs
    }
}



pub fn PackedMul2(lhs: &[BinaryField128b;2], rhs: BinaryField2b)->[BinaryField128b;2]{



    if u8::from(rhs.val()) == 0{
        [BinaryField128b::ZERO; 2]
    }
    else if u8::from(rhs.val()) == 1 {
         *lhs
    }
    else if u8::from(rhs.val()) == 2 {
        [lhs[1], lhs[0] + lhs[1]]
    }
    else {
        [lhs[0] + lhs[1], lhs[0]]
    }

}

pub fn PackedMulAlpha2(lhs: [BinaryField128b;2] )->[BinaryField128b;2]{

    [lhs[1], lhs[0] + lhs[1]]
}

pub fn packed_tower_add_32(lhs: &[BinaryField128b;32], rhs: &[BinaryField128b;32])->[BinaryField128b;32]
{
    let mut res = [BinaryField128b::ZERO; 32];
    (0..32).for_each(|i| res[i] = lhs[i] + rhs[i]);
    res
}




macro_rules! PackedTowerAlgebra {
    ($field:ty, $($level:expr),* => $( $next:expr),* ) => {
        paste!{
        $(

            pub fn [<PackedMulAlpha $level>](lhs:[$field; $level ])->[$field; $level]{
                let (a0, a1) = lhs.split_at(lhs.len()>>1);

                let a0 = <[$field; $next]>::try_from(a0).unwrap();
                let a1 = <[$field; $next]>::try_from(a1).unwrap();

                let t1 = [<PackedMulAlpha $next>](a1.clone());
                let t0 = [<packed_tower_add_ $next>](&a0, &t1);

                concat_arrays!(a1, t0)


            }

            pub fn [<packed_tower_add_ $next>](lhs: &[$field;$next], rhs: &[$field;$next])->[$field;$next]
            {
                let mut res = [$field::ZERO;$next];
                (0..$next).for_each(|i| res[i] = lhs[i] + rhs[i]);
                res
            }

            pub fn [<PackedMul $level>](lhs:&[$field ; $level], rhs: [<BinaryField $level b>])->[$field; $level]{
                let (a_0, a_1) = lhs.split_at($level>>1);
                let (b_0, b_1) = rhs.into();
                let a_0 = <[$field; $next]>::try_from(a_0).unwrap();
                let a_1 = <[$field; $next]>::try_from(a_1).unwrap();
                let z0 = [<PackedMul $next>]( &a_0, b_0);
                let z1 = [<PackedMul $next>]( &a_1, b_1);

                let t0 = [<packed_tower_add_ $next>](&a_0,&a_1);

                let t1 = [<packed_tower_add_ $next>](&z0,&z1);

                let t2 = [<PackedMulAlpha $next>](z1);
                let z3 = [<PackedMul $next>](
                    &t0,
                    (b_0+ b_1)
                );

                let z3 = [<packed_tower_add_ $next>](&z3, &t1);
                let z3 = [<packed_tower_add_ $next>](&z3, &t2);

                concat_arrays![t1,z3]

            }

        )*
        }
    };
}

PackedTowerAlgebra!{BinaryField128b, 32,16,8,4 => 16,8,4,2}



#[test]
pub fn toweralg(){
    let mut rng = thread_rng();
    let a = BinaryField32b::random(&mut rng);

    let mut b:[BinaryField128b;32] = [BinaryField128b::ZERO;32];

    (0..32).into_iter().for_each(|v| b[v] = BinaryField128b::from(((a.val()>>v) &1) as u128));

    let rand = BinaryField32b::random(&mut rng);
    println!("{:b} \n", a.val());
    for i in b.iter(){
        print!("{:?}", i.val())
    }

    println!("\n");

    println!("{:b} \n", rand.val());



    let test = a*rand;
    let res = PackedMul32(&b, rand);

    println!("{:b} \n", test.val());

    for i in res.iter(){
        print!("{:?}", i.val())
    }
    println!("\n")
}


//Evaluates a vector packed coeffs, i.e each BinaryField32b actually represents 32 coefficients of the based field at a time, at a random point.
pub fn evaluate_unpacked(poly:&Vec<BinaryField32b>, basis:&Vec<BinaryField128b>)->BinaryField128b{
    assert_eq!(poly.len()*32, basis.len());


    basis.par_iter().chunks(32).zip(poly.par_iter()).map(|(basis_chunk, packed_coeff)|
    {
    let mut acc = BinaryField128b::ZERO;

    for i in 0..32{

        if (packed_coeff.val()>>i)&1 == 1{
            acc+=basis_chunk[i]
        }
    }
    acc
    }
    ).reduce(||BinaryField128b::ZERO, |acc,g| acc+g)




}