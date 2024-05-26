use std::{collections::HashMap};

use binius_field::BinaryField32b;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sha3::{self, Digest, Keccak256};



//Implementation for Merkle Tree commitments, the hashing algorithm is set to Keccak256 and can be made generic over choice of hasher.

//The data structure to construct the merkle tree is a hashmap whoes keys represent the layer of the tree, and the vector contains the nodes in the layer
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hash(pub Vec<u8>);
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MerkleTree{
    pub data: HashMap<usize, Vec<Hash>>
}

impl MerkleTree {

    pub fn new(
        leaf_hashes:Vec<Hash>
    )->Self{
        merklize(leaf_hashes)
    }
    //The commitment in the pcs is the root of the merkle tree.
    pub fn get_root(
        &self
    )->Hash{
        self.data.get(&0).unwrap()[0].clone()
    }
    pub fn get_merkle_path(
        &self,
        leaf_index:usize
    )->Vec<Hash>{
        get_merkle_path(&self.data, leaf_index)
    }

}


//Hashing a field element
pub fn hash_field(data: &Vec<BinaryField32b>)->Hash{
    let mut hash_state = Keccak256::new();

    data.iter().for_each(|d| hash_state.update(d.val().to_le_bytes()));

    Hash(hash_state.finalize().to_vec())
}


//Hashing a byte vector
pub fn hash(data: &Vec<u8>)->Hash{
   Hash(
    Keccak256::digest(data).to_vec()
   )
}

//Hashing a concatentation of previous hashes, required for Merkle Tree construction.
pub fn hash_concatenation(data1:&Hash, data2:&Hash)->Hash{

    let val = [data1.0.clone(), data2.0.clone()].concat();

    Hash(
        Keccak256::digest(val).to_vec()
    )
}


//Constructing a merkle tree

pub fn merklize(leaf_hashes:Vec<Hash>)->MerkleTree{
    assert!(leaf_hashes.len().is_power_of_two(), "Leaf layer's size needs to be a power of 2 to construct Merkle Tree.");
    
    let mut tree:HashMap<usize, Vec<Hash>> = HashMap::new();

    let tree_depth = leaf_hashes.len().trailing_zeros() as usize;

    tree.insert(tree_depth, leaf_hashes);

    for depth in (0..tree_depth).rev(){
        let lower_layer = tree.get(&(depth+1)).unwrap();
        
        let current_layer_size = lower_layer.len()/2;
        
        //Hashing the concatenations of 2*i and 2*i + 1 of the lower layer in parallel.

        let current_layer = 
        (0..current_layer_size).into_par_iter().map(|i| hash_concatenation(&lower_layer[2*i], &lower_layer[2*i + 1])).collect();
        
        tree.insert(depth, current_layer);


    }
    
    MerkleTree{
        data:tree
    }
}

//Get a merkle path from a given leaf index to the root, i.e the siblings required to concatenate and hash so we reach the root to prove membership.
pub fn get_merkle_path(tree:&HashMap<usize, Vec<Hash>>, leaf_index:usize)->Vec<Hash>
{
    let tree_depth = tree.len();

    let mut indices = vec![leaf_index;tree_depth -1];
    let mut path = Vec::new();

    //Here essentially as we move up from leaf to root, we check if the index is odd or even at the layer, and so the index of its sibling will be to its right or left respectively
    for d in 0..tree_depth-1{
        if (leaf_index>>d)&1==0{
            indices[d] = (leaf_index>>d) + 1;
        }
        else{
            indices[d] = (leaf_index>>d) - 1;
        }
    }

    for i in 0..indices.len(){
        path.push(
            match tree.get(&(tree_depth- i - 1)) {
                Some(tree_layer)=> tree_layer[indices[i]].clone(),
                None => panic!("Tried to index layer {:?} in merkle tree, seems to be empty.", tree_depth - i)
                
            }
        )
    }

    path
}

pub fn verify_merkle_path(commitment:&Hash, leaf_hash:Hash, leaf_index:usize, merkle_path:&Vec<Hash>){

    let mut hash = leaf_hash;

    for d in 0..merkle_path.len(){
        if (leaf_index>>d)&1 == 0 {
            hash = hash_concatenation(&hash, &merkle_path[d])
        }
        else{
            hash = hash_concatenation(&merkle_path[d], &hash)
        }
    }


    assert_eq!(hash, *commitment)
}
