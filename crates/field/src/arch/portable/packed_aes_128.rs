// Copyright 2024 Ulvetanna Inc.

use super::{
	packed::{
		impl_broadcast, impl_conversion, impl_packed_extension_field, packed_binary_field_tower,
		PackedPrimitiveType,
	},
	packed_arithmetic::{alphas, impl_tower_constants},
};
use crate::{
	aes_field::{
		AESTowerField128b, AESTowerField16b, AESTowerField32b, AESTowerField64b, AESTowerField8b,
	},
	arch::{PackedStrategy, PairwiseRecursiveStrategy, PairwiseTableStrategy},
	arithmetic_traits::{impl_invert_with, impl_mul_alpha_with, impl_mul_with, impl_square_with},
};

// Define 128 bit packed AES field types
pub type PackedAESBinaryField16x8b = PackedPrimitiveType<u128, AESTowerField8b>;
pub type PackedAESBinaryField8x16b = PackedPrimitiveType<u128, AESTowerField16b>;
pub type PackedAESBinaryField4x32b = PackedPrimitiveType<u128, AESTowerField32b>;
pub type PackedAESBinaryField2x64b = PackedPrimitiveType<u128, AESTowerField64b>;
pub type PackedAESBinaryField1x128b = PackedPrimitiveType<u128, AESTowerField128b>;

// Define conversion from type to underlier
impl_conversion!(u128, PackedAESBinaryField16x8b);
impl_conversion!(u128, PackedAESBinaryField8x16b);
impl_conversion!(u128, PackedAESBinaryField4x32b);
impl_conversion!(u128, PackedAESBinaryField2x64b);
impl_conversion!(u128, PackedAESBinaryField1x128b);

// Define tower
packed_binary_field_tower!(
	PackedAESBinaryField16x8b
	< PackedAESBinaryField8x16b
	< PackedAESBinaryField4x32b
	< PackedAESBinaryField2x64b
	< PackedAESBinaryField1x128b
);

// Define extension fields
impl_packed_extension_field!(PackedAESBinaryField16x8b);
impl_packed_extension_field!(PackedAESBinaryField8x16b);
impl_packed_extension_field!(PackedAESBinaryField4x32b);
impl_packed_extension_field!(PackedAESBinaryField2x64b);
impl_packed_extension_field!(PackedAESBinaryField1x128b);

// Define broadcast
impl_broadcast!(u128, AESTowerField8b);
impl_broadcast!(u128, AESTowerField16b);
impl_broadcast!(u128, AESTowerField32b);
impl_broadcast!(u128, AESTowerField64b);
impl_broadcast!(u128, AESTowerField128b);

// Define contants
// 0xD3 corresponds to 0x10 after isomorphism from BinaryField8b to AESField
impl_tower_constants!(AESTowerField8b, u128, 0x00d300d300d300d300d300d300d300d3);
impl_tower_constants!(AESTowerField16b, u128, { alphas!(u128, 4) });
impl_tower_constants!(AESTowerField32b, u128, { alphas!(u128, 5) });
impl_tower_constants!(AESTowerField64b, u128, { alphas!(u128, 6) });

// Define multiplication
impl_mul_with!(PackedAESBinaryField16x8b @ PairwiseTableStrategy);
impl_mul_with!(PackedAESBinaryField8x16b @ PairwiseRecursiveStrategy);
impl_mul_with!(PackedAESBinaryField4x32b @ PairwiseRecursiveStrategy);
impl_mul_with!(PackedAESBinaryField2x64b @ PairwiseRecursiveStrategy);
impl_mul_with!(PackedAESBinaryField1x128b @ PairwiseRecursiveStrategy);

// Define square
impl_square_with!(PackedAESBinaryField16x8b @ PairwiseTableStrategy);
impl_square_with!(PackedAESBinaryField8x16b @ PairwiseRecursiveStrategy);
impl_square_with!(PackedAESBinaryField4x32b @ PackedStrategy);
impl_square_with!(PackedAESBinaryField2x64b @ PackedStrategy);
impl_square_with!(PackedAESBinaryField1x128b @ PairwiseRecursiveStrategy);

// Define invert
impl_invert_with!(PackedAESBinaryField16x8b @ PairwiseTableStrategy);
impl_invert_with!(PackedAESBinaryField8x16b @ PairwiseRecursiveStrategy);
impl_invert_with!(PackedAESBinaryField4x32b @ PairwiseRecursiveStrategy);
impl_invert_with!(PackedAESBinaryField2x64b @ PairwiseRecursiveStrategy);
impl_invert_with!(PackedAESBinaryField1x128b @ PairwiseRecursiveStrategy);

// Define multiply by alpha
impl_mul_alpha_with!(PackedAESBinaryField16x8b @ PairwiseTableStrategy);
impl_mul_alpha_with!(PackedAESBinaryField8x16b @ PackedStrategy);
impl_mul_alpha_with!(PackedAESBinaryField4x32b @ PackedStrategy);
impl_mul_alpha_with!(PackedAESBinaryField2x64b @ PairwiseRecursiveStrategy);
impl_mul_alpha_with!(PackedAESBinaryField1x128b @ PairwiseRecursiveStrategy);
