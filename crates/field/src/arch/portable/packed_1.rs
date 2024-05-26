// Copyright 2024 Ulvetanna Inc.

use super::packed::{
	impl_broadcast, impl_conversion, impl_ops_for_zero_height, PackedPrimitiveType,
};
use crate::{
	underlier::{UnderlierType, U1},
	BinaryField1b,
};

// Define 1 bit packed field types
pub type PackedBinaryField1x1b = PackedPrimitiveType<U1, BinaryField1b>;

// Define conversion from type to underlier
impl_conversion!(U1, PackedBinaryField1x1b);

// Define broadcast
impl_broadcast!(U1, BinaryField1b);

// Define operations for height 0
impl_ops_for_zero_height!(PackedBinaryField1x1b);
