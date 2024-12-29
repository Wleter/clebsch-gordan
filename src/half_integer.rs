use core::fmt;
use std::{iter::Sum, ops::{Add, AddAssign, Neg, Sub, SubAssign}};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HalfI32 {
    doubled: i32
}

impl std::fmt::Debug for HalfI32 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.doubled & 1 == 0 {
            write!(f, "{}", self.doubled / 2)
        } else {
            write!(f, "{}/2", self.doubled)
        }
    }
}

impl fmt::Display for HalfI32 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.doubled & 1 == 0 {
            write!(f, "{}", self.doubled / 2)
        } else {
            write!(f, "{}/2", self.doubled)
        }
    }
}

impl HalfI32 {
    pub fn from_doubled(doubled: i32) -> Self {
        Self {
            doubled,
        }
    }

    pub fn double_value(&self) -> i32 {
        self.doubled
    }

    pub fn value(&self) -> f64 {
        self.doubled as f64 / 2.
    }
}

impl Add for HalfI32 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        HalfI32::from_doubled(self.doubled + other.doubled)
    }
}

impl Sub for HalfI32 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        HalfI32::from_doubled(self.doubled - other.doubled)
    }
}

impl AddAssign for HalfI32 {
    fn add_assign(&mut self, rhs: Self) {
        self.doubled += rhs.doubled;
    }
}

impl Sum for HalfI32 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(HalfI32::from_doubled(0), |acc, x| acc + x)
    }
}

impl SubAssign for HalfI32 {
    fn sub_assign(&mut self, rhs: Self) {
        self.doubled -= rhs.doubled;
    }
}

impl Neg for HalfI32 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            doubled: -self.doubled
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HalfU32 {
    doubled: u32
}

impl std::fmt::Debug for HalfU32 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.doubled & 1 == 0 {
            write!(f, "HalfU32: {}", self.doubled / 2)
        } else {
            write!(f, "HalfU32: {}/2", self.doubled)
        }
    }
}

impl fmt::Display for HalfU32 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.doubled & 1 == 0 {
            write!(f, "{}", self.doubled / 2)
        } else {
            write!(f, "{}/2", self.doubled)
        }
    }
}

impl HalfU32 {
    pub fn from_doubled(doubled: u32) -> Self {
        Self {
            doubled,
        }
    }

    pub fn double_value(&self) -> u32 {
        self.doubled
    }

    pub fn value(&self) -> f64 {
        self.doubled as f64 / 2.
    }
}

impl Add for HalfU32 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        HalfU32::from_doubled(self.doubled + other.doubled)
    }
}

impl Sub for HalfU32 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        HalfU32::from_doubled(self.doubled - other.doubled)
    }
}

impl AddAssign for HalfU32 {
    fn add_assign(&mut self, rhs: Self) {
        self.doubled += rhs.doubled;
    }
}

impl SubAssign for HalfU32 {
    fn sub_assign(&mut self, rhs: Self) {
        self.doubled -= rhs.doubled;
    }
}

impl Sum for HalfU32 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(HalfU32::from_doubled(0), |acc, x| acc + x)
    }
}

#[macro_export]
macro_rules! half_i32 {
    // Case for whole numbers
    ($value:literal) => {
        $crate::half_integer::HalfI32::from_doubled(2 * $value)
    };
    // Case for fractions like 5/2
    ($numerator:literal / 2) => {
        $crate::half_integer::HalfI32::from_doubled(($numerator))
    };
}

#[macro_export]
macro_rules! half_u32 {
    // Case for whole numbers
    ($value:literal) => {
        $crate::half_integer::HalfU32::from_doubled(2 * $value)
    };
    // Case for fractions like 5/2
    ($numerator:literal / 2) => {
        $crate::half_integer::HalfU32::from_doubled(($numerator))
    };
}

impl From<HalfI32> for HalfU32 {
    fn from(value: HalfI32) -> Self {
        HalfU32::from_doubled(value.doubled as u32)
    }
}

impl From<HalfU32> for HalfI32 {
    fn from(value: HalfU32) -> Self {
        HalfI32::from_doubled(value.doubled as i32)
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_half_integers() {
        let spin1 = half_i32!(3);
        assert!(spin1.double_value() == 6);

        let spin2 = half_u32!(3/2);
        assert!(spin2.double_value() == 3);

        assert!((spin1 + spin2.into()).double_value() == 9);
    }
}