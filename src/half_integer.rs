macro_rules! impl_half_integers {
    ($($name:ident => $underlying:ty),*) => {
        $(
            #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
            pub struct $name {
                doubled: $underlying
            }

            impl $name {
                pub fn new(value: f64) -> Self {
                    Self {
                        doubled: (2.0 * value) as $underlying
                    }
                }

                pub fn from_doubled(doubled: $underlying) -> Self {
                    Self {
                        doubled,
                    }
                }
            
                pub fn double_value(&self) -> $underlying {
                    self.doubled
                }
            
                pub fn value(&self) -> f64 {
                    self.doubled as f64 / 2.
                }
            }

            impl std::fmt::Debug for $name {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    if self.doubled & 1 == 0 {
                        write!(f, "{}", self.doubled / 2)
                    } else {
                        write!(f, "{}/2", self.doubled)
                    }
                }
            }

            impl std::fmt::Display for $name {
                fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                    if self.doubled & 1 == 0 {
                        write!(f, "{}", self.doubled / 2)
                    } else {
                        write!(f, "{}/2", self.doubled)
                    }
                }
            }

            
            impl std::ops::Add for $name {
                type Output = Self;
                fn add(self, other: Self) -> Self {
                    $name::from_doubled(self.doubled + other.doubled)
                }
            }

            impl std::ops::Sub for $name {
                type Output = Self;
                fn sub(self, other: Self) -> Self {
                    $name::from_doubled(self.doubled - other.doubled)
                }
            }

            impl std::ops::AddAssign for $name {
                fn add_assign(&mut self, rhs: Self) {
                    self.doubled += rhs.doubled;
                }
            }

            impl std::iter::Sum for $name {
                fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                    iter.fold($name::from_doubled(0), |acc, x| acc + x)
                }
            }

            impl std::ops::SubAssign for $name {
                fn sub_assign(&mut self, rhs: Self) {
                    self.doubled -= rhs.doubled;
                }
            }
        )*
    };
}

macro_rules! impl_signed {
    ($($name:ident => $underlying:ty),*) => {
        $(
            impl std::ops::Neg for $name {
                type Output = Self;

                fn neg(self) -> Self::Output {
                    Self {
                        doubled: -self.doubled
                    }
                }
            }
        )*
    }
}

impl_half_integers!(HalfI32 => i32, HalfU32 => u32);
impl_signed!(HalfI32 => i32);

impl From<HalfI32> for HalfU32 {
    fn from(value: HalfI32) -> Self {
        HalfU32::from_doubled(value.doubled as u32)
    }
}

impl From<u32> for HalfU32 {
    fn from(value: u32) -> Self {
        HalfU32::from_doubled(2 * value)
    }
}

impl From<HalfU32> for HalfI32 {
    fn from(value: HalfU32) -> Self {
        HalfI32::from_doubled(value.doubled as i32)
    }
}

impl From<i32> for HalfI32 {
    fn from(value: i32) -> Self {
        HalfI32::from_doubled(2 * value)
    }
}

#[macro_export]
macro_rules! hi32 {
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
macro_rules! hu32 {
    // Case for whole numbers
    ($value:literal) => {
        $crate::half_integer::HalfU32::from_doubled(2 * $value)
    };
    // Case for fractions like 5/2
    ($numerator:literal / 2) => {
        $crate::half_integer::HalfU32::from_doubled(($numerator))
    };
}

#[cfg(test)]
mod test {
    use crate::half_integer::{HalfI32, HalfU32};

    #[test]
    fn test_half_integers() {
        let spin = HalfU32::new(1.5);
        assert!(spin.double_value() == 3);
        
        let spin: HalfU32 = 2.into();
        assert!(spin.double_value() == 4);

        let spin = HalfI32::new(1.5);
        assert!(spin.double_value() == 3);
        
        let spin: HalfI32 = 2.into();
        assert!(spin.double_value() == 4);

        let spin1 = hi32!(3);
        assert!(spin1.double_value() == 6);

        let spin2 = hu32!(3/2);
        assert!(spin2.double_value() == 3);

        assert!((spin1 + spin2.into()).double_value() == 9);
    }
}