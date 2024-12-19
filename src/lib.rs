mod primes;
mod rational;
pub mod half_integer;

use std::num::NonZeroUsize;

use half_integer::{HalfI32, HalfU32};
use num_rational::BigRational;
use parking_lot::Mutex;

use lru::LruCache;

use num_bigint::BigInt;
use num_traits::ToPrimitive;

use crate::primes::{factorial, PrimeFactorization};
use crate::rational::Rational;

// cache up to that many wigner_3j symbols in a LRU cache with 20_000 entries.
const WIGNER_3J_CACHE_SIZE: usize = 20_000;

type Wigner3jCacheKey = (i32, i32, i32, i32, i32);
lazy_static::lazy_static!(
    static ref CACHED_WIGNER_3J: Mutex<LruCache<Wigner3jCacheKey, f64>> = Mutex::new(
        LruCache::new(NonZeroUsize::new(WIGNER_3J_CACHE_SIZE).expect("cache size is zero"))
    );
);

// cache up to that many wigner_3j symbols in a LRU cache with 20_000 entries.
const WIGNER_6J_CACHE_SIZE: usize = 20_000;

type Wigner6jCacheKey = (u32, u32, u32, u32, u32, u32);
lazy_static::lazy_static!(
    static ref CACHED_WIGNER_6J: Mutex<LruCache<Wigner6jCacheKey, f64>> = Mutex::new(
        LruCache::new(NonZeroUsize::new(WIGNER_6J_CACHE_SIZE).expect("cache size is zero"))
    );
);

pub fn clear_wigner_3j_cache() {
    CACHED_WIGNER_3J.lock().clear();
}

/// Compute the Wigner 3j coefficient for the given `dj1`, `dj2`, `dj3`, `dm1`,
/// `dm2`, `dm3`.
pub fn wigner_3j(j1: HalfU32, j2: HalfU32, j3: HalfU32, m1: HalfI32, m2: HalfI32, m3: HalfI32) -> f64 {
    let dj1 = j1.double_value();
    let dj2 = j2.double_value();
    let dj3 = j3.double_value();
    let dm1 = m1.double_value();
    let dm2 = m2.double_value();
    let dm3 = m3.double_value();

    if dm1.unsigned_abs() > dj1 {
        panic!("invalid dj1/dm1 in wigner3j: {}/{}", dj1, dm1);
    } else if dm2.unsigned_abs() > dj2 {
        panic!("invalid dj2/dm2 in wigner3j: {}/{}", dj2, dm2);
    } else if dm3.unsigned_abs() > dj3 {
        panic!("invalid dj3/dm3 in wigner3j: {}/{}", dj3, dm3);
    }

    if (dj1 & 1 == 0) ^ (dm1 & 1 == 0) {
        panic!("invalid dj1/m1 in wigner3j: {}/{}", dj1, dm1)
    } else if (dj2 & 1 == 0) ^ (dm2 & 1 == 0) {
        panic!("invalid dj2/m2 in wigner3j: {}/{}", dj2, dm2)
    } else if (dj3 & 1 == 0) ^ (dm3 & 1 == 0) {
        panic!("invalid dj3/m3 in wigner3j: {}/{}", dj3, dm3)
    }

    if (dj1 + dj2 + dj3) & 1 == 1 {
        panic!("non compatible spins {dj1} {dj2} {dj3}")
    }

    if !triangle_condition(dj1, dj2, dj3) || dm1 + dm2 + dm3 != 0 {
        return 0.0;
    }

    let (dj1, dj2, dj3, dm1, dm2, _, mut sign) = reorder3j(dj1, dj2, dj3, dm1, dm2, dm3, 1.0);

    let total_j = (dj1 + dj2 + dj3) / 2;
    let alpha1 = (dj2 as i32 - dm1 - dj3 as i32) / 2;
    let alpha2 = (dj1 as i32 + dm2 - dj3 as i32) / 2;
    let beta1 = (dj1 + dj2 - dj3) as i32 / 2;
    let beta2 = (dj1 as i32 - dm1) / 2;
    let beta3 = (dj2 as i32 + dm2) / 2;

    // extra sign in definition: alpha1 - alpha2 = j1 + m2 - j2 + m1 = j1 - j2 + m3
    if (alpha1 - alpha2) % 2 != 0 {
        sign = -sign;
    }

    {
        let mut cache = CACHED_WIGNER_3J.lock();
        if let Some(&cached_value) = cache.get(&(alpha1, alpha2, beta1, beta2, beta3)) {
            return sign * cached_value;
        }
    }

    let s1 = triangle_coefficient(dj1, dj2, dj3);

    debug_assert!(beta2 >= 0);
    let mut s2 = factorial(beta2 as u32);

    debug_assert!((beta1 - alpha1) >= 0);
    s2 *= factorial((beta1 - alpha1) as u32);

    debug_assert!((beta1 - alpha2) >= 0);
    s2 *= factorial((beta1 - alpha2) as u32);

    debug_assert!(beta3 >= 0);
    s2 *= factorial(beta3 as u32);

    debug_assert!((beta3 - alpha1) >= 0);
    s2 *= factorial((beta3 - alpha1) as u32);

    debug_assert!((beta2 - alpha2) >= 0);
    s2 *= factorial((beta2 - alpha2) as u32);

    let (series_numerator, series_denominator) = compute_3j_series(total_j, beta1, beta2, beta3, alpha1, alpha2);

    let numerator = s1.numerator * s2;
    let mut s = Rational::new(numerator, s1.denominator);

    let series_denominator = Rational::new(PrimeFactorization::one(), series_denominator);

    // insert series denominator in the root, this improves precision compared
    // to immediately converting the full series to f64
    s *= &series_denominator;
    s *= &series_denominator;
    s.simplify();

    let result = series_numerator * s.signed_root();

    {
        let mut cache = CACHED_WIGNER_3J.lock();
        cache.put((alpha1, alpha2, beta1, beta2, beta3), result);
    }

    return sign * result;
}

/// Compute the Clebsch-Gordan coefficient <j1 m1 ; j2 m2 | j3 m3> using their
/// relation to Wigner 3j coefficients:
///
/// ```text
/// <j1 m1 ; j2 m2 | j3 m3> = (-1)^(j1 - j2 + m3) sqrt(2*j3 + 1) wigner_3j(j1, j2, j3, m1, m2, -m3)
/// ```
pub fn clebsch_gordan(j1: HalfU32, m1: HalfI32, j2: HalfU32, m2: HalfI32, j3: HalfU32, m3: HalfI32) -> f64 {
    let mut w3j = wigner_3j(j1, j2, j3, m1, m2, -m3);
    w3j *= f64::sqrt((j3.double_value() + 1) as f64);

    let sign_criterion: HalfI32 = m3 + j1.into() - j2.into();
    if (sign_criterion.double_value() / 2) & 1 == 1 {
        return -w3j;
    } else {
        return w3j;
    }
}

/// Compute the Wigner 6j coefficient for the given `j1`, `j2`, `j3`, `j4`,
/// `j5`, `j6`.
pub fn wigner_6j(j1: HalfU32, j2: HalfU32, j3: HalfU32, j4: HalfU32, j5: HalfU32, j6: HalfU32) -> f64 {
    let dj1 = j1.double_value();
    let dj2 = j2.double_value();
    let dj3 = j3.double_value();
    let dj4 = j4.double_value();
    let dj5 = j5.double_value();
    let dj6 = j6.double_value();

    if !triangle_condition(dj1, dj2, dj3) 
        || !triangle_condition(dj1, dj5, dj6)
        || !triangle_condition(dj4, dj2, dj6)
        || !triangle_condition(dj4, dj5, dj3) 
    {
        return 0.0;
    }

    let a1 = (dj1 + dj2 + dj3) / 2;
    let a2 = (dj1 + dj6 + dj5) / 2;
    let a3 = (dj2 + dj4 + dj6) / 2;
    let a4 = (dj3 + dj4 + dj5) / 2;
    let b1 = (dj1 + dj2 + dj4 + dj5) / 2;
    let b2 = (dj1 + dj3 + dj4 + dj6) / 2;
    let b3 = (dj2 + dj3 + dj5 + dj6) / 2;

    let (b1, b2, b3, a1, a2, a3, a4) = reorder6j(b1, b2, b3, a1, a2, a3, a4);

    {
        let mut cache = CACHED_WIGNER_6J.lock();
        if let Some(&cached_value) = cache.get(&(b1, b2, b3, a1, a2, a3)) {
            return cached_value;
        }
    }

    let mut ratio = triangle_coefficient(dj1, dj2, dj3);
    ratio *= triangle_coefficient(dj1, dj6, dj5);
    ratio *= triangle_coefficient(dj2, dj4, dj6);
    ratio *= triangle_coefficient(dj3, dj4, dj5);

    let (series_numerator, series_denominator) = compute_6j_series(b1 as i32, b2 as i32, b3 as i32, a1 as i32, a2 as i32, a3 as i32, a4 as i32);

    let series_denominator = Rational::new(PrimeFactorization::one(), series_denominator);
    ratio *= &series_denominator;
    ratio *= &series_denominator;
    ratio.simplify();

    let (s_num, r_num) = ratio.numerator.split_square();
    let (s_den, r_den) = ratio.denominator.split_square();

    let s = BigRational::new(s_num.as_bigint() * series_numerator, s_den.as_bigint());
    let r = Rational::new(r_num, r_den);

    let result = s.to_f64().expect("Too large value to convert to f64 in wigner6j")
        * r.numerator.as_sqrt() / r.denominator.as_sqrt();

    {
        let mut cache = CACHED_WIGNER_6J.lock();
        cache.put((b1, b2, b3, a1, a2, a3), result);
    }

    result
}

/// check the triangle condition on j1, j2, j3, i.e. `|j1 - j2| <= j3 <= j1 + j2`
fn triangle_condition(dj1: u32, dj2: u32, dj3: u32) -> bool {
    return (dj3 <= dj1 + dj2) && (dj1 <= dj2 + dj3) && (dj2 <= dj3 + dj1) && (dj1 + dj2 + dj3) % 2 == 0;
}

// reorder j1/m1, j2/m2, j3/m3 such that j1 >= j2 >= j3 and m1 >= 0 or m1 == 0 && m2 >= 0
fn reorder3j(dj1: u32, dj2: u32, dj3: u32, dm1: i32, dm2: i32, dm3: i32, mut sign: f64) -> (u32, u32, u32, i32, i32, i32, f64) {
    if dj1 < dj2 {
        return reorder3j(dj2, dj1, dj3, dm2, dm1, dm3, -sign);
    } else if dj2 < dj3 {
        return reorder3j(dj1, dj3, dj2, dm1, dm3, dm2, -sign);
    } else if dm1 < 0 || (dm1 == 0 && dm2 < 0) {
        return reorder3j(dj1, dj2, dj3, -dm1, -dm2, -dm3, -sign);
    } else {
        // sign doesn't matter if total J = j1 + j2 + j3 is even
        if (dj1 + dj2 + dj3) % 4 == 0 {
            sign = 1.0;
        }
        return (dj1, dj2, dj3, dm1, dm2, dm3, sign);
    }
}

fn reorder6j(b1: u32, b2: u32, b3: u32, a1: u32, a2: u32, a3: u32, a4: u32) -> (u32, u32, u32, u32, u32, u32, u32){
    if b1 < b2 {
        reorder6j(b2, b1, b3, a1, a2, a3, a4)
    } else if b2 < b3 {
        reorder6j(b1, b3, b2, a1, a2, a3, a4)
    } else if a1 < a2 {
        reorder6j(b1, b2, b3, a2, a1, a3, a4)
    } else if a2 < a3 {
        reorder6j(b1, b2, b3, a1, a3, a2, a4)
    } else if a3 < a4 {
        reorder6j(b1, b2, b3, a1, a2, a4, a3)
    } else {
        return (b1, b2, b3, a1, a2, a3, a4)
    }
}

fn triangle_coefficient(dj1: u32, dj2: u32, dj3: u32) -> Rational {
    let n1 = factorial((dj1 + dj2 - dj3) / 2);
    let n2 = factorial((dj1 + dj3 - dj2) / 2);
    let n3 = factorial((dj2 + dj3 - dj1) / 2);
    let numerator = n1 * n2 * n3;
    let denominator = factorial((dj1 + dj2 + dj3) / 2 + 1);

    let mut result = Rational::new(numerator, denominator);
    result.simplify();
    return result;
}

fn max(a: i32, b: i32, c: i32) -> i32 {
    std::cmp::max(a, std::cmp::max(b, c))
}

fn min(a: i32, b: i32, c: i32) -> i32 {
    std::cmp::min(a, std::cmp::min(b, c))
}

/// compute the sum appearing in the 3j symbol
fn compute_3j_series(total_j: u32, beta1: i32, beta2: i32, beta3: i32, alpha1: i32, alpha2: i32) -> (f64, PrimeFactorization) {
    let range = max(alpha1, alpha2, 0)..(min(beta1, beta2, beta3) + 1);

    let mut numerators = Vec::with_capacity(range.len());
    let mut denominators = Vec::with_capacity(range.len());
    for k in range {
        let numerator = if k % 2 == 0 {
            PrimeFactorization::one()
        } else {
            PrimeFactorization::minus_one()
        };
        numerators.push(numerator);

        debug_assert!(k >= 0);
        let mut denominator = factorial(k as u32);

        debug_assert!((k - alpha1) >= 0);
        denominator *= factorial((k - alpha1) as u32);

        debug_assert!((k - alpha2) >= 0);
        denominator *= factorial((k - alpha2) as u32);

        debug_assert!((beta1 - k) >= 0);
        denominator *= factorial((beta1 - k) as u32);

        debug_assert!((beta2 - k) >= 0);
        denominator *= factorial((beta2 - k) as u32);

        debug_assert!((beta3 - k) >= 0);
        denominator *= factorial((beta3 - k) as u32);

        denominators.push(denominator);
    }

    let denominator = common_denominator(&mut numerators, &denominators);

    let numerator = if total_j > 100 {
        // For large total J, we will overflow f64 in this sum, but performing
        // the sum with big integers is enough to recover the full precision
        let mut numerator = BigInt::from(0);
        for num in numerators {
            numerator += num.as_bigint();
        }
        numerator.to_f64().expect("not a f64")
    } else {
        let mut numerator = 0.0;
        for num in numerators {
            numerator += num.as_f64();
        }
        numerator
    };

    return (numerator, denominator);
}

fn compute_6j_series(b1: i32, b2: i32, b3: i32, a1: i32, a2: i32, a3: i32, a4: i32) -> (BigInt, PrimeFactorization) {
    let range = max(a1, a2, a3).max(a4)..(min(b1, b2, b3) + 1);

    let mut numerators = Vec::with_capacity(range.len());
    let mut denominators = Vec::with_capacity(range.len());
    for k in range {
        let mut numerator = if k % 2 == 0 {
            PrimeFactorization::one()
        } else {
            PrimeFactorization::minus_one()
        };
        numerator *= factorial(k as u32 + 1);
        numerators.push(numerator);

        assert!((k - a1) >= 0);
        let mut denominator = factorial((k - a1) as u32);

        debug_assert!((k - a2) >= 0);
        denominator *= factorial((k - a2) as u32);

        debug_assert!((k - a3) >= 0);
        denominator *= factorial((k - a3) as u32);
        
        debug_assert!((k - a4) >= 0);
        denominator *= factorial((k - a4) as u32);

        debug_assert!((b1 - k) >= 0);
        denominator *= factorial((b1 - k) as u32);

        debug_assert!((b2 - k) >= 0);
        denominator *= factorial((b2 - k) as u32);

        debug_assert!((b3 - k) >= 0);
        denominator *= factorial((b3 - k) as u32);

        denominators.push(denominator);
    }

    let denominator = common_denominator(&mut numerators, &denominators);

    let mut numerator = BigInt::from(0);
    for num in numerators {
        numerator += num.as_bigint();
    }

    return (numerator, denominator);
}

/// Given a list of numerators and denominators, compute the common denominator
/// and the rescaled numerator, putting all fractions at the same common
/// denominator
fn common_denominator(
    numerators: &mut [PrimeFactorization],
    denominators: &[PrimeFactorization]
) -> PrimeFactorization {
    debug_assert_eq!(numerators.len(), denominators.len());
    if denominators.is_empty() {
        return PrimeFactorization::one()
    }

    let mut denominator = denominators[0].clone();
    for other in denominators.iter().skip(1) {
        denominator.least_common_multiple(other);
    }

    // rescale numerators
    for (num, den) in numerators.iter_mut().zip(denominators.iter()) {
        *num *= &denominator;
        *num /= den;
    }

    return denominator;
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_ulps_eq;

    #[test]
    fn test_wigner3j() {
        // checked against sympy
        assert_ulps_eq!(
            wigner_3j(half_u32!(2), half_u32!(6), half_u32!(4), 
                      half_i32!(0), half_i32!(0), half_i32!(2)), 
            0.0
        );
        assert_ulps_eq!(
            wigner_3j(half_u32!(2), half_u32!(6), half_u32!(4), 
                      half_i32!(0), half_i32!(0), half_i32!(0)), 
            f64::sqrt(715.0) / 143.0
        );
        assert_ulps_eq!(
            wigner_3j(half_u32!(5), half_u32!(3), half_u32!(2), 
                      half_i32!(-3), half_i32!(3), half_i32!(0)), 
            f64::sqrt(330.0) / 165.0
        );
        assert_ulps_eq!(
            wigner_3j(half_u32!(5), half_u32!(3), half_u32!(2), 
                      half_i32!(-2), half_i32!(3), half_i32!(-1)), 
            -f64::sqrt(330.0) / 330.0
        );
        assert_ulps_eq!(
            wigner_3j(half_u32!(100), half_u32!(100), half_u32!(100), 
                      half_i32!(100), half_i32!(-100), half_i32!(0)), 
            2.689688852311291e-13
        );
        assert_ulps_eq!(
            wigner_3j(half_u32!(0), half_u32!(1), half_u32!(1), 
                      half_i32!(0), half_i32!(0), half_i32!(0)),
            -0.5773502691896257
        );
        // https://github.com/Luthaf/wigners/issues/7
        assert_ulps_eq!(
            wigner_3j(half_u32!(100), half_u32!(300), half_u32!(285), 
                      half_i32!(2), half_i32!(-2), half_i32!(0)), 
            0.001979165708981953
        );
    }

    #[test]
    fn test_clebsch_gordan() {
        // checked against sympy
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(2), half_i32!(0), 
                           half_u32!(6), half_i32!(0), 
                           half_u32!(4), half_i32!(1)), 
            0.0
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(1), half_i32!(1), 
                           half_u32!(1), half_i32!(1), 
                           half_u32!(2), half_i32!(2)), 
            1.0
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(2), half_i32!(2), 
                           half_u32!(1), half_i32!(-1), 
                           half_u32!(3), half_i32!(1)), 
            f64::sqrt(1.0 / 15.0)
        );
        // half spins
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(1/2), half_i32!(1/2), 
                           half_u32!(1/2), half_i32!(-1/2), 
                           half_u32!(1), half_i32!(0)), 
            f64::sqrt(0.5)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(1/2), half_i32!(1/2), 
                           half_u32!(1/2), half_i32!(-1/2), 
                           half_u32!(0), half_i32!(0)), 
            f64::sqrt(0.5)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(1/2), half_i32!(-1/2), 
                           half_u32!(1/2), half_i32!(1/2), 
                           half_u32!(1), half_i32!(0)), 
            f64::sqrt(0.5)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(1/2), half_i32!(-1/2), 
                           half_u32!(1/2), half_i32!(1/2), 
                           half_u32!(0), half_i32!(0)), 
            -f64::sqrt(0.5)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(5/2), half_i32!(3/2), 
                           half_u32!(2), half_i32!(1), 
                           half_u32!(5/2), half_i32!(5/2)), 
            -f64::sqrt(3. / 7.)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(5/2), half_i32!(3/2), 
                           half_u32!(3/2), half_i32!(1/2), 
                           half_u32!(3), half_i32!(2)), 
            f64::sqrt(1. / 12.)
        );
        assert_ulps_eq!(
            clebsch_gordan(half_u32!(5/2), half_i32!(3/2), 
                           half_u32!(3/2), half_i32!(1/2), 
                           half_u32!(2), half_i32!(2)), 
            -f64::sqrt(8. / 21.)
        );
    }

    #[test]
    fn test_wigner6j() {
        assert_ulps_eq!(
            wigner_6j(half_u32!(1), half_u32!(1), half_u32!(1), 
                      half_u32!(1), half_u32!(1), half_u32!(1)), 
            1. / 6.
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(1), half_u32!(2), half_u32!(3), 
                      half_u32!(3), half_u32!(2), half_u32!(1)), 
            f64::sqrt(14.) / 35.
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(3), half_u32!(3), half_u32!(3), 
                      half_u32!(3), half_u32!(3), half_u32!(3)), 
            -1. / 14.
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(5), half_u32!(5), half_u32!(5), 
                      half_u32!(5), half_u32!(5), half_u32!(5)), 
            1. / 52.
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(8), half_u32!(8), half_u32!(8), 
                      half_u32!(8), half_u32!(8), half_u32!(8)), 
            -0.01265208072315355
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(64), half_u32!(10), half_u32!(64), 
                      half_u32!(64), half_u32!(0), half_u32!(64)), 
            1. / 129.
        );
        assert_ulps_eq!(
            wigner_6j(half_u32!(1/2), half_u32!(1), half_u32!(1/2), 
                      half_u32!(1/2), half_u32!(0), half_u32!(1/2)), 
            0.5
        );
    }
}
