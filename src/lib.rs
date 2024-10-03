mod primes;
mod rational;

use std::num::NonZeroUsize;

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

pub fn clear_wigner_3j_cache() {
    CACHED_WIGNER_3J.lock().clear();
}

/// Compute the Wigner 3j coefficient for the given `dj1`, `dj2`, `dj2`, `dm1`,
/// `dm2`, `dm3`.
pub fn wigner_3j(dj1: u32, dj2: u32, dj3: u32, dm1: i32, dm2: i32, dm3: i32) -> f64 {
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
pub fn clebsch_gordan(dj1: u32, dm1: i32, dj2: u32, dm2: i32, dj3: u32, dm3: i32) -> f64 {
    let mut w3j = wigner_3j(dj1, dj2, dj3, dm1, dm2, -dm3);

    w3j *= f64::sqrt((dj3 + 1) as f64);
    if ((dj1 as i32 - dj2 as i32 + dm3) / 2) & 1 == 1 {
        return -w3j;
    } else {
        return w3j;
    }
}

/// check the triangle condition on j1, j2, j3, i.e. `|j1 - j2| <= j3 <= j1 + j2`
fn triangle_condition(dj1: u32, dj2: u32, dj3: u32) -> bool {
    return (dj3 <= dj1 + dj2) && (dj1 <= dj2 + dj3) && (dj2 <= dj3 + dj1);
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

fn triangle_coefficient(dj1: u32, dj2: u32, dj3: u32) -> Rational {
    let n1 = factorial((dj1 + dj2 - dj3) / 2);
    let n2 = factorial((dj1 - dj2 + dj3) / 2);
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
        assert_ulps_eq!(wigner_3j(4, 12, 8, 0, 0, 2), 0.0);
        assert_ulps_eq!(wigner_3j(4, 12, 8, 0, 0, 0), f64::sqrt(715.0) / 143.0);
        assert_ulps_eq!(wigner_3j(10, 6, 4, -6, 6, 0), f64::sqrt(330.0) / 165.0);
        assert_ulps_eq!(wigner_3j(10, 6, 4, -4, 6, -2), -f64::sqrt(330.0) / 330.0);
        assert_ulps_eq!(wigner_3j(200, 200, 200, 200, -200, 0), 2.689688852311291e-13);

        assert_ulps_eq!(wigner_3j(0, 2, 2, 0, 0, 0), -0.5773502691896257);

        // https://github.com/Luthaf/wigners/issues/7
        assert_ulps_eq!(wigner_3j(200, 600, 570, 4, -4, 0), 0.001979165708981953);
    }

    #[test]
    fn test_clebsch_gordan() {
        // checked against sympy
        assert_ulps_eq!(clebsch_gordan(4, 0, 12, 0, 8, 2), 0.0);
        assert_ulps_eq!(clebsch_gordan(2, 2, 2, 2, 4, 4), 1.0);
        assert_ulps_eq!(clebsch_gordan(4, 4, 2, -2, 6, 2), f64::sqrt(1.0 / 15.0));

        // half spins
        assert_ulps_eq!(clebsch_gordan(1, 1, 1, -1, 2, 0), f64::sqrt(0.5));
        assert_ulps_eq!(clebsch_gordan(1, 1, 1, -1, 0, 0), f64::sqrt(0.5));
        assert_ulps_eq!(clebsch_gordan(1, -1, 1, 1, 2, 0), f64::sqrt(0.5));
        assert_ulps_eq!(clebsch_gordan(1, -1, 1, 1, 0, 0), -f64::sqrt(0.5));

        assert_ulps_eq!(clebsch_gordan(5, 3, 4, 2, 5, 5), -f64::sqrt(3. / 7.));
        assert_ulps_eq!(clebsch_gordan(5, 3, 3, 1, 6, 4), f64::sqrt(1. / 12.));
        assert_ulps_eq!(clebsch_gordan(5, 3, 3, 1, 4, 4), -f64::sqrt(8. / 21.));
    }
}
