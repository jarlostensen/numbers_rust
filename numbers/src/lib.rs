/// import our big prime number list (10K) primes
mod primes;

/// based on https://web.archive.org/web/20120306040058/http://medialab.freaknet.org/martin/src/sqrt/sqrt.c
pub fn integer_square_root(n:i64) -> i64 {
    let mut op= n;
    let mut res = 0;

    // ideally do this with a bitscanreverse intrinsic...
    let mut bit = 1 << 62; // The second-to-top bit
    while bit > n {
        bit >>= 2;
    }

    // base 2 implementation of pen-and-paper square root (long division variant)
    while bit != 0 {
        if op >= res + bit {
            op -= res + bit;
            res += bit << 1;
        }
        res >>= 1;
        bit >>= 2;
    }
    res
}

/// returns true if n is a prime number. Not particularly fast.
pub fn is_prime(n:i64) -> bool {

    // some early outs
    if n == 1 { return true; }
    if n == 2 { return true; }
    if n % 2 == 0 { return false; }
    if n % 3 == 0 { return false; }

    if n > primes::K_PRIMES_RAW[primes::K_PRIMES_RAW.len()-1] as i64 {
        // n is bigger than the largest number in our prime list

        let max = integer_square_root(n);
        let mut p:i64=5;
        // brute force check if factorizable by any member of our prime list (skip past 2 and 3)
        for i in 2..primes::K_PRIMES_RAW.len() {
           if primes::K_PRIMES_RAW[i] as i64 > max {
               // no prime factor found
               return true;
           }
            if n % primes::K_PRIMES_RAW[i] as i64 == 0 {
                return false;
            }
            p = primes::K_PRIMES_RAW[i] as i64;
        }

        //TODO: can't really test any of this with 32 bit integers
        let mut di = 2;
        if (primes::K_PRIMES_RAW.len() & 1) != 0 {
            di = 4;
        }

        // straight brute force
        while p < max {
            if  n % p == 0  { return false;	}
            p += di;   // non multiples of 2 and 3 only {5,7,11,13,17,...}
            di = 6-di; // di = {2,4}
        }
    }
    else {
        // n is within the range of our prime list
        // binary search for a match
        let mut lo = 0;
        let mut hi = primes::K_PRIMES_RAW.len() - 1;
        while lo <= hi {
            let pivot = lo + (hi - lo) / 2;
            if primes::K_PRIMES_RAW[pivot] as i64 == n {
                return true;
            }
            if primes::K_PRIMES_RAW[pivot] as i64 > n {
                hi = pivot - 1;
            } else {
                lo = pivot + 1;
            }
        }
    }
    false
}

/// calculate a modulo m
/// correct according to mathematical convention (i.e. for a negative).
pub fn modulo(a:i64, m:i64) -> i64 {
    let r = a % m;
    // offset by m to make positive (iff a%m is negative)
    r + ((((r as u32) & 0x80000000)>>31) as i64)*m
}

/// extended Greatest Common Divisor algorithm
/// returns (s, t, GCD) such that a s + tb = GCD (s and t are the Bezout coefficients)
pub fn extended_gcd(a:i64, b:i64) -> (i64, i64, i64) {
    let mut s = 0;
    let mut t = 1;
    let mut r = b;
    let mut prev_s = 1;
    let mut prev_t = 0;
    let mut prev_r = a;
    while r != 0 {
        let quotient = prev_r / r;
        let tr = r;
        r = prev_r - quotient * tr;
        prev_r = tr;

        let ts = s;
        s = prev_s - quotient * ts;
        prev_s = ts;

        let tt = t;
        t = prev_t - quotient * tt;
        prev_t = tt;
    }

    // a s + t b = r
    (prev_s, prev_t, prev_r)
}

/// returns b^e (mod m)
/// not fast, but correct, and controls overflow by taking the modulo of each intermediate step.
pub fn power_mod(b:i64, e:u32, m:i64) -> Option<i64> {
    if m == 1 { return None; }
    if e == 0 { if m > 1 { return Some(1);} else { return None;} }

    let mut _c = 1;
    let mut _e_prime:u32 = 0;
    while  _e_prime != e {
        // this is an idiv every time, so don't expect performance!
        _c = (b * _c) % m;
        _e_prime = _e_prime+1;
    }
    Some(_c)
}

/// returns the inverse of x modulo m, iff x is relative prime to m, otherwise 0
/// solves the equation x*x_inv + k*m = 1, using the GCD
pub fn inv_modulo(x:i64, m:i64) -> i64 {
    let r = modulo(x,m );
    let (s,_,gcd) = extended_gcd(r, m);
    if gcd != 1 {
        // x and n have to be relative prime for there to exist an x^-1 mod n
        return 0;
    }
    modulo(s,m)
}

/// attempts to solve the congruence ax congruent b (mod m) iff a solution exists.
/// solves the equation as x congruent a_inv * b (mod m), finding a_inv using inv_modulo.
/// NOTE: this is not the most effective way to solve for large values of m, but it works
pub fn solve_linear_congruence(a:i64, b:i64, m:i64) -> Option<i64> {
    let a_inv = inv_modulo(a,m);
    if a_inv != 0 {
        return Some(modulo(b*a_inv,m));
    }
    None
}

/// ============================================================================
/// Unit tests
#[cfg(test)]
mod tests {
    use ::{modulo, extended_gcd,inv_modulo, is_prime};
    use ::{power_mod, solve_linear_congruence};
    use integer_square_root;

    #[test]
    fn test_modulo() {
        assert_eq!(6,modulo(-2,8));
        assert_eq!(6,modulo(6,11));
    }

    #[test]
    fn test_extended_gcd() {
        for a in 1..11 {
            let (s,t,gcd) = extended_gcd(a,11);
            assert_eq!(1,gcd);
            assert_eq!(a*s + 11*t,gcd);
        }
    }

    #[test]
    fn test_inv_modulo() {
        let mut inv_x = inv_modulo(4,7);
        assert_eq!(1,modulo(4 * inv_x, 7));
        assert_eq!(0,inv_modulo(4,8));
        inv_x = inv_modulo(-4,11);
        assert_eq!(1,modulo(-4 * inv_x, 11));
        inv_x = inv_modulo(-82,255);
        assert_eq!(1,modulo(-82 * inv_x, 255));
    }

    #[test]
    fn test_power_modulo() {
        assert_eq!(Some(6),power_mod(2,5,13));
        assert_eq!(Some(3),power_mod(101,100,13));
    }

    #[test]
    fn test_solve_linear_congruence() {
        assert_eq!(None, solve_linear_congruence(8,3,4));
        assert_eq!(Some(3), solve_linear_congruence(8,3,7));
        assert_eq!(3, solve_linear_congruence(100,13,7).unwrap());
    }

    #[test]
    fn test_is_prime() {
        assert_eq!(true, is_prime(2));
        assert_eq!(true, is_prime(104729));
        assert_eq!(true, is_prime(45691));
        assert_eq!(false, is_prime(8));
        assert_eq!(false, is_prime(104724));
        assert_eq!(false, is_prime(3*104729));
        assert_eq!(false, is_prime(104729*104729));
        assert_eq!(false, is_prime(104730*104730));
        assert_eq!(true, is_prime(10970467603));
        assert_eq!(true, is_prime(314173));
    }

    #[test]
    fn test_integer_squareroot() {
        assert_eq!(4,integer_square_root(16));
        assert_eq!(4,integer_square_root(17));
        assert_eq!(4,integer_square_root(18));
        assert_eq!(5,integer_square_root(25));
    }
}
