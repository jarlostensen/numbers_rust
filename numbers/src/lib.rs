
/// calculate a modulo m
/// correct according to mathematical convention.
pub fn modulo(a:i32, m:i32) -> i32 {
    let r = a % m;
    r + ((((r as u32) & 0x80000000)>>31) as i32)*m
}


#[cfg(test)]
mod tests {
    use modulo;

    #[test]
    fn test_modulo() {
        assert_eq!(6,modulo(-2,8));
        assert_eq!(6,modulo(6,11));
    }
}
