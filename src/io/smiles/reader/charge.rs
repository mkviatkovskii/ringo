use winnow::ascii::digit1;
use winnow::combinator::{opt, repeat};
use winnow::prelude::*;
use winnow::token::one_of;

/// Parses a charge, which could be:
/// * single `+` or `-` sign => return 1 or -1
/// * multiple `+` or `-` signs => return n or -n
/// * `+` or `-` sign followed by a number => return n or -n
pub fn parse_charge(input: &mut &str) -> ModalResult<i8> {
    let sign = one_of(['+', '-']).parse_next(input)?;
    let charge: i8 = if sign == '+' { 1 } else { -1 };
    let repeated: usize = repeat(0.., sign).parse_next(input)?;
    if repeated > 0 {
        return Ok(charge * (repeated + 1) as i8);
    }
    let count = opt(digit1.try_map(str::parse::<u8>)).parse_next(input)?;
    Ok(charge * count.map_or(1, |count| count as i8))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::charge::parse_charge;
    use winnow::Parser;

    #[test]
    fn parse_charge_empty() {
        assert!(parse_charge.parse_peek("").is_err());
    }

    #[test]
    fn parse_charge_a() {
        assert!(parse_charge.parse_peek("a").is_err());
    }

    #[test]
    fn parse_charge_plus() {
        assert_eq!(parse_charge.parse_peek("+").unwrap().1, 1);
    }

    #[test]
    fn parse_charge_minus() {
        assert_eq!(parse_charge.parse_peek("-").unwrap().1, -1);
    }

    #[test]
    fn parse_charge_plus_2() {
        assert_eq!(parse_charge.parse_peek("+2").unwrap().1, 2);
    }

    #[test]
    fn parse_charge_minus_2() {
        assert_eq!(parse_charge.parse_peek("-2").unwrap().1, -2);
    }

    #[test]
    fn parse_charge_plus_plus() {
        assert_eq!(parse_charge.parse_peek("++").unwrap().1, 2);
    }

    #[test]
    fn parse_charge_minus_minus() {
        assert_eq!(parse_charge.parse_peek("--").unwrap().1, -2);
    }
}
