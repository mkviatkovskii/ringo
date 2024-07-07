use nom::{
    branch::alt,
    character::complete::{char, digit1},
    combinator::map_res,
    multi::many1_count,
    IResult,
};

/// Parses a charge, which could be:
/// * single `+` or `-` sign => return 1 or -1
/// * multiple `+` or `-` signs => return n or -n
/// * `+` or `-` sign followed by a number => return n or -n
pub fn parse_charge(input: &str) -> IResult<&str, i8> {
    let (input, sign) = alt((char('+'), char('-')))(input)?;
    let charge = match sign {
        '+' => 1,
        '-' => -1,
        _ => unreachable!("charge is either + or -"),
    };
    let sign_result = many1_count(char(sign))(input);
    if sign_result.is_ok() {
        let (input, sign_count) = sign_result?;
        return Ok((input, charge * (sign_count + 1) as i8));
    }
    let digit_result = map_res(digit1, str::parse::<u8>)(input);
    return if digit_result.is_ok() {
        let (input, count) = digit_result?;
        Ok((input, charge * count as i8))
    } else {
        Ok((input, charge))
    };
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::charge::parse_charge;
    use nom::error::{Error, ErrorKind};

    #[test]
    fn parse_charge_empty() {
        assert_eq!(
            parse_charge(""),
            Err(nom::Err::Error(Error::new("", ErrorKind::Char)))
        );
    }

    #[test]
    fn parse_charge_a() {
        assert_eq!(
            parse_charge("a"),
            Err(nom::Err::Error(Error::new("a", ErrorKind::Char)))
        );
    }

    #[test]
    fn parse_charge_plus() {
        assert_eq!(parse_charge("+").unwrap().1, 1);
    }

    #[test]
    fn parse_charge_minus() {
        assert_eq!(parse_charge("-").unwrap().1, -1);
    }

    #[test]
    fn parse_charge_plus_2() {
        assert_eq!(parse_charge("+2").unwrap().1, 2);
    }

    #[test]
    fn parse_charge_minus_2() {
        assert_eq!(parse_charge("-2").unwrap().1, -2);
    }

    #[test]
    fn parse_charge_plus_plus() {
        assert_eq!(parse_charge("++").unwrap().1, 2);
    }

    #[test]
    fn parse_charge_minus_minus() {
        assert_eq!(parse_charge("--").unwrap().1, -2);
    }
}
