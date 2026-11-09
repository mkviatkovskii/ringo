#[derive(Hash, Eq, PartialEq, Debug)]
pub struct Element {
    pub atomic_number: u8,
}

impl Element {
    #[expect(
        clippy::unreadable_literal,
        clippy::too_many_lines,
        reason = "a table of weights, written as in the IUPAC table"
    )]
    pub fn atomic_weight(&self) -> f64 {
        match self.atomic_number {
            1 => 1.008,
            2 => 4.002602,
            3 => 6.94,
            4 => 9.0121831,
            5 => 10.81,
            6 => 12.011,
            7 => 14.007,
            8 => 15.999,
            9 => 18.998403163,
            10 => 20.180,
            11 => 22.990,
            12 => 24.305,
            13 => 26.9815385,
            14 => 28.085,
            15 => 30.973761998,
            16 => 32.06,
            17 => 35.45,
            18 => 39.948,
            19 => 39.0983,
            20 => 40.078,
            21 => 44.955908,
            22 => 47.867,
            23 => 50.9415,
            24 => 51.9961,
            25 => 54.938044,
            26 => 55.845,
            27 => 58.933194,
            28 => 58.6934,
            29 => 63.546,
            30 => 65.38,
            31 => 69.723,
            32 => 72.630,
            33 => 74.921595,
            34 => 78.971,
            35 => 79.904,
            36 => 83.798,
            37 => 85.4678,
            38 => 87.62,
            39 => 88.90584,
            40 => 91.224,
            41 => 92.90637,
            42 => 95.95,
            43 => 98.0,
            44 => 101.07,
            45 => 102.90550,
            46 => 106.42,
            47 => 107.8682,
            48 => 112.414,
            49 => 114.818,
            50 => 118.710,
            51 => 121.760,
            52 => 127.60,
            53 => 126.90447,
            54 => 131.293,
            55 => 132.90545196,
            56 => 137.327,
            57 => 138.90547,
            58 => 140.116,
            59 => 140.90766,
            60 => 144.242,
            61 => 145.0,
            62 => 150.36,
            63 => 151.964,
            64 => 157.25,
            65 => 158.92535,
            66 => 162.500,
            67 => 164.93033,
            68 => 167.259,
            69 => 168.93422,
            70 => 173.054,
            71 => 174.9668,
            72 => 178.49,
            73 => 180.94788,
            74 => 183.84,
            75 => 186.207,
            76 => 190.23,
            77 => 192.217,
            78 => 195.084,
            79 => 196.966569,
            80 => 200.592,
            81 => 204.38,
            82 => 207.2,
            83 => 208.98040,
            // no stable isotope: mass number of the longest-lived one, as for Tc and Pm
            84 => 209.0,
            85 => 210.0,
            86 => 222.0,
            87 => 223.0,
            88 => 226.0,
            89 => 227.0,
            90 => 232.0377,
            91 => 231.03588,
            92 => 238.02891,
            93 => 237.0,
            94 => 244.0,
            95 => 243.0,
            96 | 97 => 247.0,
            98 => 251.0,
            99 => 252.0,
            100 => 257.0,
            101 => 258.0,
            102 => 259.0,
            103 => 266.0,
            104 => 267.0,
            105 => 268.0,
            106 | 108 => 269.0,
            107 => 270.0,
            109 => 278.0,
            110 => 281.0,
            111 => 282.0,
            112 => 285.0,
            113 => 286.0,
            114 => 289.0,
            115 => 290.0,
            116 => 293.0,
            117 | 118 => 294.0,
            // the `*` wildcard atom
            0 => 0.0,
            _ => panic!("Unsupported atomic number {}", self.atomic_number),
        }
    }

    /// The symbol, such as "C" or "Cl"; "*" for the wildcard.
    pub fn symbol(&self) -> &'static str {
        SYMBOLS[usize::from(self.atomic_number)]
    }

    /// The element for a symbol such as "C" or "Cl"; "*" is the wildcard, atomic number 0.
    pub fn from_symbol(symbol: &str) -> Option<Element> {
        let index = SYMBOLS.iter().position(|&s| s == symbol)?;
        u8::try_from(index)
            .ok()
            .map(|atomic_number| Element { atomic_number })
    }
}

/// Element symbols by atomic number; index 0 is the wildcard.
const SYMBOLS: [&str; 119] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

#[cfg(test)]
mod test {
    use crate::model::element::Element;

    #[test]
    fn test_element() {
        let element = Element { atomic_number: 1 };
        assert_eq!(element.atomic_weight(), 1.008);
    }

    #[test]
    fn test_every_element_has_a_symbol_and_weight() {
        for (symbol, atomic_number) in [("*", 0), ("C", 6), ("Cl", 17), ("U", 92), ("Og", 118)] {
            assert_eq!(
                Element::from_symbol(symbol),
                Some(Element { atomic_number }),
                "{symbol}"
            );
        }
        assert_eq!(Element::from_symbol("Xx"), None);
        assert_eq!(Element::from_symbol("c"), None);
        // used to panic above bismuth (83)
        for atomic_number in 0..=118 {
            let weight = Element { atomic_number }.atomic_weight();
            assert!(atomic_number == 0 || weight > 1.0, "{atomic_number}");
        }
    }
}
