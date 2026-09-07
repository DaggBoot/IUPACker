from .entities import MotifPattern, BondReq, AtomCond

# --- Reusable Condition Constants ---

COND_NO_CHARGE = AtomCond("charge", 0)
COND_HAS_H = AtomCond("has_h", True)
COND_NO_H = AtomCond("has_h", False)
COND_BOND_SUM_1 = AtomCond("bond_sum", 1)
COND_BOND_SUM_2 = AtomCond("bond_sum", 2)
COND_AROMATIC_FALSE = AtomCond("aromatic", False)

# --- Functional Group Patterns ---

CARBOXYLIC_ACID = MotifPattern(
    name="carboxylic acid",
    priority=100,
    suffix="oic acid",
    prefix="carboxy-",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
        ),
        BondReq(
            symbol="O",
            order=1,
            count=1,
            conditions=[
                COND_HAS_H,
                COND_BOND_SUM_2
            ],
        ),
    ],
    excludes=[],
)

SULFONIC_ACID = MotifPattern(
    name="sulfonic acid",
    priority=90,
    suffix="sulfonic acid",
    prefix="sulfo-",
    inline=False,
    center_symbol="S",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE,
        COND_NO_H
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=2,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_NO_H
            ],
        ),
        BondReq(
            symbol="O",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_HAS_H
            ],
        )
    ],
    excludes=[],
)


ANHYDRIDE = MotifPattern(
    name="anhydride",
    priority=80,
    suffix="anoic anhydride",
    prefix="oxy",
    inline=True,
    center_symbol="O",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE,
    ],
    bonds=[
        BondReq(
            symbol="C",
            order=1,
            count=2,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_NO_H
            ],
            future_req=[
                BondReq(
                    symbol="O",
                    order=2,
                    count=1,
                    conditions=[
                        COND_NO_H,
                        COND_BOND_SUM_2
                    ]
                )
            ]
        )
    ]
)

ESTERS = MotifPattern(
    name="ester",
    priority=75,
    suffix="oate",
    prefix=None,
    inline=True,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE,
        COND_NO_H
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
        ),
        BondReq(
            symbol="O",
            order=1,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
            future_req=[
                BondReq(
                    symbol="C",
                    order=1,
                    count=1,
                    conditions=[
                        COND_NO_CHARGE,
                        COND_AROMATIC_FALSE,
                    ]
                )
            ]
        ),
    ]
)

ACYL_HALIDE = MotifPattern(
    name="acyl halide",
    priority=70,
    suffix="oyl",
    prefix="halocarbonyl",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
        ),
        BondReq(
            symbol="FClIBr",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_BOND_SUM_1
            ]
        )
    ]
)

AMIDE = MotifPattern(
    name="amide",
    priority=65,
    suffix="amide",
    prefix="carbamoyl",
    inline=True,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
        ),
        BondReq(
            symbol="N",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE
            ]
        )
    ]
)

NITRILE = MotifPattern(
    name="nitrile",
    priority=60,
    suffix="cyano",
    prefix="nitrile",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="N",
            order=3,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_NO_H
            ]
        )
    ]
)

ALDEHYDE = MotifPattern(
    name="aldehyde",
    priority=55,
    suffix="al",
    prefix="formyl",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE,
        COND_HAS_H
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ],
        )
    ]
)

KETONE = MotifPattern(
    name="ketone",
    priority=55,
    suffix="oxo",
    prefix="one",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE,
        COND_NO_H
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=[
                COND_NO_H,
                COND_BOND_SUM_2
            ]
        ),
        BondReq(
            symbol="C",
            order=1,
            count=2,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE
            ]
        )
    ]
)

ALCOHOL = MotifPattern(
    name="alcohol",
    priority=10,
    suffix="ol",
    prefix="hydroxy",
    inline=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="O",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_HAS_H
            ]
        )
    ]
)

AMINE = MotifPattern(
    name="amine",
    priority=5,
    suffix="amine",
    prefix="amino",
    inline=True,
    center_symbol="N",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ]
)

ALL_PATTERNS = [CARBOXYLIC_ACID, SULFONIC_ACID, ANHYDRIDE, ESTERS, ACYL_HALIDE,
                AMIDE, NITRILE, ALDEHYDE, KETONE, ALCOHOL, AMINE]

# --- Carbon Chain Length Pattern ---

MULT_PREFIXES = {
    0: "",
    1: "",
    2: "di",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
    7: "hepta",
    8: "octa",
    9: "nona",
    10: "deca",
}

SIMPLE_PREFIXES = {
    1: "meth",
    2: "eth",
    3: "prop",
    4: "but",
    5: "pent",
    6: "hex",
    7: "hept",
    8: "oct",
    9: "non",
    10: "dec",
    11: "undec",
    12: "dodec",
    13: "tridec",
    14: "tetradec",
    15: "pentadec",
    16: "hexadec",
    17: "heptadec",
    18: "octadec",
    19: "nonadec",
    20: "eicos",
}

ONES_CONSTRUCTOR_PREFIXES = {
    1: "hen",
    2: "do",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
    7: "hepta",
    8: "octa",
    9: "nona",
}

TENS_CONSTRUCTOR_PREFIXES = {
    1: "dec",
    2: "icos",
    3: "triacont",
    4: "tetracont",
    5: "pentacont",
    6: "hexacont",
    7: "heptacont",
    8: "octacont",
    9: "nonacont",
}

SATS = {
        1: "an",
        2: "en",
        3: "yn",
    }

HETEROATOM_PRIORITY = {

}
