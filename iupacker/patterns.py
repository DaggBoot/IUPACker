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
    detached_suffix="carboxylic acid",
    prefix="carboxy-",
    inline=False,
    always_terminal=True,
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
    detached_suffix=None,
    prefix="sulfo-",
    inline=False,
    always_terminal=False,
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
    detached_suffix=None,
    prefix="oxy",
    inline=True,
    always_terminal=False,
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

ESTER = MotifPattern(
    name="ester",
    priority=75,
    suffix="oate",
    detached_suffix=None,
    prefix=None,
    inline=True,
    always_terminal=False,
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

ACYL_F = MotifPattern(
    name="acyl fluoride",
    priority=70,
    suffix="oyl fluoride",
    detached_suffix="carbonyl fluoride",
    prefix="fluorocarbonyl",
    inline=False,
    always_terminal=True,
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
            symbol="F",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_BOND_SUM_1
            ]
        )
    ]
)

ACYL_CL = MotifPattern(
    name="acyl chloride",
    priority=70,
    suffix="oyl chloride",
    detached_suffix="carbonyl chloride",
    prefix="chlorocarbonyl",
    inline=False,
    always_terminal=True,
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
            symbol="Cl",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_BOND_SUM_1
            ]
        )
    ]
)

ACYL_BR = MotifPattern(
    name="acyl bromide",
    priority=70,
    suffix="oyl bromide",
    detached_suffix="carbonyl bromide",
    prefix="bromocarbonyl",
    inline=False,
    always_terminal=True,
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
            symbol="Br",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_BOND_SUM_1
            ]
        )
    ]
)

ACYL_I = MotifPattern(
    name="acyl iodide",
    priority=70,
    suffix="oyl iodide",
    detached_suffix="carbonyl iodide",
    prefix="iodocarbonyl",
    inline=False,
    always_terminal=True,
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
            symbol="I",
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
    detached_suffix="carboxamide",
    prefix="carbamoyl",
    inline=True,
    always_terminal=False,
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
    suffix="nitrile",
    detached_suffix="carbonitrile",
    prefix="cyano",
    inline=False,
    always_terminal=True,
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
    detached_suffix="carbaldehyde",
    prefix="formyl",
    inline=False,
    always_terminal=True,
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
    suffix="one",
    detached_suffix=None,
    prefix="oxo",
    inline=False,
    always_terminal=False,
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
        )
    ]
)

ALCOHOL = MotifPattern(
    name="alcohol",
    priority=10,
    suffix="ol",
    detached_suffix=None,
    prefix="hydroxy",
    inline=False,
    always_terminal=False,
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
    detached_suffix=None,
    prefix="amino",
    inline=True,
    always_terminal=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="N",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
                COND_HAS_H
            ],
        )
    ]
)

FLUORIDE = MotifPattern(
    name="fluoride",
    priority=0,
    suffix="",
    detached_suffix=None,
    prefix="fluoro",
    inline=False,
    always_terminal=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="F",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
            ],
        )
    ]
)

CHLORIDE = MotifPattern(
    name="chloride",
    priority=0,
    suffix="",
    detached_suffix=None,
    prefix="chloro",
    inline=False,
    always_terminal=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="Cl",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
            ],
        )
    ]
)

BROMIDE = MotifPattern(
    name="bromide",
    priority=0,
    suffix="",
    detached_suffix=None,
    prefix="bromo",
    inline=False,
    always_terminal=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="Br",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
            ],
        )
    ]
)

IODIDE = MotifPattern(
    name="iodide",
    priority=0,
    suffix="",
    detached_suffix=None,
    prefix="iodo",
    inline=False,
    always_terminal=False,
    center_symbol="C",
    center_conditions=[
        COND_NO_CHARGE,
        COND_AROMATIC_FALSE
    ],
    bonds=[
        BondReq(
            symbol="I",
            order=1,
            count=1,
            conditions=[
                COND_NO_CHARGE,
                COND_AROMATIC_FALSE,
            ],
        )
    ]
)

ALL_PATTERNS = [CARBOXYLIC_ACID, SULFONIC_ACID, ANHYDRIDE, ESTER, ACYL_F, ACYL_CL, ACYL_BR, ACYL_I,
                AMIDE, NITRILE, ALDEHYDE, KETONE, ALCOHOL, AMINE, FLUORIDE, CHLORIDE, BROMIDE, IODIDE]

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

COMPLEX_MULT_PREFIXES = {
    1: "",
    2: "bis",
    3: "tris",
    4: "tetrakis",
    5: "pentakis",
    6: "hexakis",
    7: "heptakis",
    8: "octakis",
    9: "nonakis",
    10: "decakis",
    11: "undecakis",
    12: "dodecakis",
    13: "tridecakis",
    14: "tetradecakis",
    15: "pentadecakis",
    16: "hexadecakis",
    17: "heptadecakis",
    18: "octadecakis",
    19: "nonadecakis",
    20: "icosakis",
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
