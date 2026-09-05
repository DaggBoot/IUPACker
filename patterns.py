from entities import MotifPattern, BondReq, AtomCond

# --- Reusable Condition Constants ---

COND_NO_CHARGE = AtomCond("charge", 0)
COND_HAS_H = AtomCond("has_h", True)
COND_NO_H = AtomCond("has_h", False)
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

NITRILE = MotifPattern(
    name="nitrile",
    priority=None,
    suffix=None,
    prefix=None,
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

ALL_PATTERNS = [CARBOXYLIC_ACID, SULFONIC_ACID, NITRILE, ALCOHOL, ANHYDRIDE]

# --- Carbon Chain Length Pattern ---

ALKYL_PREFIXES = {
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
}

SATS = {
        1: "an",
        2: "en",
        3: "yn",
    }

HETEROATOM_PRIORITY = {

}
