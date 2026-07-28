from sys import prefix

from mypy.checker import conditional_types

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
    priority=0,
    suffix=None,
    prefix=None,
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

ALL_PATTERNS = [CARBOXYLIC_ACID, SULFONIC_ACID, NITRILE, ALCOHOL]

# --- Carbon Chain Length Pattern ---

PREFIXES = {
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

STEMS = {
        "alkane": "an",
        "alkene": "en",
        "alkyne": "yn",
    }

HETEROATOM_PRIORITY = {

}
