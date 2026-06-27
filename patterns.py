from entities import MotifPattern, BondReq, AtomCond

# --- Reusable Condition Constants ---

COND_NO_CHARGE = [AtomCond("charge", 0)]
COND_HAS_H = [AtomCond("has_h", True)]
COND_NO_H = [AtomCond("has_h", False)]
COND_BOND_SUM_2 = [AtomCond("bond_sum", 2)]
COND_AROMATIC_FALSE = [AtomCond("aromatic", False)]

# --- Carboxylic Acid Pattern ---

CARBOXYLIC_ACID = MotifPattern(
    name="carboxylic_acid",
    priority=100,
    suffix="oic acid",
    prefix="carboxy-",
    center_symbol="C",
    center_conditions=COND_NO_CHARGE + COND_AROMATIC_FALSE,  # Must be neutral, not aromatic
    bonds=[
        BondReq(
            symbol="O",
            order=2,
            count=1,
            conditions=COND_NO_H + COND_BOND_SUM_2,  # Carbonyl O: no H, exactly 2 bonds
        ),
        BondReq(
            symbol="O",
            order=1,
            count=1,
            conditions=COND_HAS_H + COND_BOND_SUM_2,  # Hydroxyl O: has H, exactly 2 bonds
        ),
    ],
    excludes=[],  # No exclusions needed for carboxylic acid
)
