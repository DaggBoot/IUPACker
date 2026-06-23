"""
Contains the motif engine to match MotifPatterns to a Molecule.
"""
from typing import Optional
from entities import MotifPattern, BondReq, AtomCond, Atom, Molecule


class MotifMatch:
    """Result of a successful pattern match under a MotifPattern for a molecule.

    Instance Attributes:
        - pattern: The MotifPattern that was matched.
        - center_idx: The index of the center of the MotifPattern in the molecule.
        - matched_atoms: A list of indices for all the matched atoms for the MotifPattern.
        - atom_role: A mapping that takes the chemical role and maps it to the index of the atom under said role.
    """
    pattern: MotifPattern
    center_idx: int
    matched_atoms: list[int]
    atom_role: dict[str, int]  # chemical role (eg: carbonyl) -> atom_idx

    def __init__(self, pattern: MotifPattern, center_idx: int, matched_atoms: list[int]):
        self.pattern = pattern
        self.center_idx = center_idx
        self.matched_atoms = matched_atoms
        self.atom_roles = {}


class MotifEngine:
    """The Matching engine, takes a molecule and matches it against a specific pattern.

    Instance Attributes:
        - molecule: The molecule that we wish to find patterns in.
        - _matches: A private attribute, holds all the matches so far.
    """
    molecule: Molecule
    _matches: list[MotifMatch]

    def __init__(self, molecule: Molecule) -> None:
        self.molecule = molecule
        self._matches = []

    def match_pattern(self, pattern: MotifPattern) -> list[MotifMatch]:
        """Returns a list of MotifMatches, detailing the matches found of the input pattern in the molecule stored.

        Parameters:
            - pattern: MotifPattern.
                The exact pattern to match against the molecule. Asssumed to be a valid pattern.
        """
        matches = []

        candidates = self._find_candidates(pattern)

        for atom in candidates:
            match = self._match_single(atom, pattern)
            if match:
                matches.append(match)

        return matches

    def _find_candidates(self, pattern: MotifPattern) -> list[Atom]:
        """Returns a list of atoms which are candidates for a central atom to the pattern."""
        candidates = []

        for atom in self.molecule:
            if atom.element.symbol != pattern.center_symbol:
                continue

            for cond in pattern.center_conditions:
                if self._condition_single(atom, cond):
                    candidates.append(atom)

        return candidates

    def _match_single(self, atom: Atom, pattern: MotifPattern) -> Optional[MotifMatch]:
        """Returns a MotifMatch if there is a match for the MotifPattern from this atom, None otherwise.

        Parameters:
            - atom: Atom
                An Atom, assumed to be the central atom in the MotifPattern
            - pattern: MotifPattern
                The pattern that is trying to be matched within the molecule.
        """
        bond_matches = self._match_bonds(atom, pattern.bonds)

    def _match_bonds(self, atom: Atom, bonds: list[BondReq]) -> Optional[list[tuple[int, float]]]:

    @staticmethod
    def _condition_single(atom: Atom, cond: AtomCond) -> bool:
        """Returns if the atom follows the input condition.

        Parameters:
            - atom: Atom
                The atom that is being checked.
            - cond: AtomCond
                The condition the atom is being checked against.
        """
        prop = getattr(atom, cond.property, None)

        if cond.oper == "eq":
            return prop == cond.value
        elif cond.oper == "ne":
            return prop != cond.value
        elif cond.oper == "lt":
            return prop < cond.value
        elif cond.oper == "gt":
            return prop > cond.value
        elif cond.oper == "le":
            return prop <= cond.value
        elif cond.oper == "ge":
            return prop >= cond.value
        else:
            return False
