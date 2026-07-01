"""
Contains the motif engine to match MotifPatterns to a Molecule.
"""
from typing import Optional, Any
from entities import MotifPattern, AtomCond, Atom, Molecule


class MotifMatch:
    """Result of a successful pattern match under a MotifPattern for a molecule.

    Instance Attributes:
        - pattern: The MotifPattern that was matched.
        - center_idx: The index of the center of the MotifPattern in the molecule.
        - matched_atoms: A list of indices for all the matched atoms for the MotifPattern.
    """
    pattern: MotifPattern
    matched_atoms: list[int]
    center_idx: int = 0

    def __init__(self, pattern: MotifPattern, center_idx: int, matched_atoms: list[int]):
        self.pattern = pattern
        self.center_idx = center_idx
        self.matched_atoms = matched_atoms


class MotifEngine:
    """The Matching engine, takes a molecule and matches it against a specific pattern.

    Instance Attributes:
        - molecule: The molecule that we wish to find patterns in.
        - _matches: A private attribute, holds all the matches so far.

    Class Constants:
        - CONDITION_HANDLER: A mapping from a condition property to its type (function, attribute, etc.)
    """
    molecule: Molecule
    _matches: list[MotifMatch]

    CONDITION_HANDLER: dict[str, Any]

    def __init__(self, molecule: Molecule) -> None:
        self.molecule = molecule
        self._matches = []

        self.CONDITION_HANDLER = {"charge": lambda atom: atom.charge,
                                  "aromatic": lambda atom: atom.is_aromatic,
                                  "has_h": lambda atom: atom.get_hydrogen_count() >= 1,
                                  "bond_sum": lambda atom: atom.bond_sum()}

    def match_all(self, patterns: list[MotifPattern]) -> list[MotifMatch]:
        """Returns a list of MotifMatches, detailing the matches found within the molecule of the input patterns.

        Parameters:
            - pattern: list of MotifPatterns.
                A list containing evey pattern to match against the molecule. All patterns are assumed to be valid.
        """
        all_matches = []

        for pattern in patterns:
            matches = self.match_pattern(pattern)
            all_matches.extend(matches)

        all_matches.sort(key=lambda m: m.pattern.priority, reverse=True)

        return self.resolve_overlaps(all_matches)


    def match_pattern(self, pattern: MotifPattern) -> list[MotifMatch]:
        """Returns a list of MotifMatches, detailing the matches found within the molecule of the input pattern.

        Parameters:
            - pattern: MotifPattern.
                The exact pattern to match against the molecule. Assumed to be a valid pattern.
        """
        matches = []

        candidates = self._find_candidates(pattern)

        for atom in candidates:
            match = self._match_single(atom, pattern)
            if match:
                match.center_idx = atom.idx
                matches.append(match)

        self._matches = matches
        return matches

    def _find_candidates(self, pattern: MotifPattern) -> list[Atom]:
        """Returns a list of atoms which are candidates for a central atom to the pattern.

        Parameters:
            - pattern: MotifPattern
                The exact pattern to match against the molecule. Asssumed to be a valid pattern.
        """
        candidates = []

        for atom in self.molecule:
            if atom.element.symbol != pattern.center_symbol:
                continue

            if all(self._condition_single(atom, cond) for cond in pattern.center_conditions):
                candidates.append(atom)

        return candidates

    def _match_single(self, atom: Atom, pattern: MotifPattern, visited: set[int] = None) -> Optional[MotifMatch]:
        """Returns a MotifMatch if there is a match for the MotifPattern from this atom, None otherwise.

                Parameters:
                    - atom: Atom
                        The atom being matched to some portion of the MotifPattern.
                    - pattern: MotifPattern
                        The pattern that is trying to be matched within the molecule.
                    - visited: set of integers
                        The set of visitied atom indices in order to prevent infinited regress.
        """
        if not visited:
            visited = set()

        if atom.idx in visited:
            return None

        visited.add(atom.idx)
        used = set()  # Tracks consumed neighbours
        matches = []  # Tracks any matched neighbours

        for req in pattern.bonds:
            found_count = 0

            for neighbour_idx, order in atom.bonds.items():
                if neighbour_idx in used:
                    continue  # Disregard consumed neighbours.

                neighbour = self.molecule[neighbour_idx]
                if neighbour.element.symbol == req.symbol and order == req.order:
                    if all(self._condition_single(neighbour, cond) for cond in req.conditions):
                        # If it passes all atom condititions under that role, then we have a match!
                        used.add(neighbour_idx)
                        matches.append(neighbour_idx)
                        found_count += 1

                        # Checks for any further bonds that need to be made.
                        if req.future_req:
                            future_match = self._match_single(neighbour, pattern, visited)
                            if not future_match:
                                return None

                            # Merge the matches
                            matches.extend(idx for idx in future_match.matched_atoms if idx not in matches)

            if found_count != req.count:
                return None

        return MotifMatch(pattern=pattern, matched_atoms=matches, center_idx=0)

    def _condition_single(self, atom: Atom, cond: AtomCond) -> bool:
        """Returns if the atom follows the input condition.

        Parameters:
            - atom: Atom
                The atom that is being checked.
            - cond: AtomCond
                The condition the atom is being checked against.
        """
        prop = self.CONDITION_HANDLER[cond.property](atom)

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

    @staticmethod
    def resolve_overlaps(matches: list[MotifMatch]) -> list[MotifMatch]:
        """Returns a list of MotifMatches with all overlaping matches resolved. Resolution is done by removing the match
        from whose pattern has a lower priority.

        Parameters:
            - matches: list of MotifMatchs
                The list of MotifMatchs to be resolved. The list must be sorted descending based on the priority of the
                match's pattern.
        """
        used_atoms = set()
        resolved = []

        for match in matches:
            if any(idx in used_atoms for idx in match.matched_atoms):
                continue

            used_atoms.update(match.matched_atoms)
            resolved.append(match)

        return resolved
