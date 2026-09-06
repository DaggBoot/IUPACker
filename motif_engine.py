"""
Contains the motif engine to match MotifPatterns to a Molecule.
"""
from typing import Optional, Any
from entities import MotifPattern, AtomCond, Atom, Molecule
import itertools

# FOR TESTING ONLY --------------------
from smiles_parser import SMILESParser
import patterns
# -------------------------------------


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

        return tuple(self.resolve_overlaps(all_matches))

    def match_pattern(self, pattern: MotifPattern) -> list[MotifMatch]:
        """Returns a list of MotifMatches, detailing the matches found within the molecule of the input pattern.

        Parameters:
            - pattern: MotifPattern.
                The exact pattern to match against the molecule. Assumed to be a valid pattern.
        """
        matches = []

        candidates = self._find_candidates(pattern)

        for atom in candidates:
            atom_matches = self._match_single(atom, pattern)
            for match in atom_matches:
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

    def _match_single(self, atom: Atom, pattern: MotifPattern, visited: set[int] = None) -> Optional[tuple[MotifMatch]]:
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
            return []

        visited.add(atom.idx)

        # First, you gotta find all atoms that qualify as parts of a functional group.
        req_options = []
        for req in pattern.bonds:
            qualifying = [
                neighbour_idx for neighbour_idx, order in atom.bonds.items()
                if self.molecule[neighbour_idx].element.symbol in req.symbol
                and order == req.order
                and all(self._condition_single(self.molecule[neighbour_idx], cond) for cond in req.conditions)
            ]

            # If you don't have enough for 1 bond req, then you feasibly cant find the func group
            if len(qualifying) < req.count:
                return []

            # Generates all combinations of the qualifying atoms (of group siz req.count) for a requirment
            req_options.append((req, list(itertools.combinations(qualifying, req.count))))

        # Now, we try every single combination of choices for each bond requirment
        results = []
        for combo_choice in itertools.product(*(options for _, options in req_options)):
            chosen = set()
            overlap = False

            # Check if any atom is used in multiple bond requirements
            for atom in combo_choice:
                if chosen & set(atom):
                    overlap = True
                    break
                chosen.update(atom)

            if overlap:
                continue  # Skip this combination to prevent double-counting of atom

            matched_atoms = list(chosen)
            ok = True

            # Handle the future_reqs
            for (req, _), combo in zip(req_options, combo_choice):
                if req.future_req:
                    for neighbour_idx in combo:
                        temp_pattern = MotifPattern(
                            name="x",
                            priority=None,
                            inline=False,
                            prefix=None,
                            suffix=None,
                            center_symbol=req.symbol,
                            center_conditions=req.conditions,
                            bonds=req.future_req
                        )
                        future_matches = self._match_single(self.molecule[neighbour_idx], temp_pattern, visited.copy())

                        if not future_matches:
                            ok = False
                            break

                        # merge any matched atoms from future matches to avoid any duplicates
                        matched_atoms.extend(
                            idx for idx in future_matches[0].matched_atoms if idx not in matched_atoms
                        )
                    if not ok:
                        break

            if ok:
                results.append(MotifMatch(pattern=pattern, matched_atoms=matched_atoms, center_idx=0))

        return results

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
            if any(idx in used_atoms for idx in match.matched_atoms + [match.center_idx]):
                continue

            used_atoms.update(match.matched_atoms + [match.center_idx])
            resolved.append(match)

        return resolved


if __name__ == "__main__":
    smiles = "CCCCN(CCC)CCCC"
    parser = SMILESParser()
    molecule = parser.parse(smiles)
    engine = MotifEngine(molecule)
    groups = engine.match_all(patterns.ALL_PATTERNS)
    for match in groups:
        print(match.pattern.name)
