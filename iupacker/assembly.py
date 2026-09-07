"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from __future__ import annotations
from typing import Optional, NamedTuple, Union
from .motif_engine import MotifEngine, MotifMatch
from .entities import Molecule, Ring
from .smiles_parser import SMILESParser
from . import patterns


def generate_name(smiles: str) -> str:
    """Returns the chemical name of the inputted chemcial formula as per IUPAC nomenclature.

    Parameters:
        - smiles: string
            The molecule's chemical formula, assumed to be in the SMILES formula under the OpenSMILES documentation.
            Will raise exceptions in the case of improper SMILES formatting or chemcially invalid formulae.
    """
    parser = SMILESParser()
    molecule = parser.parse(smiles)
    namer = IUPACker(molecule)
    name = namer.generate()
    for match in namer._groups:
        print(match.pattern.name)

    return name
    print(name)


def generate_name_from_mol(molecule: Molecule):
    """Returns the chemical name of the inputted Molecule as per IUPAC nomenclature.

    Parameters:
        - mol: Molecule
            An assumed chemically valid Molecule Graph object.
    """
    namer = IUPACker(molecule)
    return namer.generate()


class _Candidate(NamedTuple):
    """TODO:"""
    chain: list[int]
    elem: str
    cyclic: Union[bool, Ring] = False


class _Substituent(NamedTuple):
    """TODO:"""
    root: int
    chain: tuple[int, ...]
    cyclic: Optional[Ring]
    groups: tuple[MotifMatch, ...]
    nested: tuple[_Substituent, ...]
    locant: int = None


class IUPACker:
    """Generates IUPAC names from molecular graphs

    This class handles the complete IUPAC naming process for organic molecules.
    It detects functional groups, identifies the principal characteristic group finds the parent chain, determines
    substituents, applies correct numbering, and assembles the final IUPAC name following IUPAC priority rules.

    The naming process follows these steps:
        1. Detect all functional groups in the molecule
        2. Identify the highest priority group (principal characteristic group)
        3. Find the longest chain containing the principal group (parent chain)
        4. Identify all substituents (branches and other functional groups)
        5. Number the parent chain to give lowest locants
        6. Assemble the final IUPAC name with appropriate suffixes and prefixes

    Instance Attributes:
        - _molecule: The Molecule object to be named
        - _groups: Immutable collection of all MotifMatch objects that correspond to all the functional groups present.
                   Arranged descending based on priority.
        - _parent_chain: The atom indices of the chosen parent chain
        - _subs: Immutable collection of all subsitients to the parent chain
    """
    _molecule: Molecule
    _groups: tuple[MotifMatch, ...]
    _parent_chain: tuple[int, ...]
    _subs: tuple[_Substituent, ...]

    _SENIOR_ELEMENTS = ("N", "P", "Si", "B", "O", "S", "C")

    def __init__(self, molecule: Molecule = None):
        self._molecule = molecule
        self._groups = ()
        self._parent_chain = ()
        self._subs = ()

    def generate(self) -> str:
        """Returns the chemical name of self._molecule as per IUPAC nomenclature of the Molecule.

        TODO: Rest of the stuff
        """

        engine = MotifEngine(self._molecule)
        self._groups = engine.match_all(patterns.ALL_PATTERNS)
        princip_groups = [group for group in self._groups if group.pattern == self._groups[0].pattern]

        if princip_groups:
            princip_candidate = self._find_parent_chain(princip_groups)
            print(princip_candidate)
        else:
            princip_candidate = self._find_parent_chain_no_p()

        self._parent_chain = tuple(self._number_parent(princip_candidate))
        self._subs = self._find_substituents()

        print((str(len(self._parent_chain)) + str(self._parent_chain) + str(len(self._subs)) + str(self._subs)))

        central = self._name_parent(princip_groups, isinstance(princip_candidate.cyclic, Ring))

        return central

    # Name construction

    def _name_parent(self, princip_groups: list[MotifMatch], cyclic: bool = False) -> str:
        name = ""

        if cyclic:
            name = "cyclo"

        try:
            name += patterns.SIMPLE_PREFIXES[len(self._parent_chain)]
        except KeyError:
            raise ValueError(f"Unsupported chain length: {len(self._parent_chain)}")

        doubles = []
        triples = []
        princips = []
        for i, idx in enumerate(self._parent_chain):
            atom = self._molecule[idx]

            if idx in {group.center_idx for group in princip_groups}:
                princips.append(i + 1)

            if i < len(self._parent_chain) - 1:
                if atom.bonds.get(self._parent_chain[i + 1], 1) == 2:
                    doubles.append(i + 1)

                elif atom.bonds.get(self._parent_chain[i + 1], 1) == 3:
                    triples.append(i + 1)

            elif cyclic:
                if atom.bonds.get(self._parent_chain[0], 1) == 2:
                    doubles.append(i + 1)

                elif atom.bonds.get(self._parent_chain[0], 1) == 3:
                    triples.append(i + 1)

        try:
            if doubles:
                mult = patterns.MULT_PREFIXES[len(doubles)]
                locants = ",".join(map(str, doubles))

                if mult == "" or mult[0] in {"a", "e", "i", "o", "u"}:
                    name += "-" + locants + "-" + mult + "en"
                else:
                    name += "a-" + locants + "-" + mult + "en"

            if triples:
                mult = patterns.MULT_PREFIXES[len(triples)]
                locants = ",".join(map(str, triples))

                if mult == "" or mult[0] in {"a", "e", "i", "o", "u"}:
                    name += "-" + locants + "-" + mult + "yn"
                else:
                    name += "a-" + locants + "-" + mult + "yn"

            if not doubles and not triples:
                name += "an"

        except KeyError:
            raise ValueError(f"Unsupported mult for unsaturation. Doubles: {len(doubles)} Triples: {len(triples)}")

        if not princip_groups:
            return name + "e"

        mult = patterns.MULT_PREFIXES[len(princips)]
        princip_locants = ",".join(map(str, princips))
        return name + "-" + princip_locants + "-" + mult + princip_groups[0].pattern.suffix

    # Shared Low-Level Helper Functions

    def _group_attachment(self, group: MotifMatch, chain_elem: str):
        """The atom(s) that count as "this candidate chain of element  contains
        group X".
        """
        center_idx = group.center_idx
        center_elem = self._molecule[center_idx].element.symbol

        if center_elem == chain_elem:
            return {center_idx}

        matched = set(group.matched_atoms)
        external = {n for n in self._molecule[center_idx].bonds if n not in matched}

        if len(external) == 1:
            return {center_idx} | external

        return {center_idx}

    def _ring_element(self, ring_atoms: list[int]) -> str:
        """Returns the most senior element PRESENT in the ring (checked in self._SENIOR_ELEMENTS
        order)
        """
        symbols = {self._molecule[idx].element.symbol for idx in ring_atoms}
        for element in self._SENIOR_ELEMENTS:
            if element in symbols:
                return element
        return next(iter(symbols))

    def _local_groups(self, atoms: set[int]) -> tuple[MotifMatch, ...]:
        """Returns every detected functional-group match whose atoms overlap "atoms"."""
        return tuple(group for group in self._groups if set(group.matched_atoms) & atoms)

    def _chain_candidates(self, idx: int, element: str, visited: set[int],
                          allow_pivot: bool = True, allow_change: bool = False) -> Union[list[list[int]], Ring]:
        """TODO:"""

        for ring in self._molecule.atom_rings:
            if idx in ring.atoms:
                return ring

        curr = self._molecule[idx]
        visited.add(idx)

        chains = self._follow_chain(idx, element, visited, allow_change)
        if not allow_pivot or curr.element.symbol != element or not chains or len(chains) < 2:
            return chains

        by_subtree = {}
        for chain in chains:
            first_step = chain[1]
            best = by_subtree.get(first_step)
            if best is None or self._score_chain(_Candidate(chain, element)) > self._score_chain(
                    _Candidate(best, element)):
                by_subtree[first_step] = chain

        final_chains = list(by_subtree.values())

        # idx as an interior atom, joining the two strongest subtrees.
        if len(final_chains) >= 2:
            top_two = sorted(final_chains, key=lambda c: self._score_chain(_Candidate(c, element)), reverse=True)[:2]
            merged = list(reversed(top_two[0][1:])) + [idx] + top_two[1][1:]
            final_chains.append(merged)

        return final_chains

    def _ring_system_atoms(self, ring: Ring) -> set[int]:
        """Every atom in ring's whole fused/bridged/spiro ring system, not just this one
        component ring -- found by walking .fused transitively."""
        seen = {ring}
        stack = [ring]
        while stack:
            curr = stack.pop()
            for other in curr.fused:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        return {atom for r in seen for atom in r.atoms}

    # Parent Chain Searchers

    def _find_parent_chain_no_p(self) -> _Candidate:
        """Returns a list of indices of the atoms in the parent backbone connected to the principal group's
        center atom.
        """
        visited = set()
        for ring in self._molecule.atom_rings:
            visited.update(set(ring.atoms))

        candidates = []
        seen = set()

        for atom in self._molecule:
            if atom.element.symbol not in self._SENIOR_ELEMENTS or atom.idx in visited:
                continue

            result = self._chain_candidates(atom.idx, atom.element.symbol, visited.copy())

            if isinstance(result, list):
                for chain in result:
                    if len(chain) > 1 or atom.element.symbol == "C":
                        key = tuple(chain)

                        if key not in seen:
                            seen.add(key)
                            candidates.append(_Candidate(chain, atom.element.symbol))

        for ring in self._molecule.atom_rings:
            elem = self._ring_element(ring.atoms)
            key = tuple(sorted(ring.atoms))
            if key not in seen:
                seen.add(key)
                candidates.append(_Candidate(ring.atoms, elem, ring))

        best = self._select_best_parent(candidates)
        return best

    def _find_parent_chain(self, princip_groups: list[MotifMatch]) -> _Candidate:
        """Returns a list of indices of the atoms in the parent backbone connected to the principal group's
        center atom.

        Parameters:
            - principal_idx: list of MotifMatches
                Motif Matches of the principal group that will have its parent backbone found.
        """
        visited = set()
        candidates = []
        seen = set()

        for group in princip_groups:
            if not group.pattern.inline:
                visited.update(set(group.matched_atoms))

        for ring in self._molecule.atom_rings:
            visited.update(set(ring.atoms))

        for group in princip_groups:
            principal_idx = group.center_idx

            for element in self._SENIOR_ELEMENTS:
                result = self._chain_candidates(principal_idx, element, visited.copy())

                if isinstance(result, list):
                    for chain in result:
                        if len(chain) > 1 or element == "C":
                            key = tuple(chain)

                            if key not in seen:
                                seen.add(key)
                                candidates.append(_Candidate(chain, element))

                elif isinstance(result, Ring) and self._ring_element(result.atoms) == element:
                    key = ("ring", tuple(sorted(result.atoms)))
                    if key not in seen:
                        seen.add(key)
                        candidates.append(_Candidate(result.atoms, element, result))

        princip_idxs = {group.center_idx for group in princip_groups}
        return self._select_best_parent(candidates, princip_idxs)

    # Subsitient Searchers

    def _substituent_hunting(self, parent_atoms: set[int], exclude: set[int]) -> list[tuple[int, int]]:
        """Returns substituent roots"""
        roots = []

        functional_atoms = set()
        for group in self._groups:
            functional_atoms.update(group.matched_atoms)

        for parent_idx in parent_atoms:
            for neighbour_idx in self._molecule[parent_idx].bonds:
                if (neighbour_idx not in parent_atoms and neighbour_idx not in exclude
                        and neighbour_idx not in functional_atoms):
                    roots.append((parent_idx, neighbour_idx))
        return roots

    def _name_substituent(self, root_idx: int, excluded: set[int], locant: int) -> _Substituent:
        """TODO"""
        root_elem = self._molecule[root_idx].element.symbol
        result = self._chain_candidates(root_idx, root_elem, excluded.copy(), allow_pivot=True, allow_change=False)

        if isinstance(result, Ring):
            own_chain = result.atoms
            cyclic = result
        else:
            own_chain = (
                max(result, key=lambda c: self._score_chain(_Candidate(c, root_elem)))
                if result else [root_idx]
            )
            cyclic = None

        subgraph = set(own_chain)
        local_groups = self._local_groups(subgraph)

        nested = tuple(
            self._name_substituent(child_root, excluded | subgraph, own_chain.index(parent_idx) + 1)
            for parent_idx, child_root in self._substituent_hunting(subgraph, excluded)
            if parent_idx not in excluded
        )

        return _Substituent(root_idx, tuple(own_chain), cyclic, local_groups, nested, locant)

    def _find_substituents(self) -> tuple[_Substituent, ...]:
        """Names every substituent branching directly off the parent chain."""
        return tuple(
            self._name_substituent(child_root, set(self._parent_chain), self._parent_chain.index(parent_idx) + 1)
            for parent_idx, child_root in self._substituent_hunting(self._parent_chain, set())
            if child_root not in self._parent_chain
        )

    # Locant Numbering

    def _number_parent(self, candidate: _Candidate) -> tuple[int]:
        if candidate.cyclic:
            return self._number_ring(candidate.chain)
        else:
            return self._number_chain(candidate.chain)

    def _number_chain(self, chain: list[int]):
        """TODO"""
        forward = chain
        reverse = chain[::-1]

        forward_profile = self._locant_profile(forward)
        reverse_profile = self._locant_profile(reverse)

        if forward_profile < reverse_profile:
            return forward
        else:
            return reverse

    def _number_ring(self, ring: list[int]):
        best_chain = None
        best_profile = None

        for start in range(len(ring)):
            # clockwise being purely a naming convention, not really matching the technical direction we are travelling
            clockwise = ring[start:] + ring[:start]
            profile = self._locant_profile(clockwise)
            if best_profile is None or profile < best_profile:
                best_profile = profile
                best_chain = clockwise

            anti_clockwise = clockwise[0:1] + clockwise[:0:-1]
            profile = self._locant_profile(anti_clockwise)
            if best_profile is None or profile <= best_profile:
                best_profile = profile
                best_chain = anti_clockwise

        return best_chain

    def _locant_profile(self, chain: list[int]) -> tuple[int, int, int, int, int, int]:
        """TODO"""
        hetero_loc = []
        hydrogen_loc = []
        princip_loc = []
        unsat_loc = []
        sub_loc = []
        sub_order_loc = []
        princip_idxs = [group.center_idx for group in self._groups if group.pattern == self._groups[0].pattern]
        subs_idxs = [parent_idx for parent_idx, _ in self._substituent_hunting(chain, set())]

        for i, idx in enumerate(chain):
            atom = self._molecule[idx]

            # a) Heteroatoms
            if atom.element.symbol in {"N", "P", "Si", "B", "O", "S"}:
                hetero_loc.append(i + 1)

            # b) Skip the indicated hydrogen for now TODO

            # c) Principle Groups
            if idx in princip_idxs:
                princip_loc.append(i + 1)

            # d) Multiple Bonds
            if i < len(chain) - 1:
                if atom.bonds.get(chain[i + 1], 1) >= 2:
                    unsat_loc.append(i + 1)

            # e) Substituents
            if idx in subs_idxs:
                sub_loc.append(i + 1)

            # f) Skip the Substituents order for now TODO
        return hetero_loc, hydrogen_loc, princip_loc, unsat_loc, sub_loc, sub_order_loc

    # Chain Walker

    def _follow_chain(self, curr_idx: int, element: str, visited: set = None, allow_change: bool = False) \
            -> list[list[int]]:
        """Follow a chain of "element" until it ends and returns said chain. Uses a DFS approach.

        Parameters:
            - curr_idx: integer
                Index of the current atom in the chain.
            - element: string
                Symbol of the element that makes up a chain.
        """
        if visited is None:
            visited = set()

        curr = self._molecule[curr_idx]
        visited.add(curr_idx)
        branches = []

        for neighbour_idx, _ in curr.bonds.items():
            if neighbour_idx in visited:
                continue

            neighbour = self._molecule[neighbour_idx]
            if neighbour.element.symbol == element or allow_change:
                visited.add(neighbour_idx)
                branches.extend(self._follow_chain(neighbour_idx, element, visited.copy(), allow_change=allow_change))

        if branches:
            return [[curr_idx] + chain for chain in branches] \
                if self._molecule[curr_idx].element.symbol == element or allow_change\
                else [chain for chain in branches]

        elif self._molecule[curr_idx].element.symbol == element or allow_change:
            return [[curr_idx]]

    # Seniority and Scoring

    def _select_best_parent(self, candidates: list[_Candidate], princip_idxs: Optional[set[int]] = None) \
            -> Optional[_Candidate]:
        """TODO:"""
        if not candidates:
            return None

        chains = []
        for candidate in candidates:
            chains.extend(candidate.chain)

        if princip_idxs:
            max_count = max(self._princip_count(candidate) for candidate in candidates)

            candidates = [
                candidate for candidate in candidates
                if self._princip_count(candidate) == max_count
            ]

        for element in self._SENIOR_ELEMENTS:
            element_candidates = [candidate for candidate in candidates if candidate.elem == element]
            if not element_candidates:
                continue
            rings = [candidate for candidate in element_candidates if candidate.cyclic]
            element_candidates = rings or element_candidates
            return self._score_chains(element_candidates)

        return None

    def _score_chains(self, chains: list[list[_Candidate]]) -> Optional[list[int]]:
        """"""
        if not chains:
            return None

        answer = max(chains, key=self._score_chain)

        if isinstance(answer.cyclic, Ring):
            if answer.cyclic.fused:
                chain = self._ring_system_atoms(answer.cyclic)
            else:
                chain = answer.cyclic.atoms
            answer = _Candidate(chain, answer.elem, answer.cyclic)

        return answer

    def _score_chain(self, candidate: _Candidate) -> tuple[int, ...]:
        """TODO"""
        chain = candidate.chain
        if not chain:
            return 0, 0, 0, 0

        length = len(chain)
        double_bonds = 0
        triple_bonds = 0
        heteroatom_sum = 0
        senior_hetero = 0
        heteroatom_prio = {"N": 6, "O": 5, "S": 4, "P": 3, "Si": 2, "B": 1}

        for i in range(len(chain) - 1):
            atom = self._molecule[chain[i]]
            bond_order = atom.bonds.get(chain[i + 1])

            if atom.element.symbol in heteroatom_prio:
                hetero_value = heteroatom_prio[atom.element.symbol]
                heteroatom_sum += hetero_value
                if hetero_value > senior_hetero:
                    senior_hetero = hetero_value

            if bond_order == 3:
                triple_bonds += 1
            elif bond_order == 2:
                double_bonds += 1

        multiple_bonds = double_bonds + triple_bonds

        if not candidate.cyclic:
            return length, multiple_bonds, double_bonds, len(self._substituent_hunting(candidate.chain, set()))

        else:
            return (senior_hetero, len(candidate.cyclic.fused), length, heteroatom_sum, multiple_bonds, double_bonds,
                    len(self._substituent_hunting(candidate.chain, set())))

    def _princip_count(self, candidate: _Candidate) -> int:
        chain_set = set(candidate.chain)
        return sum(1 for group in self._groups if chain_set & self._group_attachment(group, candidate.elem))


if __name__ == "__main__":
    mol = "CCCC"
    print(mol)
    generate_name(mol)
