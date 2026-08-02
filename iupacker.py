"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from typing import Optional, NamedTuple, Union
from motif_engine import MotifEngine, MotifMatch
from entities import Molecule, Ring
from smiles_parser import SMILESParser
import patterns


def generate_name(smiles: str) -> str:
    """Returns the chemical name of the inputted chemcial formula as per IUPAC nomenclature.

    Parameters:
        - smiles: string
            The molecule's chemical formula, assumed to be in the SMILES formula under the OpenSMILES documentation.
            Will raise exceptions in the case of improper SMILES formatting or chemcially invalid formulae.
    """
    parser = SMILESParser()
    mol = parser.parse(smiles)
    namer = IUPACker(mol)
    name = namer.generate()
    for match in namer._groups:
        print(match.pattern.name)

    print(name)


def generate_name_from_mol(mol: Molecule):
    """Returns the chemical name of the inputted Molecule as per IUPAC nomenclature.

    Parameters:
        - mol: Molecule
            An assumed chemically valid Molecule Graph object.
    """
    namer = IUPACker(mol)
    return namer.generate()


class _Candidate(NamedTuple):
    """TODO:"""
    chain: list[int]
    elem: str


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
        - _groups: List of all MotifMatch objects that correspond to all the functional groups present.
                   Arranged descending based on priority.
        - _parent_chain: The atom indices of the chosen parent chain
    """
    _molecule: Molecule
    _groups: list[MotifMatch]
    _parent_chain: list[int]

    _SENIOR_ELEMENTS = ("N", "P", "Si", "B", "O", "S", "C")

    def __init__(self, molecule: Molecule = None):
        self._molecule = molecule
        self._groups = []
        self._parent_chain = []

    def generate(self) -> str:
        """Returns the chemical name of self._molecule as per IUPAC nomenclature of the Molecule.

        TODO: Rest of the stuff
        """
        engine = MotifEngine(self._molecule)
        self._groups = engine.match_all(patterns.ALL_PATTERNS)
        princip_groups = [group for group in self._groups if group.pattern == self._groups[0].pattern]
        if princip_groups:
            princip_chain = self._find_parent_chain(princip_groups)
        else:
            princip_chain = self._find_parent_chain_no_p()

        return str(len(princip_chain)) + str(princip_chain)

    def _find_parent_chain_no_p(self) -> list[int]:
        """Returns a list of indices of the atoms in the parent backbone connected to the principal group's
        center atom.
        """
        chains = {"N": [], "P": [], "Si": [], "B": [], "O": [], "S": [], "C": []}
        visited = set()

        for atom in self._molecule:
            if atom.element.symbol in chains:
                posssible_chain = self._follow_chain(atom.idx, atom.element.symbol, set())

                if posssible_chain:
                    chains[atom.element.symbol].extend(posssible_chain)
                    visited.update(set(posssible_chain[0]))

                visited.add(atom.idx)

        for elem, chain in chains.items():
            if chain:
                return _Candidate(max(chain, key=len), elem)

        return None

    def _find_parent_chain(self, princip_groups: list[MotifMatch]) -> list[int]:
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
            visited.update(set(group.matched_atoms))

        for ring in self._molecule.atom_rings:
            visited.update(set(ring.atoms))

        for group in princip_groups:
            principal_idx = group.center_idx

            for element in self._SENIOR_ELEMENTS:
                print(element)
                chains = self._chains_through(principal_idx, element, visited.copy())

                if isinstance(chains, list):
                    for chain in chains:
                        if len(chain) > 1 or element == "C":
                            key = tuple(chain)

                            if key not in seen:
                                seen.add(key)
                                candidates.append(_Candidate(chain, element))

                elif isinstance(chains, Ring):
                    seen.add(Ring)
                    candidates.append(Ring)

        princip_idxs = {group.center_idx for group in princip_groups}
        return self._select_best_parent(candidates, princip_idxs).chain

    def _select_best_parent(self, candidates: list[_Candidate], princip_idxs: Optional[set[int]] = None) \
            -> Optional[_Candidate]:
        """TODO:"""
        if not candidates:
            return None

        if princip_idxs:
            max_count = max(self._princip_count(candidate.chain, princip_idxs) for candidate in candidates)
            print(max_count)
            candidates = [
                candidate for candidate in candidates
                if self._princip_count(candidate.chain, princip_idxs) == max_count
            ]

        for element in self._SENIOR_ELEMENTS:
            element_candidates = [candidate for candidate in candidates if candidate.elem == element]
            if element_candidates:
                return self._score_chains(element_candidates)

        return None

    def _chains_through(self, idx: int, element: str, visited=None) -> Union[list[list[int]], Ring]:
        """TODO:"""
        curr = self._molecule[idx]
        visited = visited
        visited.add(idx)

        for ring in self._molecule.atom_rings:
            if idx in ring.atoms:
                return ring

        chains = self._follow_chain(idx, element, visited)
        print(chains)
        if curr.element.symbol != element or len(chains) < 2:
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

    def _follow_chain(self, curr_idx: int, element: str, visited=None) -> list[list[int]]:
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
            if neighbour.element.symbol == element:
                visited.add(neighbour_idx)
                branches.extend(self._follow_chain(neighbour_idx, element, visited.copy()))

        if branches:
            return [[curr_idx] + chain for chain in branches] if self._molecule[curr_idx].element.symbol == element \
                else [chain for chain in branches]

        elif self._molecule[curr_idx].element.symbol == element:
            return [[curr_idx]]

    def _score_chains(self, chains: list[list[int]]) -> Optional[list[int]]:
        """"""
        if not chains:
            return None

        answer = max(chains, key=self._score_chain)
        return answer

    def _score_chain(self, candidate: _Candidate) -> tuple[int, int, int]:
        """"""
        chain = candidate.chain
        if not chain:
            return 0, 0, 0

        length = len(chain)
        double_bonds = 0
        triple_bonds = 0

        for i in range(len(chain) - 1):
            atom = self._molecule[chain[i]]
            bond_order = atom.bonds.get(chain[i + 1])

            if bond_order == 3:
                triple_bonds += 1
            elif bond_order == 2:
                double_bonds += 1

        multiple_bonds = double_bonds + triple_bonds
        return length, multiple_bonds, double_bonds

    def _princip_count(self, chain: list[int], idxs) -> int:
        return sum(1 for idx in chain if (idx in idxs or any(b_idx in idxs for b_idx in self._molecule[idx].bonds)))


if __name__ == "__main__":
    mol = "CCCCCCCC(C(S(=O)(=O)(O))CS(=O)(=O)(O))(CCCS(=O)(=O)(O))"
    print(mol)
    generate_name(mol)
