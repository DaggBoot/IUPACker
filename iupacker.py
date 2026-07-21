"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from operator import concat

from motif_engine import MotifEngine, MotifMatch
from entities import Molecule, Atom
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
        elems = ["N", "P", "Si", "B", "O", "S", "C"]
        visited = set()

        for atom in self._molecule:
            if atom.element.symbol in {"N", "P", "Si", "B", "O", "S", "C"}:
                posssible_chain = self._follow_chain(atom.idx, atom.element.symbol, visited)

                if posssible_chain:
                    chains[atom.element.symbol].extend(posssible_chain)
                    visited.update(set(posssible_chain[0]))

                visited.add(atom.idx)

        for _, chain in chains.items():
            if chain:
                return max(chain, key=len)

        return None

    def _find_parent_chain(self, princip_groups: list[MotifMatch]) -> list[int]:
        """Returns a list of indices of the atoms in the parent backbone connected to the principal group's
        center atom.

        Parameters:
            - principal_idx: list of MotifMatches
                Motif Matches of the principal group that will have its parent backbone found.
        """
        chains = {"N": [], "P": [], "Si": [], "B": [], "O": [], "S": [], "C": []}
        visited = set()
        func_idxs= []

        for group in princip_groups:
            visited.update(set(group.matched_atoms))

        for group in princip_groups:
            print(group.pattern.name)
            principal_idx = group.center_idx
            func_idxs.append(principal_idx)
            principal = self._molecule[principal_idx]

            for element in {"N", "P", "Si", "B", "O", "S", "C"}:
                chain = self._follow_chain(principal_idx, element, visited.copy())
                if chain:
                    chains[element].extend(chain)
                print(chain)

        print(chains)

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
            visited.add(neighbour_idx)
            if neighbour.element.symbol == element:
                branches.extend(self._follow_chain(neighbour_idx, element, visited.copy()))

        if branches:
            return [[curr_idx] + chain for chain in branches]

        elif self._molecule[curr_idx].element.symbol == element:
            return [[curr_idx]]

    def _score_chains(self, chains: dict[str, list[list[int]]]) -> list[int]:
        """"""

    def _score_chain(self, chain: list[int]) -> int:
        """"""
        if not chain:
            return 0

        score = 0
        for idx in chain:
            atom = self._molecule[idx]



if __name__ == "__main__":
    mol = "C(=O)(O)C(C(C(=O)(O))CS(=O)(=O)O)CC(=O)O"
    print(mol)
    generate_name(mol)
