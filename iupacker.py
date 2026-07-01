"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from operator import concat

from motif_engine import MotifEngine, MotifMatch
from entities import Molecule
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
        princip_group = self._groups[0]
        print(princip_group.pattern)
        princip_chain = self._find_parent_chain(princip_group.center_idx)

        return str(len(princip_chain))

    def _find_parent_chain(self, principal_idx: int) -> list[int]:
        """Returns a list of indices of the carbons in the longest carbon chain connected to the principal group's
        center atom.

        Parameters:
            - principal_idx: integer
                Index of the central atom of the principal group that will have its longest carbon chain found.
        """
        principal = self._molecule[principal_idx]

        if principal.element.symbol == "C":
            return self._follow_chain(principal_idx, {principal_idx})

        carbon_neighbours = [idx for idx, _ in principal.bonds.items() if self._molecule[idx].element.symbol == "C"]

        if not carbon_neighbours:
            return []

        if len(carbon_neighbours) == 1:
            return self._follow_chain(carbon_neighbours[0])

        branches = []
        for neighbour_idx in carbon_neighbours:
            branch = self._follow_chain(neighbour_idx)
            branches.append(branch)

        # Now that we have the longest branches from each carbon neighbour, we now must consider how to construct the
        # longest chain by connecting the branches.

        branches.sort(key=len, reverse=True)

        if len(branches) >= 2:
            return branches[0] + [principal_idx] + branches[1]
        else:
            return branches[0] + [principal_idx]

    def _follow_chain(self, curr_idx: int, visited=None) -> list[int]:
        """Follow a chain of carbons until it ends and returns said chain. Uses a DFS approach.

        Parameters:
            - prev_idx: integer
                Index of the previous atom in the chain.
            - curr_idx: integer
                Index of the current atom in the chain.
        """
        if visited is None:
            visited = set()

        curr = self._molecule[curr_idx]
        branches = []

        for neighbour_idx, _ in curr.bonds.items():
            if neighbour_idx in visited:
                continue

            neighbour = self._molecule[neighbour_idx]
            visited.add(neighbour_idx)
            if neighbour.element.symbol == "C":
                branches.append(self._follow_chain(neighbour_idx, visited))

        if branches:
            longest_chain = max(branches, key=len)
            # noinspection PyTypeChecker
            return [curr_idx] + longest_chain
        else:
            return [curr_idx]


if __name__ == "__main__":
    generate_name("CCCCCCCCCCCC(CS(=O)(=O)O)CC(=O)O")
