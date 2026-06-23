"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from typing import Optional
from entities import Atom, Molecule
from smiles_parser import SMILESParser


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
    return namer.generate()


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
        - _groups: A list of FunctionalGroup objects found in the molecule, ordered by priority.
        - _parent_chain: The atom indices of the chosen parent chain
        - _substituents: Functional groups not part of the parent chain
    """
    _molecule: Molecule
    # _groups: list[FunctionalGroup]
    _parent_chain: list[int]
    # _subsituents: list[FunctionalGroup]

    def __init__(self, molecule: Molecule = None):
        self._molecule = molecule
        self._groups = []
        self._parent_chain = []
        self._subsituents = []

    def generate(self) -> str:
        """Returns the chemical name of self._molecule as per IUPAC nomenclature of the Molecule.
        """
        return "placeholder"

    #  CURRENTLY UNDER ARCHITECURAL CHANGE. DUE TO BE CHANGED.

    # def _locate_func_groups(self) -> Optional[FunctionalGroup]:
    #     """Finds all functional groups (that were implemented, see Documentation for more information), adds them to
    #     self._groups. Returns highest priority group.
    #     """
    #     for atom in self._molecule:
    #         if atom.element.symbol == 'C':
    #             self.is_carboxylic_acid(atom)
    #         if atom.element.symbol == 'S':
    #             self.is_sulfonic_acid(atom)
    #
    #     if self._groups:
    #         return self._groups[0]
    #     else:
    #         return None
    #
    # def is_carboxylic_acid(self, atom: Atom) -> None:
    #     """Checks if atom is part of a Carboxylic acid functional group and adds said group to self._groups.
    #
    #     Parameters:
    #         - atom: Atom
    #             Assumed to be carbon atom, will not work as intended otherwise.
    #     """
    #     oh_idx = None
    #     carbonyl_idx = None
    #     alkyl_h_connector = atom.get_hydrogen_count() == 1
    #
    #     for neighbour_idx, order in atom.bonds.items():
    #         neighbour = self._molecule[neighbour_idx]
    #
    #         if neighbour.element.symbol == "O" and neighbour.bond_sum() <= 2:
    #
    #             # Checks for the C=O
    #             if order == 2:
    #                 carbonyl_idx = neighbour_idx
    #
    #             # Checks for the C-OH
    #             elif order == 1 and neighbour.get_hydrogen_count() == 1:
    #                 oh_idx = neighbour_idx
    #         else:
    #             if neighbour.element.symbol == "C":
    #                 alkyl_h_connector = True
    #
    #     if carbonyl_idx is not None and oh_idx is not None and alkyl_h_connector and atom.charge == 0:
    #         # Technically the atom.charge check is a very niche and unlikely check, but it's there for my sanity.
    #         self._groups.append(FunctionalGroup(
    #             "carboxylic acid", atom.idx, [atom.idx, carbonyl_idx, oh_idx]))
    #
    # def is_sulphonic_acid(self, atom: Atom):
    #     """Checks if atom is part of a sulphonic acid functional group and adds said group to self._groups.
    #     (Calls _is_sulfonic_acid)
    #     """
    #     self.is_sulfonic_acid(atom)
    #
    # def is_sulfonic_acid(self, atom: Atom):
    #     """Checks if atom is part of a sulfonic acid functional group and adds said group to self._groups.
    #
    #     Parameters:
    #         - atom: Atom
    #             Assumed to be carbon atom, will not work as intended otherwise.
    #     """
    #     oh_idx = None
    #     sulfonyl_idxs = []
    #     alkyl_connector = False
    #
    #     for neighbour_idx, order in atom.bonds.items():
    #         neighbour = self._molecule[neighbour_idx]
    #
    #         if neighbour.element.symbol == "O" and neighbour.bond_sum() <= 2:
    #
    #             if order == 2:
    #                 sulfonyl_idxs.append(neighbour_idx)
    #
    #             elif order == 1 and neighbour.get_hydrogen_count() == 1:
    #                 oh_idx = neighbour_idx
    #         else:
    #             if neighbour.element.symbol == "C":
    #                 alkyl_connector = True
    #
    #     if len(sulfonyl_idxs) == 2 and oh_idx is not None and alkyl_connector:
    #         self._groups.append(FunctionalGroup(
    #             "sulfonic acid", atom.idx, [atom.idx, sulfonyl_idxs[0], sulfonyl_idxs[1], oh_idx]))
    #
    # def is_anhydrides(self, atom: Atom):
    #     """Checks if atom is part of an anhydride functional group and adds said group to self._groups.
    #
    #     Parameters:
    #         - atom: Atom
    #             Assumed to be oxygen atom, will not work as intended otherwise.
    #     """
    #     carbonyl_idxs = []
    #     for neighbour_idx, order in atom.bonds.items():
    #         neighbour = self._molecule[neighbour_idx]
    #
    #         if neighbour.element.symbol == "C" and neighbour.bond_sum() == 4:
    #
    #             carbonyl_idx = self.has_carbonyl(neighbour)
    #             if carbonyl_idx:
    #                 carbonyl_idxs.append(carbonyl_idx)
    #
    # def has_carbonyl(self, atom: Atom) -> list[int]:
    #     """Checks if atom is part of a carbonyl group and returns a list of indices from the group.
    #
    #     Parameters:
    #         - atom: Atom
    #             Assumed to be carbon atom, will not work as intended otherwise.
    #     """
    #     carbonyl_idxs = []
    #     for neighbour_idx, order in atom.bonds.items():
    #         neighbour = self._molecule[neighbour_idx]
    #
    #         if neighbour.element.symbol == "O" and neighbour.bond_sum() <= 2 and order == 2:
    #             carbonyl_idxs.append(neighbour_idx)
    #
    #     return carbonyl_idxs
    #
    # def has_hydroxyl(self, atom: Atom) -> list[int]:
    #     """Checks if atom is part of a hydroxyl group and returns a list of indices from the group.
    #
    #     Parameters:
    #         - atom: Atom
    #             Assumed to be carbon atom, will not work as intended otherwise.
    #     """
    #     hydroxyl_idxs = []
    #     for neighbour_idx, order in atom.bonds.items():
    #         neighbour = self._molecule[neighbour_idx]
    #
    #         if (neighbour.element.symbol == "O" and neighbour.bond_sum() <= 2
    #                 and order == 1 and neighbour.get_hydrogen_count() == 1):
    #             hydroxyl_idxs.append(neighbour_idx)
    #
    #     return hydroxyl_idxs


if __name__ == "__main__":
    generate_name("CS(=O)(=O)O")
