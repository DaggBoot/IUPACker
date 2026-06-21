"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from typing import Tuple
from entities import Atom, Molecule, FunctionalGroup
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
    _groups: list[FunctionalGroup]
    _parent_chain: list[int]
    _subsituents: list[FunctionalGroup]

    def __init__(self, molecule: Molecule = None):
        self._molecule = molecule

    def generate(self) -> str:
        """Returns the chemical name of self._molecule as per IUPAC nomenclature of the Molecule.
        """

    def _locate_func_groups(self) -> FunctionalGroup:
        """Finds all functional groups (that were implemented, see Documentation for more information), adds them to
        self._groups. Returns highest priority group.
        """
        for atom in self._molecule:
            if atom.element.symbol == 'C':
                self._is_carboxylic_acid(atom)

        return self._groups[0]

    def _is_carboxylic_acid(self, atom: Atom) -> None:
        """Checks if atom is part of a Carboxylic acid functional group and adds said group to self._groups.

        Parameters:
            - atom: Atom
                Assumed to be carbon atom, will not work as intended otherwise.
        """
        oh_idx = None
        carbonyl_idx = None
        for neighbour_idx, order in atom.bonds.items():
            neighbour = self._molecule[neighbour_idx]

            if neighbour.element.symbol == "O" and neighbour.bond_sum() == 2:

                # Checks for the C=O
                if order == 2:
                    carbonyl_idx = neighbour_idx

                # Checks for the C-OH
                elif order == 1 and (neighbour.get_hydrogen_count() == 1 or neighbour.exp_h == 1):
                    oh_idx = neighbour_idx

        if carbonyl_idx and oh_idx:
            self._groups.append(FunctionalGroup(
                "carboxylic acid", atom.idx, [atom.idx, carbonyl_idx, oh_idx]))
