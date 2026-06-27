"""
Contains the namer of a chemical molecule under IUPAC nomenclature.
"""
from motif_engine import MotifEngine
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
    print(namer.generate())


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

    def generate(self) -> list[str]:
        """Returns the chemical name of self._molecule as per IUPAC nomenclature of the Molecule.
        """
        engine = MotifEngine(self._molecule)
        matches = engine.match_pattern(patterns.CARBOXYLIC_ACID)
        return [match.pattern.name for match in matches]


if __name__ == "__main__":
    generate_name("O=C(CCC)O")
