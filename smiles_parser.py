"""
Contains the main parser to convert a SMILES formula to a Molecule object.
"""
from typing import Optional
from entities import _Element, Atom, Molecule
import re


class SMILESParser:
    """Parser Object, returns a molecule from an input of a SMILES formula.

    Instance Attributes:
        - molecule: The Molecule object formed so far from the SMILES formula.
        - pos: The current position in SMILES formula.
        - smiles: The original SMILES formula.
    """
    molecule: Optional[Molecule]
    pos: int
    tokens: list[str]

    _pending_bond_order: int

    def __init__(self):
        self.molecule = None
        self.pos = 0
        self.smiles = []

        self._pending_bond_order = 0

    def _tokenize(self, smiles: str) -> None:
        """Convert a smiles string into a tokenized list"""
        self.tokens = re.findall(
            r'\[[^\]]+\]|'
            r'%[0-9]{2}|'
            r'[A-Z][a-z]?|'
            r'[a-z]|'
            r'[0-9]|'
            r'[=#()]', smiles
        )

    def parse(self, smiles: str):
        self._tokenize(smiles)

        while self.pos < len(self.smiles):
            token = self.tokens[self.pos]

            if token in {"=", "#"}:
                self._pending_bond_order = {"=": 2, "#": 3}[token]
