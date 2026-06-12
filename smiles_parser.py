"""
Contains the main parser to convert a SMILES formula to a Molecule object.
"""
from typing import Optional, Tuple
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
    current: int
    tokens: list[str]

    order: int
    branch: list[int]
    rings: dict[int, int]   # digit to atom index

    def __init__(self):
        self.molecule = None
        self.pos = 0
        self.current = -1
        self.tokens = []

        self.order = 1
        self.branch = []
        self.rings = {}

    def _tokenize(self, smiles: str) -> None:
        """Convert a smiles string into a tokenized list, which is stored in self.tokens attribute

        Parameters:
            - smiles: string
                A chemical formula in SMILES form, assumed to be chemcially accurate.
        """
        self.tokens = re.findall(
            r'\[[^\]]+\]|'   # Bracket atoms: [CH3], [NH4+], [cH]
            r'%[0-9]{2}|'           # Two-digit rings: %12
            r'[0-9]|'               # Ring closures: 1, 2, 3
            r'[A-Z][a-z]?|'         # Aliphatic atoms: C, N, O, Cl, Br
            r'[a-z]|'               # Aromatic atoms: c, n, o, s, p
            r'[-=#$:\.()]',         # Bonds and branches: -, =, #, $, :, ., (, )
            smiles
        )

    def _reset(self) -> None:
        """Resets all instance atteibutes in preperation for the parsing of a new SMILES formula."""
        self.molecule = Molecule()
        self.pos = 0
        self.current = -1
        self.tokens = []

        self.order = 1
        self.branch = []
        self.rings = {}

    def parse(self, smiles: str) -> Molecule:
        """A naive parser, Returns a molecule based on the """
        self._reset()
        self._tokenize(smiles)

        while self.pos < len(self.tokens):
            token = self.tokens[self.pos]

            if token == "(":
                self._handle_branch_start()
            elif token == ")":
                self._handle_branch_end()
            elif token.isdigit():
                self._handle_rings(int(token))
            elif token.startswith("%"):
                self._handle_rings(int(token[1:]))
            elif token in {"-", "=", "#", "$", ":", "."}:
                self.order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}[token]
                self.pos += 1
            else:
                self._handle_atom(token)

        return self.molecule

    def _handle_branch_start(self) -> None:
        """Opens a branch and stores the parent atom"""
        if self.current == -1:
            raise SyntaxError("Branch begun before any parent atom")
        self.branch.append(self.current)
        self.pos += 1

    def _handle_branch_end(self) -> None:
        """Closes a branch and connects it to the parent atom"""
        if not self.branch:
            raise SyntaxError("Unmatched closing parenthesis")
        self.current = self.branch.pop()
        self.pos += 1

    def _handle_rings(self, num: int) -> None:
        """Handles ring opening and closure as per the SMILES convention

        Parameters:
            - num: integer
                The index of the atom associated with this number in the SMILES formula.
        """
        if self.current == -1:
            raise SyntaxError(f"Ring closure for {num} before any parent atom")

        if num in self.rings:
            # For ring closure
            ring_origin = self.rings[num]
            if ring_origin == self.current:
                raise SyntaxError(f"Ring closure for {num} connects atom to itself")

            self.molecule.add_bonds(num, ring_origin, self.order)
            self.bond = 1
            del self.rings[num]
        else:
            # For ring opening
            self.rings[num] = self.current

        self.pos += 1

    def _handle_atom(self, token: str) -> None:
        """Parse and add an atom to the molecule

        Parameters:
            - token: string
                Symbol that represents the atom to be added to the molecule.
        """
        # Get the properties of the token
        symbol, is_aromatic, charge = self._parse_atomic_token(token)

        idx = self.molecule.add_atom(symbol, is_aromatic)
        if charge != 0:
            self.molecule.get_atom().charge = charge

        if self.current > -1:
            self.molecule.add_bonds(self.current, idx, self.order)
            self.order = 1

        self.current = idx
        self.pos += 1

    def _parse_atomic_token(self, token: str) -> Tuple[str, bool, int]:
        """Returns the chemical symbol, the aromaticity and charge of an atom token

        Parameters:
            - token: string
                Symbol that represents the atom to be added to the molecule.
        """
        is_aromatic = False
        charge = 0
        symbol = None

        if token.startswith("["):
            inner = token[1:-1]
            if re.match(r'^(\d+)', inner):
                raise NotImplementedError

        else:
            if token.islower() and token in {"c", "n", "s", "o", "p"}:
                is_aromatic = True
                symbol = token.upper()
            else:
                symbol = token

        return symbol, is_aromatic, charge


if __name__ == "__main__":
    _Element.load_data("Periodic.json")
    parser = SMILESParser()

    test_cases = [
        "C",
        "CC",
        "C=C",
        "CCO",
        "CC(C)C",
        "c1ccccc1",
        "C1CCCC1",
        "C=CC=C"
    ]

    for smiles in test_cases:
        print(smiles)
        print(parser.parse(smiles))
