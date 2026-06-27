"""
Contains the main parser to convert a SMILES formula to a Molecule object.
"""
from typing import Optional
from entities import _Element, Molecule
import re


class SMILESParser:
    """Parser Object, returns a molecule from an input of a SMILES formula.

    Instance Attributes:
        - molecule: The Molecule object formed so far from the SMILES formula.
        - pos: The current position in SMILES formula.
        - smiles: The original SMILES formula.

    TODO: MUST UPDATE
    """
    molecule: Optional[Molecule]
    pos: int
    current: int
    tokens: list[str]

    order: int
    branch: list[int]
    rings: dict[int, int]   # digit to atom index

    def __init__(self):
        _Element.load_data("periodic.json")
        self.molecule = None
        self.pos = 0
        self.current = -1
        self.tokens = []

        self.order = 1
        self.branch = []
        self.rings = {}

    def parse(self, smiles: str) -> Molecule:
        """Returns the Molecule graph object described by the SMILES formula inputted. Validates if the resultant object
        is a chemically-valid molecule, and raises an exception if not.

        Parameters:
            - smiles: string
            A chemical formula in SMILES form.
        """
        self._syntax_validate(smiles)
        self._naive_parse(smiles)
        self.molecule.valence_validate()
        return self.molecule

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

    def _naive_parse(self, smiles: str) -> Molecule:
        """Returns the Molecule graph object described by the SMILES formula inputted, assuming the formula describes a
        chemically valid molecule.

        Parameters:
            - smiles: string
            A chemical formula in SMILES form, assumed to be chemcially accurate.

        Exceptions:
            - Will raise an error in the case of unbalanced parenthesis, self-bonded rings,
            and isotopes (due to not being implmented yet).
        """
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
                The integer that uniquely identifies which ring an end-atom is connected to.
        """
        if self.current == -1:
            raise SyntaxError(f"Ring closure for {num} before any parent atom")

        if num in self.rings:
            # For ring closure
            ring_origin = self.rings[num]
            if ring_origin == self.current:
                raise SyntaxError(f"Ring closure for {num} connects atom to itself")

            self.molecule.add_bonds(self.current, ring_origin, self.order)
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
        symbol, is_aromatic, charge, isotope, exp_h = self._parse_atomic_token(token)

        idx = self.molecule.add_atom(symbol, is_aromatic, charge, isotope, exp_h)

        if self.current > -1:
            self.molecule.add_bonds(self.current, idx, self.order)
            self.order = 1

        self.current = idx
        self.pos += 1

    @staticmethod
    def _parse_atomic_token(token: str) -> tuple[str, bool, int, Optional[int], Optional[int]]:
        """Returns the chemical symbol and properties of a token.

        Parameters:
            - token: string
                Symbol that represents the atom to be added to the molecule.
        """
        is_aromatic = False
        charge = 0
        isotope = None
        exp_h = None

        if token.startswith("["):
            inner = re.match(r'^(\d+)?([A-z][a-z]?|[a-z])(H\d?)?([+-]+\d*)?', token[1:-1])
            if not inner:
                raise ValueError(f"Invaid bracket atom: {token}")

            isotope, token, h_str, charge_str = inner.groups()

            if charge_str:
                if charge_str.startswith("+"):
                    charge = int(charge_str[1:]) if len(charge_str) > 1 else 1
                if charge_str.startswith("-"):
                    charge = -int(charge_str[1:]) if len(charge_str) > 1 else -1
            else:
                charge = 0

            if h_str:
                exp_h = int(h_str[1:]) if len(h_str) > 1 else 1
            else:
                exp_h = None

        if token.islower() and token in {"c", "n", "s", "o", "p"}:
            is_aromatic = True
            symbol = token.upper()
        else:
            symbol = token

        return symbol, is_aromatic, charge, isotope, exp_h

    @staticmethod
    def _syntax_validate(smiles: str) -> None:
        """Validates the syntax of the SMILES formula, ensuring it has mathcing digits and brackets. Raises an exception
        for any unmatched digits or brackets.

        Parameters:
            - smiles: string
            A chemical formula in SMILES form.
        """
        para_balance = 0
        sqr_balance = 0
        ring_digits = set()
        i = 0

        while i < len(smiles):
            ch = smiles[i]

            # parentheses
            if ch == '(':
                para_balance += 1
            elif ch == ')':
                para_balance -= 1
                if para_balance < 0:
                    raise SyntaxError("Unmatched closing parenthesis")

            # brackets
            elif ch == '[':
                sqr_balance += 1
            elif ch == ']':
                sqr_balance -= 1
                if sqr_balance < 0:
                    raise SyntaxError("Unmatched closing bracket")

            # ring closures (single digit)
            elif ch.isdigit() and sqr_balance < 1:
                if ch in ring_digits:
                    ring_digits.remove(ch)  # Closing
                else:
                    ring_digits.add(ch)  # Opening

            # two-digit ring closures
            elif ch == '%' and sqr_balance < 1:
                if i + 2 >= len(smiles):
                    raise SyntaxError("Incomplete % ring closure")
                digit = smiles[i + 1:i + 3]
                if digit in ring_digits:
                    ring_digits.remove(digit)
                else:
                    ring_digits.add(digit)
                i += 2

            i += 1

        if para_balance != 0:
            raise SyntaxError("Unclosed parentheses")
        if sqr_balance != 0:
            raise SyntaxError("Unclosed brackets")
        if ring_digits:
            raise SyntaxError(f"Unmatched ring closures: {ring_digits}")


if __name__ == "__main__":
    parser = SMILESParser()

    # INVALID_SMILES = [
    #     "C12C",
    #     "C(C",
    #     "C)CC",
    #     "C1CC",
    #     "C(C)(C)(C)(C)(C)"
    # ]
    #
    # test_cases = [
    #     "CCC(=O)O",
    #     "C(C)(C)N",
    #     "C1CC(CC)C1",
    #     "c%11cccc%11",
    #     "C#CC=O",
    # ]

    print(parser.parse("[nH+]1ccccc1"))

    # for smiles in test_cases:
    #     print(smiles)
    #     print(parser.parse(smiles))
    #
    # for smiles in INVALID_SMILES:
    #     try:
    #         mol = parser.parse(smiles)
    #         print(f"✓ {smiles}")
    #     except Exception as e:
    #         print(f"✗ {smiles}: {e}")
