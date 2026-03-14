"""
File to handle all objects required to represent a hydrocarbon as a Graph.
"""
from __future__ import annotations
from typing import Optional, Any


class _Element:
    """Immutable representation of a chemcial element

    Class Attributes:
        - _elements: A dictionary containing the data of each element
        - _instances: A dictionary containing all exisiting instances of elements.

    Instance Attributes:
        - symbol: Name of this element.
        - valence: The number of valence electrons of the _Element instance.
    """
    _elements: dict[str, dict] = {}
    _instances: dict[str, _Element] = {}

    symbol: str
    valence: list[int]

    @classmethod
    def load_data(cls, file_path: str) -> None:
        """Loads the periodic element data needed to set up the private element classes"""
        import json
        with open(file_path) as f:
            cls._elements = json.load(f)

    @classmethod
    def get(cls, symbol: str):
        """Gets the right instance of a given element"""
        if symbol not in cls._elements:
            raise ValueError(f"Unknown element: {symbol}")

        if symbol not in cls._instances:
            cls._instances[symbol] = cls(symbol)
        return cls._instances[symbol]

    def __init__(self, symbol):
        self.symbol = symbol
        self.valence = _Element._elements[symbol]["valence"]


class Atom:
    """An atom in a molecule

    Instance Attributes:
        - element: The element is this atom comprised of.
        - hydrogens: The number of hydrogens is this connected to.
        - bonds: A list of tuples containing the other connected atoms and their bond order.
    """
    element: _Element
    hydrogens: int
    bonds: list[tuple[Atom, int]]

    def __init__(self, element: _Element, hydrogens: int = 0, bonds: list[tuple[Atom, int]] = None):
        self.element = element
        self.hydrogens = hydrogens
        self.bonds = bonds or []

    def add_bond(self, other_atom: Atom, order: int = 1) -> None:
        """Adds a bond to another atom"""
        for bonded_atom, _ in self.bonds:
            if bonded_atom is other_atom:
                return
        self.bonds.append((other_atom, order))


class Molecule:
    """An organic molecule

    Instance Attributes:
        -_atoms: A collection of atoms that form a graph which represents the molecule.
    """
    formula: str
    pos: int
    tokens: list[str]

    _atoms: list[Atom]
    _pending_bond_order: int

    def __init__(self, formula: str) -> None:
        self.formula = formula
        self._pending_bond_order = 1
        self.pos = 0
        self.tokens = self._extract(formula)
        self._atoms = self._parse_chain(None)

    def _parse_chain(self, parent_atom: Optional[Atom] = None) -> list[Atom]:
        """Parse atoms until a ")" or the end. Returns atoms with their bomds in this chain."""
        chain_atoms = []

        while self.pos < len(self.tokens):
            token = self.tokens[self.pos]

            if token == 'C':
                atom = self._parse_carbon()
                self._connect_to_parent(atom, parent_atom)
                chain_atoms.append(atom)
                parent_atom = atom

            elif token == 'O':
                atom, becomes_parent = self._parse_oxygen()
                self._connect_to_parent(atom, parent_atom)
                chain_atoms.append(atom)

                if becomes_parent:
                    parent_atom = atom
                else:
                    parent_atom = None

            elif token == 'N':
                atom, becomes_parent = self._parse_nitrogen()
                self._connect_to_parent(atom, parent_atom)
                chain_atoms.append(atom)

                if becomes_parent:
                    parent_atom = atom
                else:
                    parent_atom = None

            elif token == "S":
                atom, becomes_parent = self._parse_sulfur()
                self._connect_to_parent(atom, parent_atom)
                chain_atoms.append(atom)

                if becomes_parent:
                    parent_atom = atom
                else:
                    parent_atom = None

            elif token in {"Br", "F", "Cl", "I"}:
                atom, becomes_parent = self._parse_halogen(token)
                self._connect_to_parent(atom, parent_atom)
                chain_atoms.append(atom)

                if becomes_parent:
                    parent_atom = atom
                else:
                    parent_atom = None

            elif token == '(':
                self.pos += 1
                branch_atoms = self._parse_chain(parent_atom)
                chain_atoms.extend(branch_atoms)
                self.pos += 1

            elif token == ')':
                break

            elif token in ['-', '=', '#']:
                self._pending_bond_order = {'-': 1, '=': 2, '#': 3}[token]
                self.pos += 1

            else:
                self.pos += 1

        return chain_atoms

    def _parse_carbon(self) -> Atom:
        carb = Atom(_Element.get("C"))
        self.pos += 1

        if self.pos < len(self.tokens) and self.tokens[self.pos] == "H":
            self.pos += 1
            if self.pos < len(self.tokens) and self.tokens[self.pos].isdigit():
                carb.hydrogens = int(self.tokens[self.pos])
            else:
                carb.hydrogens = 1

        return carb

    def _parse_oxygen(self) -> tuple[Atom, bool]:
        oxy = Atom(_Element.get("O"))
        if self._pending_bond_order != 1:
            return oxy, False

        self.pos += 1

        if self.tokens[self.pos] == "H":
            oxy.hydrogens = 1
            return oxy, False

        elif self.tokens[self.pos] in {"C", "O", "("}:
            return oxy, True

        return oxy, False

    def _parse_nitrogen(self) -> tuple[Atom, bool]:
        nit = Atom(_Element.get("N"))
        if self._pending_bond_order != 1:
            return nit, False

        self.pos += 1

        if self.pos < len(self.tokens) and self.tokens[self.pos] == "H":
            self.pos += 1
            if self.pos < len(self.tokens) and self.tokens[self.pos].isdigit():
                nit.hydrogens = int(self.tokens[self.pos])
            else:
                nit.hydrogens = 1

        elif self.tokens[self.pos + 1] in {"C", "O", "("}:
            return nit, True

        return nit, False

    # def _parse_sulfur(self) -> tuple[Atom, bool]:
    #     sul = Atom(_Element.get("S"))
    #
    #     return sul

    def _parse_halogen(self, token) -> tuple[Atom, bool]:
        hal = Atom(_Element.get(token))
        self.pos += 1

        return hal, False

    def _connect_to_parent(self, atom: Atom, parent_atom: Atom):
        """Connect atom to parent with pending bond order"""
        if parent_atom:
            self._add_bonds(parent_atom, atom, self._pending_bond_order)
            self._pending_bond_order = 1

    @staticmethod
    def _extract(string: str) -> list[str]:
        """Extracts the element characters from a string"""
        import re
        return [el for el in re.findall(r'[A-Z][a-z]?|\d+|[#=()-]', string)]

    @staticmethod
    def _add_bonds(atom1: Atom, atom2: Atom, order: int = 1) -> None:
        """Adds the bonds between atoms in the Molecule"""
        atom1.add_bond(atom2, order)
        atom2.add_bond(atom1, order)

    def __str__(self):
        result = f"{self.formula} \n"
        for atom in self._atoms:
            result += f"{atom.element.symbol} Connected to \n"
            for bonded in atom.bonds:
                result += f"   {bonded[0].element.symbol} with order {bonded[1]} \n"
            result += f"    Has {atom.hydrogens} hydrogens\n \n"
        return result


if __name__ == "__main__":
    pass
