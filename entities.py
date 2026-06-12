"""
Contains entities that represent a molecule.
"""
from __future__ import annotations
from typing import Optional, Any


class ValenceError(Exception):
    """Custom Exception for when the valence balance of a Molecule is invalid"""
    atom: str
    hydrogens: int
    valences: list[int]

    def __init__(self, atom: str, hydrogens: int):
        super().__init__()
        self.atom = atom
        self.hydrogens = hydrogens

    def __str__(self):
        return f"Invalid, Atom w/ index {self.atom.idx} has more bonds than its {self.atom.valences} available valences"


class _Element:
    """Immutable representation of a chemcial element

        Class Attributes:
            - _elements: A dictionary containing the data of each element
            - _instances: A dictionary containing all exisiting instances of elements.

        Instance Attributes:
            - symbol: Name of this element.
            - valence: The possible valences for this element.
        """
    _elements: dict[str, dict] = {}
    _instances: dict[str, _Element] = {}

    symbol: str
    valences: list[int]

    @classmethod
    def load_data(cls, file_path: str) -> None:
        """Loads the periodic element data needed to set up the private element classes"""
        import json
        with open(file_path) as f:
            cls._elements = json.load(f)

    @classmethod
    def get(cls, symbol: str) -> _Element:
        """Gets the right instance of a given element"""
        if symbol not in cls._elements:
            raise ValueError(f"Unknown element: {symbol}")

        if symbol not in cls._instances:
            cls._instances[symbol] = cls(symbol)
        return cls._instances[symbol]

    def __init__(self, symbol):
        self.symbol = symbol
        self.valences = _Element._elements[symbol]["valences"]


class Atom:
    """An atom in a molecule

    Instance Attributes:
        - element: The element is this atom comprised of.
        - idx: The index of this atom for fast referencing.
        - bonds: A list of tuples containing the other connected atoms' index and their bond order.
        - is_aromatic: A boolean that notes if the atom is a part of an aromatic ring.
        - charge: Stores the charge of this atom.
    """
    element: _Element
    idx: int
    bonds: dict[int, int]  # Atom index -> Order
    is_aromatic: bool
    charge: int
    isotope: Optional[int]

    def __init__(self, element: _Element, idx: int, is_aromatic: bool = False):
        self.element = element
        self.idx = idx
        self.bonds = {}
        self.is_aromatic = is_aromatic
        self.charge = 0
        self.isotope = None

    def add_bond(self, other_idx: int, order: int = 1) -> None:
        """Adds a bond to another atom"""
        if other_idx not in self.bonds:
            self.bonds[other_idx] = order

    def get_hydrogen_count(self) -> int:
        """Returns the number of Hydrogens attached to this atom.

        Precondition:
            - Relies on the bonds to be checmially valid.
        """
        bond_sum = sum(self.bonds.values())
        h_count = 0
        for valence in self.element.valences:
            h_count = valence - bond_sum - self.charge
            if h_count >= 0:
                return h_count

        raise ValenceError(self.element.symbol, -h_count)


class Molecule:
    """ The graph representation for an Organic Molecule

    Instance Attributes:
        -_atoms: A collection of Atom objects (vertices) that form a graph which represents the molecule.
    """

    _atoms: list[Atom]

    def __init__(self):
        self._atoms = []

    def __str__(self):
        result = ""
        for atom in self._atoms:
            result += f"{atom.element.symbol}{atom.idx} Connected to \n"
            for bond in atom.bonds:
                result += f"   {self._atoms[bond].element.symbol}{bond} with order {atom.bonds[bond]} \n"
        return result

    def add_atom(self, symbol: str, is_aromatic: bool = False) -> int:
        """Adds atom to the Molecule object and returns its index"""
        idx = len(self._atoms)
        atom = Atom(_Element.get(symbol), idx, is_aromatic)
        self._atoms.append(atom)
        return idx

    def add_bonds(self, idx1: int, idx2: int, order: int = 1):
        """Adds an edge between the two atoms referenced by index, with a weight based on bond order"""
        if not (0 <= idx1 < len(self._atoms) and 0 <= idx2 < len(self._atoms)):
            raise IndexError

        self._atoms[idx1].add_bond(idx2, order)
        self._atoms[idx2].add_bond(idx1, order)

    def get_atom(self, idx: int) -> Atom:
        """Returns atom object based on index."""
        return self._atoms[idx]
