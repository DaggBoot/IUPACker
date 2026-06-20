"""
Contains entities that represent a molecule.
"""
from __future__ import annotations
from typing import Optional, Any


class ValenceError(Exception):
    """Custom Exception for when the valence balance of a Molecule is invalid"""
    atom: Atom

    def __init__(self, atom: Atom):
        super().__init__()
        self.atom = atom

    def __str__(self):
        return (f"Invalid, Atom ({self.atom.element.symbol}) w/ index {self.atom.idx} "
                f"has more bonds than its {self.atom.element.valences} available valences")


class FunctionalGroup:
    """A functional group in a molecule

    Instance Attributes:
        - group_type: The type of functional group.
        - priority: The IUPAC priority number (higher = higher priority for suffix selection).
        - suffix: The suffix used in IUPAC naming.
        - prefix: The prefix used in IUPAC naming.
        - attachment_atom: The index of the atom where this group attaches to the parent chain.
        - locant: The position of this group on the parent chain (set during numbering).
        - atoms: A list of atom indices that make up this functional group (for detection).
        - is_principal: Whether this group is the principal characteristic group (becomes suffix).
    """
    _groups: dict[str, dict]

    group_type: str
    priority: int
    suffix: str
    prefix: str
    attachment_atom: int
    locant: int | None
    atoms: list[int]
    is_principal: bool

    @classmethod
    def load_data(cls, file_path: str) -> None:
        """Loads the periodic element data needed to set up the private element classes"""
        import json
        with open(file_path) as f:
            cls._groups = json.load(f)

    @classmethod
    def get_group_data(cls, group_type: str) -> dict:
        """Get functional group data by functional group type, from external JSON file.

        Paramaters:
            - group_type:
                The functional group whose data has been queried.
        """
        if not cls._groups:
            cls.load_data("func_groups.json")

        if group_type not in cls._groups:
            raise ValueError(f"Unknown functional group: {group_type}")

        return cls._groups[group_type]

    def __init__(self, group_type: str, attachment_atom: int, atoms: list[int]):
        data = self.get_group_data(group_type)

        self.group_type = group_type
        self.attachment_atom = attachment_atom
        self.atoms = atoms
        self.priority = data["priority"]
        self.suffix = data["suffix"]
        self.prefix = data["prefix"]
        self.locant = None
        self.is_principal = False


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
        """Gets the exisiting instance of a given element, raises an error if the element is not real.

        Parameters:
            - symbol: string
                Chemical Symbol from the periodic table that represents the element in question.
        """
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
        - is_aromatic: A boolean that notes if the atom is a part of an aromatic ring, default is False.
        - charge: Stores the charge of this atom, set to 0 if unspecified.
        - isotope: Stores the proton number if this atom is an unusual isotope, None otherwise.
        - ex_h: Stores the number of hydrogens implicitly if specified, None otherwise.
    """
    element: _Element
    idx: int
    bonds: dict[int, int]  # Atom index -> Order
    is_aromatic: bool
    charge: int
    isotope: Optional[int]
    exp_h: Optional[int]

    def __init__(self, element: _Element, idx: int, is_aromatic: Optional[bool] = False,
                 charge: Optional[int] = 0, isotope: Optional[int] = None, exp_h: Optional[int] = None):
        self.element = element
        self.idx = idx
        self.bonds = {}
        self.is_aromatic = is_aromatic
        self.charge = charge
        self.isotope = isotope
        self.exp_h = exp_h

    def add_bond(self, other_idx: int, order: int = 1) -> None:
        """Adds a bond to another atom"""
        if other_idx not in self.bonds:
            self.bonds[other_idx] = order

    def get_hydrogen_count(self) -> int:
        """Returns the number of Hydrogens attached to this atom.

        Precondition:
            - Relies on the bonds to be checmially valid.
        """
        if self.exp_h:
            return self.exp_h

        bond_sum = sum(self.bonds.values())

        for valence in self.element.valences:
            h_count = valence - bond_sum - self.charge
            if h_count >= 0:
                return h_count

        raise ValenceError(self)


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
                result += f"   {self._atoms[bond].element.symbol}{bond} with order {atom.bonds[bond]}\n"
            result += f"{' & is aromatic' if atom.is_aromatic else ''}"\
                      f"{' & has charge ' + str(atom.charge) if atom.charge != 0 else ''}"\
                      f"{' & is an isotope' if atom.isotope else ''}\n"
        return result

    def add_atom(self, symbol: str, is_aromatic: bool = False, charge: Optional[int] = 0,
                 isotope: Optional[int] = None, exp_h: Optional[int] = None) -> int:
        """Adds atom to the Molecule object and returns its index.

        Parameters:
            - symbol: string
                Chemical symbol of the element of atom to be added.
            - is_aromatic: boolean
                The Aromaticity of the atom.
            - charge: integer
                Stores the charge of this atom, set to 0 if unspecified.
            - isotope: integer | None
                Stores the proton number if this atom is an unusual isotope, None otherwise.
            - ex_h: integer | None
                Stores the number of hydrogens implicitly if specified, None otherwise.
        """
        idx = len(self._atoms)
        atom = Atom(_Element.get(symbol), idx, is_aromatic, charge, isotope, exp_h)
        self._atoms.append(atom)
        return idx

    def add_bonds(self, idx1: int, idx2: int, order: int = 1):
        """Adds an edge between the two atoms referenced by index, with a weight based on bond order

        Parameters:
            - idx1: integer
                Index of the first Atom to be bonded.
            - idx2: integer
                Index of the second Atom to be bonded.
            - order: integer (default 1)
                The chemical order of the bond.

        Exceptions:
            - idx1 != idx2
            - idx1 and idx2 exist
            - the bond doesn't exist yet.
        """
        if not (0 <= idx1 < len(self._atoms) and 0 <= idx2 < len(self._atoms)):
            raise IndexError

        if idx2 in self._atoms[idx1].bonds:
            raise ValueError(f"Bond already exists between {idx1} and {idx2}")

        if idx1 == idx2:
            raise ValueError(f"Cannot bond atom to itself: {idx1}")

        self._atoms[idx1].add_bond(idx2, order)
        self._atoms[idx2].add_bond(idx1, order)

    def get_atom(self, idx: int) -> Atom:
        """Returns atom object based on index."""
        return self._atoms[idx]

    def valence_validate(self) -> bool:
        """Returns true if the molecule's component atoms do not violate their valences.
        Raises an exception otherwise"""
        for atom in self._atoms:
            _ = atom.get_hydrogen_count()

        return True
