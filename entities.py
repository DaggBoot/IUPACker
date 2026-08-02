"""
Contains entities that represent a molecule.
"""
from __future__ import annotations

from collections import deque, defaultdict
from typing import Optional, Any
from dataclasses import dataclass, field


class ValenceError(Exception):
    """Custom Exception for when the valence balance of a Molecule is invalid"""
    atom: Atom

    def __init__(self, atom: Atom):
        super().__init__()
        self.atom = atom

    def __str__(self):
        return (f"Invalid, Atom ({self.atom.element.symbol}) w/ index {self.atom.idx} "
                f"has more bonds than its {self.atom.element.valences} available valences")


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

        for valence in self.element.valences:
            h_count = valence - sum(self.bonds.values()) - self.charge
            if h_count >= 0:
                return h_count

        raise ValenceError(self)

    def bond_sum(self) -> int:
        """Returns the sum of the bond orders"""
        return sum(self.bonds.values()) + self.get_hydrogen_count()


class Ring:
    """A ring in a molecule.

    Instance Attributes:
        - atoms: list of the integer idxs of the atoms in the ring
        - bonds: a frozenset of frozenset, encoding all the conds between the atoms in atoms.
    """
    atoms: list[int]

    def __init__(self, atoms: list[int]):
        self.atoms = atoms

    @property
    def bonds(self) -> frozenset[frozenset[int]]:
        """Returns a frozenset of all the bonds that make up this cycle.
        """
        n = len(self.atoms)
        return frozenset(
            frozenset((self.atoms[i], self.atoms[(i + 1) % n])) for i in range(n)
        )


class Molecule:
    """ The graph representation for an Organic Molecule

    Instance Attributes:
        - _atoms: A collection of Atom objects (vertices) that form a graph which represents the molecule.
        - atom_rings: A collection of Ring objects that form a ring within this Molecule.
    """

    _atoms: list[Atom]
    atom_rings: list[Ring]

    def __init__(self):
        self._atoms = []
        self.atom_rings = []

    def __str__(self) -> str:
        result = ""
        for atom in self._atoms:
            result += f"{atom.element.symbol}{atom.idx} Connected to \n"
            for bond in atom.bonds:
                result += f"   {self._atoms[bond].element.symbol}{bond} with order {atom.bonds[bond]}\n"
            result += f"{' & is aromatic' if atom.is_aromatic else ''}"\
                      f"{' & has charge ' + str(atom.charge) if atom.charge != 0 else ''}"\
                      f"{' & is an isotope' if atom.isotope else ''}\n"
        return result

    def __iter__(self):
        """Iterate over all atoms in a Molecule"""
        return iter(self._atoms)

    def __len__(self) -> int:
        """Returns the integer number of atoms in the Molecule."""
        return len(self._atoms)

    def __getitem__(self, idx: int) -> Atom:
        """Returns Atom corresponding to input index.

        Parameters:
            - idx: integer
                Assumed to be an existing integer index.
        """
        return self._atoms[idx]

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

    def valence_validate(self) -> bool:
        """Returns true if the molecule's component atoms do not violate their valences.
        Raises an exception otherwise"""
        for atom in self._atoms:
            _ = atom.get_hydrogen_count()

        return True

    def compute_rings(self, closure_bonds: list[tuple[int, int]]) -> list[Ring]:
        """Returns a list of all rings in the molecule through Ring objects, and updates the self.atom_rings attribute

        Parameters:
            closure_bonds: list of tuples of integer pairs
                - contains a connecting bond in a ring.
        """
        if not closure_bonds:
            self.atom_rings = []
            return self.atom_rings

        all_bonds = set()
        for atom in self:
            for neighbour_idx in atom.bonds:
                all_bonds.add(frozenset((atom.idx, neighbour_idx)))

        candidates = []
        for a, b in closure_bonds:
            paths = self._shortest_path_excluding(a, b, frozenset((a, b)))
            if paths:
                candidates.extend([Ring(path) for path in paths])

        target_rank = len(all_bonds) - len(self._atoms) + 1
        bond_index = {bond: i for i, bond in enumerate(all_bonds)}

        self.atom_rings = self._select_independent_rings(candidates, bond_index, target_rank)
        return self.atom_rings

    def _shortest_path_excluding(self, start: int, end: int, visited_bond: frozenset[int]) -> Optional[list[int]]:
        """Returns the shortest paths (as a list of lists of atom indices, start to end inclusive)
        between `start` and `end`, without using the bond `visited_bond`. Returns None otherwise"""

        distance = {start: 0}
        parents = defaultdict(list)

        queue = deque([start])

        while queue:
            curr = queue.popleft()

            for neighbour_idx in self._atoms[curr].bonds:
                if frozenset((curr, neighbour_idx)) == visited_bond:
                    continue

                new_distance = distance[curr] + 1

                if neighbour_idx not in distance:
                    distance[neighbour_idx] = new_distance
                    parents[neighbour_idx].append(curr)
                    queue.append(neighbour_idx)

                elif distance[neighbour_idx] == new_distance:
                    parents[neighbour_idx].append(curr)

        if end not in distance:
            return []

        paths = []

        self._backtrack_shortest_paths(start, end, parents, [end], paths)

        return paths

    def _backtrack_shortest_paths(self, start: int, node: int, parents: dict[int, list[int]],
                                  path: list[int], paths: list[list[int]]) -> None:
        """Recursively reconstructs all shortest paths."""
        if node == start:
            paths.append(list(reversed(path)))
            return

        for parent in parents[node]:
            self._backtrack_shortest_paths(start, parent, parents, path + [parent], paths)

    @staticmethod
    def _select_independent_rings(candidates: list[Ring], bond_index: dict[frozenset[int], int],
                                  target_rank: int) -> list[Ring]:
        """Returns a greedily selected linearly independent (in a vector space over a finite field of size 2)
        subset of candidates with the smallest rings first, stopping at `target_rank` independent rings.
        Uses Gaussian elimination, considering cycles (aka rings) as bit-vectors.
        """
        candidates = sorted(candidates, key=lambda x: len(x.atoms))
        basis: list[tuple[int, int]] = []
        selected = []

        for ring in candidates:
            if len(selected) == target_rank:
                break

            mask = 0
            for bond in ring.bonds:
                mask |= 1 << bond_index[bond]

            for pivot, pivot_mask in basis:
                if mask & (1 << pivot):
                    mask ^= pivot_mask

            if mask == 0:
                continue

            pivot = mask.bit_length() - 1
            basis.append((pivot, mask))
            selected.append(ring)

        return selected


@dataclass
class BondReq:
    """The required bond for the identification of a MotifPattern, from its central atom.

    Instance Attributes:
        - symbol: Chemical formula for the other atom in the bond.
        - order: The order of the bond.
        - count: The number of such bonds we allow. Defaults to 1.
        - conditions: The conditions to further filter valid atoms. Defaulted to empty.
        - future_req: BondRequirements for the atom in this bond
    """
    symbol: str
    order: float
    count: int = 1
    conditions: list[AtomCond] = field(default_factory=list)
    future_req: list[BondReq] = field(default_factory=list)


@dataclass
class AtomCond:
    """Conditions to enforce on an atom to accept it as a specific MotifPattern.

    Instance Attributes:
        - property: Any property that must have an enforced value.
            Those supported are:
                - "charge": The integer formal charge.
                - "aromatic": Boolean on aromaticity.
                - "has_h": Boolean for if it has a hydrogen bonded.
                - "bond_sum": The positive integer sum of bond orders.

        - value: The value that we wish to enforce onto the property under the operation described by oper.

        - oper: The operator that describes the relationship beteen the property and the value.
            Those supported are:
                - "eq": equals, i.e. the == operator.
                - "ne": not equals, i.e. the != operator.
                - "gt": greater than, i.e. the > operator.
                - "lt": less than, i.e. the < operator.
                - "ge": greater than or equal, i.e. the >= operator.
                - "le": less than or equal to, i.e. the <= operator.
    """
    property: str
    value: Any
    oper: str = "eq"


@dataclass
class MotifPattern:
    """A complete functional group pattern.

    Instance Attributes:
        - name: The name of a MotifPattern
        - priority: How it ranks in terms of other MotifPatterns (used in case of overlaps).
        - suffix: Used for the IUPACker chemical molecule namer, gives the suffix string for the name.
        - prefix: Used for the IUPACker chemical molecule namer, gives the prefix string for the name.

        - center_symbol: The symbol of the atom that this MotifPattern is based on.
        - center_conditions: The conditions that classify the center atom as part of the MotifPattern.
        - bonds: The list of bond requirments to classify the center as a part of the MotifPattern.
        - excludes: Other MotifPatterns it must ensure are not found. (Used in case of conflcits between patterns)
    """
    name: str
    priority: Optional[int]
    suffix: Optional[str]
    prefix: Optional[str]

    center_symbol: str
    center_conditions: list[AtomCond] = field(default_factory=list)
    bonds: list[BondReq] = field(default_factory=list)
    excludes: list[MotifPattern] = field(default_factory=list)
