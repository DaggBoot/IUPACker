# (Work In Progress) IUPACker: IUPAC Chemical Name Generator

A Python project that generates systematic **IUPAC names** from **SMILES** strings, by parsing them into a molecular graph and applying IUPAC nomenclature rules algorithmically.

The goal of this project is to explore how chemical nomenclature rules — parent selection, ring perception, functional group priority, numbering — can be implemented as graph algorithms rather than pattern-matched by hand.

---

## Overview

A SMILES string is parsed into a **molecular graph**, where atoms are vertices and bonds are edges. From there, the pipeline:

1. **Detects rings** in the graph, including fused, bridged, and spiro systems.
2. **Detects functional groups** by matching structural patterns against the graph.
3. **Selects the parent structure** (a chain or a ring system) by applying the IUPAC seniority cascade — principal group presence, principal group count, senior element, ring-vs-chain seniority, chain length, and unsaturation.
4. **Identifies substituents** branching off the parent, recursively, including substituents linked through a bridging heteroatom (e.g. an ether oxygen).
5. **Numbers the parent** by choosing the locant direction that satisfies the numbering criteria in order (heteroatom locants, principal group locants, unsaturation locants, substituent locant set, alphabetical citation order).

Step 6 — assembling the final name string (multiplying prefixes, alphabetized citation, suffixes) — is not yet implemented; see **Current Status** below.

---

## Architecture

```
SMILES string
     │
     ▼
SMILESParser            (smiles_parser.py)
     │  tokenizes, builds atom/bond graph, tracks ring-closure bonds
     ▼
Molecule                (entities.py)
     │  Molecule.compute_rings(): finds an independent ring basis via
     │  shortest-path-excluding-each-bond candidates + Gaussian elimination
     │  over GF(2); groups fused/bridged/spiro rings into ring systems
     ▼
MotifEngine              (motif_engine.py)
     │  matches MotifPatterns (functional groups) against the molecule;
     │  a center atom with more qualifying neighbours than a pattern
     │  requires yields multiple separate matches (e.g. a geminal diol)
     ▼
IUPACker                 (namer.py)
     │  _find_parent_chain / _find_parent_chain_no_p: parent selection
     │  _name_substituent: recursive substituent detection
     │  _number_chain: locant-direction selection
     ▼
(name string assembly -- not yet implemented)
```

---

## Features

### Parsing (`smiles_parser.py`, `entities.py`)
- Tokenizes and parses SMILES, including branches, ring-closure digits (single and `%NN` two-digit form), bracket atoms (charge, isotope, explicit hydrogens), aromatic lowercase atoms, and bond-order symbols (`-`, `=`, `#`, `$`, `:`).
- Builds a `Molecule` graph of `Atom` objects with bond dictionaries; validates valence.

### Ring detection (`entities.py`)
- Finds a linearly independent set of rings spanning the molecule's full cycle space (one ring per unit of cyclomatic complexity), using every bond in the molecule as a candidate source and Gaussian elimination over GF(2) to select an independent basis — correctly handles cases where a naive "shortest path per closure bond" approach would under-count independent rings (e.g. cubane, bridged bicyclics).
- Classifies which rings are fused, bridged, or spiro to one another (`Ring.fused`), and groups a whole connected ring system into one parent candidate rather than treating each component ring separately.

### Functional group detection (`motif_engine.py`)
- Declarative `MotifPattern`/`BondReq`/`AtomCond` definitions describe a functional group as a center atom plus required neighbouring bonds and conditions.
- Correctly produces multiple separate matches when a center atom has more qualifying neighbours than a requirement's count (e.g. two `-OH` groups on the same carbon yield two matches, not zero or one merged match).
- Resolves overlapping matches by pattern priority.

### Parent structure selection (`namer.py`)
- Applies the IUPAC seniority cascade for choosing the parent: contains the principal characteristic group → maximum number of principal groups → senior element → ring senior to chain of the same element → chain length → unsaturation.
- Handles chains where a principal-group atom sits mid-chain rather than at an end (e.g. propan-2-ol), via a subtree-pivot-merge search.
- Correctly restricts single-atom candidates to carbon (per the rule that a lone heteroatom, like a fully-substituted sulfonic acid sulfur, is never a valid standalone parent on its own).

### Substituent detection (`namer.py`)
- Recursively finds and structures every substituent branching off the parent (and off nested substituents), excluding ring atoms and atoms already claimed by a detected functional group.
- Correctly absorbs a substituent linked through a bridging heteroatom with no same-element neighbour of its own (e.g. an ether oxygen) into a single combined substituent, rather than splitting it into a heteroatom "substituent" with the real alkyl chain incorrectly nested underneath — this generalizes to multi-atom bridges (e.g. `-O-CH2CH2-O-`) automatically.

### Numbering (`namer.py`)
- `_number_chain` selects the correct locant direction for a parent chain by applying, in order: heteroatom locants, principal-group locants, unsaturation locants, the full substituent locant set, and alphabetical order of citation for ties.
- Locant assignment for substituents uses the parent-chain attachment point, not the substituent's own internal atom indices.

---

## Example Usage

```python
from smiles_parser import SMILESParser
from namer import IUPACker

parser = SMILESParser()
molecule = parser.parse("CC(C)CCCCC")

namer = IUPACker(molecule)
namer.generate()

print(namer._parent_chain)   # the parent, correctly directed for lowest locants
print(namer._subs)           # structured substituent tree
```

Full IUPAC name strings (e.g. `"2-methylheptane"`) are not produced yet — `generate()` currently returns/exposes the resolved parent chain and structured substituent data, not final text. See **Current Status**.

---

## Project Structure

```
iupacker/
│
├── entities.py         # _Element, Atom, Molecule, Ring, MotifPattern/BondReq/AtomCond
├── smiles_parser.py     # SMILESParser: SMILES -> Molecule
├── motif_engine.py      # MotifEngine: functional group pattern matching
├── namer.py             # IUPACker: parent selection, substituents, numbering
├── patterns.py          # MotifPattern definitions for supported functional groups
└── periodic.json        # Element data (valences)
```

---

## Technologies Used
- Python
- Graph algorithms (cycle detection, Gaussian elimination over GF(2), DFS/BFS traversal)
- Recursive parsing and structure-building

---

## Current Status

**Work in Progress.**

**Implemented:**
- SMILES parsing, including rings and branches
- Ring detection, including fused/bridged/spiro classification
- Functional group pattern matching, including multiple groups on one atom
- Parent chain and parent ring-system selection via the full seniority cascade
- Recursive substituent detection, including bridging-heteroatom substituents
- Locant-direction selection for a parent **chain** (heteroatom, principal group, unsaturation, substituent-set, and citation-order criteria)

**Not yet implemented:**
- Numbering for **rings** (monocyclic and fused/bridged systems need separate numbering algorithms from chains — an open design question currently being worked through)
- Von Baeyer (bicyclic+) and spiro nomenclature and numbering
- "Indicated hydrogen" locants and "hydro-" prefix locants
- Substituent name **text** generation (currently only a placeholder exists, used internally for alphabetical-order tiebreaking)
- Multiplying prefixes (`di-`, `tri-`, `bis-`, `tris-`) for repeated identical substituents
- Final name string assembly (locants + prefixes + parent name + suffix)
- Hantzsch-Widman and retained names for heterocycles

---

## Motivation

This project was created to explore the intersection of:
- **chemistry**
- **graph theory**
- **algorithm design**

---

## Future Improvements
- Ring numbering and von Baeyer/spiro nomenclature
- Full name-text generation and assembly
- Broader functional group coverage in `patterns.py`
- Visualization of molecular graphs

---

## Author
Shivanshu Vel Rajeev
