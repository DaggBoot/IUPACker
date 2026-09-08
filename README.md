# (Work In Progress) IUPACker: IUPAC Chemical Name Generator

A Python project that generates systematic **IUPAC names** from **SMILES** strings, by parsing them into a molecular graph and applying IUPAC nomenclature rules algorithmically.

The goal of this project is to explore how chemical nomenclature rules — parent selection, ring perception, functional group priority, numbering, and name assembly — can be implemented as graph algorithms rather than pattern-matched by hand.

---

## Overview

A SMILES string is parsed into a **molecular graph**, where atoms are vertices and bonds are edges. From there, the pipeline:

1. **Detects rings** in the graph, including fused, bridged, and spiro systems.
2. **Detects functional groups** by matching structural patterns against the graph.
3. **Selects the parent structure** (a chain or a ring system) by applying the IUPAC seniority cascade — principal group presence, principal group count, senior element, ring-vs-chain seniority, chain length, and unsaturation — including deciding, for principal groups whose defining atom is carbon (carboxylic acid, aldehyde, nitrile, amide), whether every instance can be absorbed into the parent chain or whether the name must switch to a detached suffix form (`-carboxylic acid`, `-carbaldehyde`, etc.) applied uniformly across all instances.
4. **Identifies substituents** branching off the parent, recursively, including substituents linked through a bridging heteroatom (e.g. an ether oxygen), substituents that are themselves rings, and non-principal functional groups (halogens, a non-selected carbonyl, etc.) sitting on the parent or on any substituent.
5. **Numbers the parent** by choosing the locant direction that satisfies the numbering criteria in order (heteroatom locants, principal group locants, unsaturation locants, substituent locant set, alphabetical citation order), for both chains and simple monocyclic rings, and numbers substituent chains/rings so the point of attachment gets the lowest possible locant.
6. **Assembles the final name string**: alphabetized, grouped, and multiplied substituent/non-principal-group prefixes; parenthesized composite substituent names; the parent stem and unsaturation infix; and the correct principal-group suffix — while omitting locants wherever a position is structurally forced rather than genuinely competing with an alternative (e.g. `ethene` not `eth-1-ene`, `cyclohexanecarboxylic acid` not `cyclohexane-1-carboxylic acid`, `cyanomethyl` not `1-cyanomethyl`).

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
     │  _find_parent_chain / _find_parent_chain_no_p: parent selection,
     │      including chain-absorption vs. detached-suffix resolution
     │  _construct_substituent / _name_substituent: recursive substituent
     │      detection and recursive name-string assembly
     │  _number_chain / _number_ring / _number_substituent_ring: locant
     │      direction selection
     │  _name_parent: final name string assembly
     ▼
IUPAC name string (e.g. "3-methylheptane", "cyclohexanecarboxylic acid",
                         "heptane-1,4,7-tricarboxylic acid")
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
- Carbon-defined suffix classes (carboxylic acid, aldehyde, nitrile, amide) carry both a chain-absorbed suffix (`-oic acid`, `-al`, `-nitrile`, `-amide`) and a detached suffix (`-carboxylic acid`, `-carbaldehyde`, `-carbonitrile`, `-carboxamide`) in `patterns.py`, since the correct form depends on whether the group's carbon can be counted as part of the parent.

### Parent structure selection (`namer.py`)
- Applies the IUPAC seniority cascade for choosing the parent: contains the principal characteristic group → maximum number of principal groups → senior element → ring senior to chain of the same element → chain length → unsaturation.
- Handles chains where a principal-group atom sits mid-chain rather than at an end (e.g. propan-2-ol), via a subtree-pivot-merge search.
- Correctly restricts single-atom candidates to carbon (per the rule that a lone heteroatom, like a fully-substituted sulfonic acid sulfur, is never a valid standalone parent on its own).
- For carbon-defined principal groups with more than two instances, or with an instance that is structurally a branch relative to the others, determines whether a single chain can absorb every instance and, if not, forces all instances into the detached-suffix form uniformly (e.g. a tricarboxylic acid becomes `heptane-1,4,7-tricarboxylic acid` rather than mixing an absorbed `-dioic acid` with a stray `carboxy-` prefix).

### Substituent detection and naming (`namer.py`)
- Recursively finds and structures every substituent branching off the parent (and off nested substituents), excluding ring atoms and atoms already claimed by any detected functional group (center and matched atoms alike).
- Correctly absorbs a substituent linked through a bridging heteroatom with no same-element neighbour of its own (e.g. an ether oxygen) into a single combined substituent, rather than splitting it into a heteroatom "substituent" with the real alkyl chain incorrectly nested underneath — this generalizes to multi-atom bridges (e.g. `-O-CH2CH2-O-`) automatically.
- Recursively renders each substituent (and nested substituent) to a name string, including ring substituents (e.g. `cyclopropylmethyl`), non-principal functional groups sitting on a substituent (e.g. `cyanomethyl`), and correct enclosing marks (`()`, `[]`, `{}`, cycled by nesting depth) whenever a substituent name is itself composite.
- Omits a substituent's own attachment or internal locant wherever it is structurally forced (a one-atom substituent, an attachment point on an otherwise-unadorned ring), and correctly numbers substituent chains/rings so the attachment point gets the lowest possible locant when it is not forced to be first.

### Numbering (`namer.py`)
- `_number_chain` and `_number_ring` select the correct locant direction for a parent chain or monocyclic ring by applying, in order: heteroatom locants, principal-group locants, unsaturation locants, the full substituent locant set, and alphabetical order of citation for ties.
- `_number_substituent_ring` numbers a substituent ring so its point of attachment is locant 1, per the substituent-numbering rule, then resolves direction by the same profile-based criteria among whatever else is on the ring.
- Locant assignment for substituents uses the parent-chain attachment point, not the substituent's own internal atom indices.

### Name assembly (`namer.py`)
- Collects non-principal functional-group prefixes and named substituents into one alphabetized, grouped list, applying the correct multiplying prefix (`di-`/`tri-`/... for simple substituents, `bis-`/`tris-`/... for composite ones) and joining locants with commas.
- Applies the chain-absorbed or detached suffix form for the principal group depending on parent selection's absorb/detach decision, with the correct multiplying prefix and locant set (or their omission) for each.
- Omits locants specifically where the position is structurally forced rather than merely convenient — a two-atom unsaturated chain (`ethene`), a sole substituent on an otherwise plain ring, a sole detached-suffix instance on an otherwise plain ring, and a one-atom substituent's internal group — while always retaining them wherever a genuine numbering choice exists.

---

## Example Usage

```python
from smiles_parser import SMILESParser
from namer import IUPACker

parser = SMILESParser()
molecule = parser.parse("CC(C)CCCCC")

namer = IUPACker(molecule)
name = namer.generate()

print(name)             # "2-methylheptane"
print(namer._parent_chain)   # the parent, correctly directed for lowest locants
print(namer._subs)           # structured, named substituent tree
```

---

## Project Structure

```
iupacker/
│
├── entities.py         # _Element, Atom, Molecule, Ring, MotifPattern/BondReq/AtomCond
├── smiles_parser.py     # SMILESParser: SMILES -> Molecule
├── motif_engine.py      # MotifEngine: functional group pattern matching
├── namer.py             # IUPACker: parent selection, substituents, numbering, name assembly
├── patterns.py          # MotifPattern definitions for supported functional groups
└── periodic.json        # Element data (valences)
```

---

## Current Status

**Work in Progress.**

**Implemented:**
- SMILES parsing, including rings and branches
- Ring detection, including fused/bridged/spiro classification
- Functional group pattern matching, including multiple groups on one atom, and chain-absorbed vs. detached suffix forms for carbon-defined principal groups
- Parent chain and parent ring-system selection via the full seniority cascade, including the absorb-vs-detach decision for carbon-defined principal groups present in more than two instances or in a branching arrangement
- Recursive substituent detection and full recursive name-string generation, including bridging-heteroatom substituents, ring substituents, and non-principal functional groups on any substituent
- Locant-direction selection for a parent **chain** and for simple **monocyclic rings** (heteroatom, principal group, unsaturation, substituent-set, and citation-order criteria)
- Substituent-level numbering (attachment point gets the lowest possible locant, for both chain and ring substituents)
- Full name-string assembly: multiplying prefixes (`di-`/`tri-`/`bis-`/`tris-`), alphabetized and grouped citation order, enclosing-mark nesting, and locant omission for structurally forced positions

**Not yet implemented / known gaps:**
- Aromatic ring naming (benzene-derived and other aromatic systems)
- Numbering and nomenclature for heteroatoms sitting *inside* a ring or chain skeleton (replacement/Hantzsch-Widman nomenclature) — the current namer does not treat these as anything other than "not yet handled"
- Von Baeyer (bicyclic+), bridged, and spiro ring **numbering and naming** (ring *detection* and fused-system grouping exist; the numbering/naming rules for these systems are not implemented)
- "Indicated hydrogen" locants and "hydro-"/"dehydro-" prefix locants
- Detached-suffix handling for principal groups present in more than two instances where a mixed absorbed/detached treatment might otherwise be considered (currently resolved by forcing all instances to detached mode uniformly, which is correct per IUPAC rules, but the underlying candidate search for this case has only been stress-tested on a small set of examples)
- Retained/common names (e.g. formic acid, acetone) — only fully systematic names are produced
- A formal automated regression test suite — correctness has been validated against a hand-picked set of SMILES strings run manually, not via a repeatable test harness, so there is real risk of regressions as new features are added

---

## Motivation

This project was created to explore the intersection of:
- **chemistry**
- **graph theory**
- **algorithm design**

---

## Future Improvements
- Aromatic ring naming
- Heteroatom-in-ring/chain (replacement) nomenclature
- Von Baeyer/bridged/spiro ring numbering and naming
- A proper automated regression test suite covering the functional-group, parent-selection, substituent, numbering, and name-assembly edge cases already discovered during development
- Broader functional group coverage in `patterns.py`
- Visualization of molecular graphs

---

## Author
Shivanshu Vel Rajeev
