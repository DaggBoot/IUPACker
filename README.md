# (Work In Progress) IUPACker: IUPAC Chemical Name Generator

A Python project that generates systematic **IUPAC names** from molecular graph representations of hydrocarbons.

The goal of this project is to explore how chemical nomenclature rules can be implemented algorithmically by representing molecules as graphs and applying traversal and parsing algorithms.

---

## Overview

Organic molecules can be represented as **graphs**, where:

- **Atoms** are vertices
- **Bonds** are edges

Using this representation, the program analyzes the structure of a molecule and determines its **systematic IUPAC name** by applying rules such as:

- Identifying the **longest carbon chain**
- Detecting **substituent groups**
- Determining **substituent positions**
- Constructing the correct **IUPAC name format**

This project focuses primarily on **hydrocarbon molecules**.

---

## Features

- Represents molecules as **graph structures** using `Atom` and `Molecule` classes.
- Supports parsing of:
  - Carbon chains
  - Oxygen, Nitrogen, Sulfur atoms
  - Halogens: F, Cl, Br, I
- Handles branching using parentheses `( )`.
- Automatically assigns **bond orders** from `-`, `=`, `#`.

---
## Example Usage
```
python
from entities import _Element, Molecule

# Load periodic element data
_Element.load_data("Periodic.json")

# Create a molecule
m = Molecule("CH3CH(CH3)CH2CHBrCOOH")

# Print molecule structure
print(m)
```
### Output:
```
CH3CH(CH3)CH2CHBrCOOH 
C Connected to 
   C with order 1 
    Has 3 hydrogens
 
C Connected to 
   C with order 1 
   C with order 1 
   C with order 1 
    Has 1 hydrogens
 
C Connected to 
   C with order 1 
    Has 3 hydrogens
 
C Connected to 
   C with order 1 
   C with order 1 
    Has 2 hydrogens
 
C Connected to 
   C with order 1 
   Br with order 1 
    Has 1 hydrogens
 
Br Connected to 
   C with order 1 
    Has 0 hydrogens
 
C Connected to 
   O with order 1 
    Has 0 hydrogens
 
O Connected to 
   C with order 1 
   O with order 1 
    Has 0 hydrogens
 
O Connected to 
   O with order 1 
    Has 1 hydrogens
```
Example molecular graph:
```
    CH3    Br
    |      |
CH3-CH-CH2-CH-C=O
              |
              OH
```
---

## Project Structure


iupacker/
│
├── entities.py         # Core classes: _Element, Atom, Molecule
├── Periodic.json       # Element data (atomic number, valence)
└── main.py             # Optional entry point / tests


---

## Technologies Used

- Python
- Graph algorithms
- Recursive traversal
- Parsing and tokenization

---

## Current Status

**Work in Progress**

The current implementation focuses on:

- Hydrocarbon chain detection
- Basic substituent identification
- Core graph representation

Future work includes:

- Support for more functional groups
- More complex branching structures
- Improved input parsing

---

## Motivation

This project was created to explore the intersection of:

- **chemistry**
- **graph theory**
- **algorithm design**

Systematic chemical naming provides an interesting challenge because it requires translating a set of **formal scientific rules into algorithms** operating on graph structures.

---

## Future Improvements

- Functional group detection
- Ring structures
- Visualization of molecular graphs

---

## Author

Shivanshu Vel Rajeev
