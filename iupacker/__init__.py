"""
IUPACker - Chemical name generator from SMILES strings.
"""

from .assembly import generate_name, generate_name_from_mol
from .smiles_parser import SMILESParser
from .motif_engine import MotifEngine

__all__ = [
    'generate_name',
    'generate_name_from_mol',
    'SMILESParser',
    'MotifEngine',
]

__version__ = "0.4.0"
