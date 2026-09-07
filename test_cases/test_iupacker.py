"""
Test Suite for the iupacker.py module.
"""
import sys
import os
import pytest
from iupacker import generate_name

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.mark.parametrize("smiles, expected", [
    ("CC", "ethane"),
    ("CCC", "propane"),
    ("CCCC", "butane"),
    ("CCCCC", "pentane"),
    ("CCCCCC", "hexane"),
    ("CCCCCCC", "heptane"),
    ("CCCCCCCC", "octane"),
    ("CCCCCCCCC", "nonane"),
    ("CCCCCCCCCC", "decane"),
])
def test_alkanes(smiles, expected):
    assert generate_name(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ("C=C", "eth-1-ene"),
    ("CC=C", "prop-1-ene"),
    ("C=CC", "prop-1-ene"),
    ("CC=CC", "but-2-ene"),
])
def test_alkenes(smiles, expected):
    assert generate_name(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ("C#C", "eth-1-yne"),
    ("CC#C", "prop-1-yne"),
    ("C#CC", "prop-1-yne"),
    ("CC#CC", "but-2-yne"),
])
def test_alkynes(smiles, expected):
    assert generate_name(smiles) == expected
