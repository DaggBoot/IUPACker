from entities import _Element, Molecule


def main():
    # Load periodic table data
    _Element.load_data("Periodic.json")

    # Example molecules to test the parser
    molecules = [
        "CH3N(CH3)COOH",
        "CH3CH(OH)CH3",
        "CH3CH=O"
    ]

    for formula in molecules:
        print(f"Parsing molecule: {formula}")
        mol = Molecule(formula)
        print(mol)
        print("-" * 40)


if __name__ == "__main__":
    main()
