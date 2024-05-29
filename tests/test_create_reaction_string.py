from chembalancer.chembalancer import create_reaction_string


def test_create_reaction_string():
    # Test case 1: Reaction with simple reactants and products
    reactants_1 = ['CH4', 'O2']
    products_1 = ['CO2', 'H2O']
    expected_output_1 = "CH4.O2>>CO2.H2O"
    assert create_reaction_string(reactants_1, products_1) == expected_output_1

    # Test case 2: Reaction with more complex reactants and products
    reactants_2 = ['C2H6', 'O2']
    products_2 = ['CO2', 'H2O', 'C2H5OH']
    expected_output_2 = "C2H6.O2>>CO2.H2O.C2H5OH"
    assert create_reaction_string(reactants_2, products_2) == expected_output_2

    # Test case 3: Reaction with a single reactant and a single product
    reactants_3 = ['H2']
    products_3 = ['H2O']
    expected_output_3 = "H2>>H2O"
    assert create_reaction_string(reactants_3, products_3) == expected_output_3

    # Test case 4: Reaction with an empty list of reactants
    reactants_4 = []
    products_4 = ['H2O']
    expected_output_4 = ">>H2O"
    assert create_reaction_string(reactants_4, products_4) == expected_output_4

    # Test case 5: Reaction with an empty list of products
    reactants_5 = ['H2']
    products_5 = []
    expected_output_5 = "H2>>"
    assert create_reaction_string(reactants_5, products_5) == expected_output_5

    # Test case 6: Reaction with reactants and products as empty lists
    reactants_6 = []
    products_6 = []
    expected_output_6 = ">>"
    assert create_reaction_string(reactants_6, products_6) == expected_output_6

    # Test case 7: Reaction with invalid input (integers instead of strings)
    reactants_7 = [1, 2]
    products_7 = [3, 4]
    try:
        create_reaction_string(reactants_7, products_7)
    except TypeError as e:
        assert str(e) == "sequence item 0: expected str instance, int found"
