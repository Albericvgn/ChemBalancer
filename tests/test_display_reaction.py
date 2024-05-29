from chembalancer.chembalancer import display_reaction

def test_display_reaction():
    # Test case 1: Combustion of methane (CH4)
    reactants_1 = [(1, 'CH4'), (2, 'O2')]
    products_1 = [(1, 'CO2'), (2, 'H2O')]
    expected_output_1 = "CH4 + 2 O2 → CO2 + 2 H2O"
    assert display_reaction(reactants_1, products_1) == expected_output_1

    # Test case 2: Reaction with multiple reactants and products
    reactants_2 = [(2, 'CO2'), (3, 'H2O')]
    products_2 = [(2, 'C2H5OH'), (4, 'CO2')]
    expected_output_2 = "2 CO2 + 3 H2O → 2 C2H5OH + 4 CO2"
    assert display_reaction(reactants_2, products_2) == expected_output_2

    # Test case 3: Reaction with a single reactant and multiple products
    reactants_3 = [(1, 'CH4')]
    products_3 = [(2, 'CO2'), (2, 'H2O')]
    expected_output_3 = "CH4 → 2 CO2 + 2 H2O"
    assert display_reaction(reactants_3, products_3) == expected_output_3

    # Test case 4: Reaction with only one reactant and one product
    reactants_4 = [(1, 'H2')]
    products_4 = [(1, 'H2O')]
    expected_output_4 = "H2 → H2O"
    assert display_reaction(reactants_4, products_4) == expected_output_4

    # Test case 5: Reaction with reactants but no products
    reactants_5 = [(2, 'O2')]
    products_5 = []
    try:
        display_reaction(reactants_5, products_5)
    except ValueError as e:
        assert str(e) == "Both reactants and products need at least one component."

    # Test case 6: Reaction with products but no reactants
    reactants_6 = []
    products_6 = [(1, 'H2O')]
    try:
        display_reaction(reactants_6, products_6)
    except ValueError as e:
        assert str(e) == "Both reactants and products need at least one component."

    # Test case 7: Reaction with invalid component format in reactants
    reactants_7 = [(2, 'O2'), 'H2']
    products_7 = [(1, 'H2O')]
    try:
        display_reaction(reactants_7, products_7)
    except ValueError as e:
        assert str(e) == f"Invalid component format: {('H2',)}. Expected a tuple of (coefficient, molecule)."

    # Test case 8: Reaction with invalid component format in products
    reactants_8 = [(2, 'O2')]
    products_8 = [(1, 'H2O'), (2,)]  # Missing coefficient in 'H2O'
    try:
        display_reaction(reactants_8, products_8)
    except ValueError as e:
        assert str(e) == f"Invalid component format: {(2,)}. Expected a tuple of (coefficient, molecule)."

    # Test case 9: Reaction with an empty reactant list
    reactants_9 = []
    products_9 = [(1, 'H2O')]
    try:
        display_reaction(reactants_9, products_9)
    except ValueError as e:
        assert str(e) == "Both reactants and products need at least one component."

    # Test case 10: Reaction with an empty product list
    reactants_10 = [(2, 'O2')]
    products_10 = []
    try:
        display_reaction(reactants_10, products_10)
    except ValueError as e:
        assert str(e) == "Both reactants and products need at least one component."
