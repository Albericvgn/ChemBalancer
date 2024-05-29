def balance_chemical_equation(reactant_smiles, product_smiles):
    reactant_counts = [count_atoms(smiles) for smiles in reactant_smiles]
    product_counts = [count_atoms(smiles) for smiles in product_smiles]

    reactant_elements = set(sum([list(counts.keys()) for counts in reactant_counts], []))
    product_elements = set(sum([list(counts.keys()) for counts in product_counts], []))

    if reactant_elements != product_elements:
        missing_in_products = reactant_elements - product_elements
        missing_in_reactants = product_elements - reactant_elements
        error_message = "Element mismatch found: "
        if missing_in_products:
            error_message += f"Elements {missing_in_products} are in reactants but not in products. "
        if missing_in_reactants:
            error_message += f"Elements {missing_in_reactants} are in products but not in reactants."
        raise ValueError(error_message)

    elements = sorted(reactant_elements.union(product_elements))
    A_reactants = setup_matrix(elements, reactant_counts)
    A_products = setup_matrix(elements, product_counts)
    A = np.concatenate([A_reactants, -A_products], axis=1)

    integer_coefficients = solve_ilp(A)
    if integer_coefficients is None or not integer_coefficients:
        raise ValueError("Failed to solve the balance equation. The system may be underdetermined or inconsistent.")

    reactant_coeffs = integer_coefficients[:len(reactant_smiles)]
    product_coeffs = integer_coefficients[len(reactant_smiles):]

    reactant_data = [(coeff, get_molecular_formula(smiles)) for coeff, smiles in zip(reactant_coeffs, reactant_smiles)]
    product_data = [(coeff, get_molecular_formula(smiles)) for coeff, smiles in zip(product_coeffs, product_smiles)]

    return reactant_data, product_data


    
def setup_matrix(elements, counts):
    # Assuming counts is a list of dictionaries where each dictionary is the atomic count for a molecule
    matrix = []
    for count in counts:
        row = [count.get(element, 0) for element in elements]
        matrix.append(row)
    
    # Ensure matrix is 2D
    matrix = np.array(matrix)
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)  # Reshape to 2D if it's inadvertently 1D

    return matrix