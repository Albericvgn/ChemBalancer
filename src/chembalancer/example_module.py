"""Example module to get you started."""


def hello_smiles(smiles: str) -> str:
    """
    Return a greeting string that incorporates the given smiles.

    Parameters
    ----------
    smiles : str
        A text string representing a SMILES (Simplified
        Molecular Input Line Entry System) notation or any string.

    Returns
    -------
    str
        A greeting message incorporating the input smiles.

    Examples
    --------
    >>> hello_smiles("C(=O)O")
    'Hello, C(=O)O!'
    """
    return f"Hello, {smiles}!"


if __name__ == "__main__":
    print(hello_smiles("C(=O)O"))


"""DEBUT DE NOTRE CODE A NOUS
""""

def get_smiles_from_name(name):
    """Fetch the SMILES string of a molecule by its common name from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        return smiles
    else:
        return "No data found or error occurred."

def count_atoms(smiles):
    """ Count atoms in a SMILES string. """
    mol = Chem.MolFromSmiles(smiles)
    atom_counts = defaultdict(int)
    if mol:
        for atom in mol.GetAtoms():
            atom_counts[atom.GetSymbol()] += 1
    return dict(atom_counts)

def solve_ilp(A):
    """ Solve the integer linear programming problem to find stoichiometric coefficients. """
    num_vars = A.shape[1]
    prob = pulp.LpProblem("Balancing_Chemical_Equation", pulp.LpMinimize)
    x_vars = [pulp.LpVariable(f'x{i}', lowBound=1, cat='Integer') for i in range(num_vars)]
    prob += pulp.lpSum(x_vars)
    for i in range(A.shape[0]):
        prob += pulp.lpDot(A[i, :], x_vars) == 0
    solver = pulp.PULP_CBC_CMD(msg=False)
    prob.solve(solver)
    if pulp.LpStatus[prob.status] == 'Optimal':
        return [int(pulp.value(var)) for var in x_vars]
    else:
        raise RuntimeError("Failed to find a valid solution.")

def get_molecular_formula(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return Chem.rdMolDescriptors.CalcMolFormula(molecule)
    else:
        return "Invalid SMILES string"

def balance_chemical_equation(reactant_smiles, product_smiles):
    """ Balance a chemical equation given reactants and products as SMILES strings. """
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
    reactant_coeffs = integer_coefficients[:len(reactant_smiles)]
    product_coeffs = integer_coefficients[len(reactant_smiles):]
    
    reactant_data = [(coeff, get_molecular_formula(smiles)) for coeff, smiles in zip(reactant_coeffs, reactant_smiles)]
    product_data = [(coeff, get_molecular_formula(smiles)) for coeff, smiles in zip(product_coeffs, product_smiles)]

    return reactant_data, product_data
    
def setup_matrix(elements, compounds):
    """ Create a stoichiometry matrix for the elements and compounds. """
    matrix = []
    for element in elements:
        row = [compound.get(element, 0) for compound in compounds]
        matrix.append(row)
    return np.array(matrix, dtype=int)

def display_reaction(reactants, products):
    """Format and display the chemical reaction."""
    def format_component(component):
        try:
            coefficient, molecule = component
            return f"{coefficient} {molecule}" if coefficient != 1 else molecule
        except ValueError:
            raise ValueError(f"Invalid component format: {component}. Expected a tuple of (coefficient, molecule).")

    if not reactants or not products:
        raise ValueError("Both reactants and products need at least one component.")

    try:
        reactants_str = ' + '.join(format_component(r) for r in reactants)
        products_str = ' + '.join(format_component(p) for p in products)
        return f"{reactants_str} → {products_str}"
    except ValueError as e:
        print(e)
        return None  # or handle differently

def create_reaction_string(reactants, products):
    reactants_str = '.'.join(reactants)
    products_str = '.'.join(products)
    return f"{reactants_str}>>{products_str}"

def display_svg(svg):
    """Display SVG in Streamlit using markdown with unsafe HTML."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = f"<img src='data:image/svg+xml;base64,{b64}'/>"
    st.markdown(html, unsafe_allow_html=True)
