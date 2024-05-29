import streamlit as st
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from collections import defaultdict
import pulp
from rxnmapper import RXNMapper
from rxn_insight.reaction import Reaction
from rxn_insight.utils import draw_chemical_reaction, curate_smirks, get_similarity, get_fp
from IPython.display import SVG
import time
import requests
import base64
from io import BytesIO
from chemicals import CAS_from_any, Tb, Tm, Tc, Hfs, Hfl, Hfg, S0s, S0l, S0g
from collections import defaultdict
from rdkit import Chem
import requests
import base64


def get_smiles_from_name(name):
    """Fetch the SMILES string of a molecule by its common name from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise HTTPError for bad responses
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        return smiles
    except requests.exceptions.HTTPError as errh:
        # Handle HTTPError
        return f"HTTP Error: {errh}"
    except requests.exceptions.RequestException as err:
        # Handle other RequestExceptions
        return f"Request Exception: {err}"
    except (KeyError, IndexError):
        # Handle missing data or incorrect JSON format
        return "No data found or error occurred."



def count_atoms(smiles):
    """ Count atoms in a SMILES string. """
    mol = Chem.MolFromSmiles(smiles)
    atom_counts = defaultdict(int)
    if mol:
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            atom_counts[atom.GetSymbol()] += 1
    return dict(atom_counts)

def solve_ilp(A):
    """ Solve the integer linear programming problem to find stoichiometric coefficients. """
    num_vars = A.shape[1]
    prob = pulp.LpProblem("Balancing_Chemical_Equation", pulp.LpMinimize)
    
    # Define variables with lower bound starting from 1
    x_vars = [pulp.LpVariable(f'x{i}', lowBound=1, cat='Integer') for i in range(num_vars)]
    
    # Objective function
    prob += pulp.lpSum(x_vars)
    
    # Constraints
    for i in range(A.shape[0]):
        prob += pulp.lpDot(A[i, :], x_vars) == 0
    
    # Solve the problem
    solver = pulp.PULP_CBC_CMD(msg=True)  # Enable logging from the solver
    prob.solve(solver)
    
    print(f"Status: {pulp.LpStatus[prob.status]}")  # Print the status of the problem
    
    if pulp.LpStatus[prob.status] == 'Optimal':
        solution = [int(pulp.value(var)) for var in x_vars]
        print(f"Solution: {solution}")  # Print the solution found
        
        # Check if solution is not just zeros
        if all(x == 0 for x in solution):
            return None
        
        return solution
    else:
        return None



def get_molecular_formula(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is not None:
        return Chem.rdMolDescriptors.CalcMolFormula(molecule)
    else:
        return "Invalid SMILES string"


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
    return html


def compound_state(compound, temp):
    CAS_compound = CAS_from_any(compound)
    boiling_p = Tb(CAS_compound)
    melting_p = Tm(CAS_compound)

    if float(temp) <= float(melting_p):
        return 'solid'
    elif float(temp) >= float(boiling_p):
        return 'gas'
    else:
        return 'liquid'

def enthalpy(coeff, compound, state):
    Cas_compound=CAS_from_any(compound)
    if state == 'solid': 
        if Hfs(Cas_compound)== None:
            return 0
        else: 
            return float(coeff) * Hfs(Cas_compound)
    elif state == 'liquid':
        if Hfl(CAS_from_any(compound))== None:
            return 0
        else:
            return float(coeff) * Hfl(Cas_compound)
    else: 
        if Hfg(CAS_from_any(compound))== None:
            return 0
        else:
            return float(coeff) * Hfg(Cas_compound)
                                    
def entropy(coeff, compound, state):
    Cas_compound=CAS_from_any(compound)
    if state == 'solid': 
        if S0s(Cas_compound)== None:
            return 0
        else: 
            return float(coeff) * S0s(Cas_compound)
    elif state == 'liquid':
        if S0l(CAS_from_any(compound))== None:
            return 0
        else:
            return float(coeff) * S0l(Cas_compound)
    else: 
        if S0g(CAS_from_any(compound))== None:
            return 0
        else:
            return float(coeff) * S0g(Cas_compound)

                                    
