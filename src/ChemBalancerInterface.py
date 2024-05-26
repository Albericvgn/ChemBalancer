#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
from chembalancer import *

'''
Creates default session states for potential user inputs which allows the application to remember the values in the event that user wishes to rerun the script.
'''

if 'reactant_input' not in st.session_state:
    st.session_state['reactant_input'] = ""
if 'product_input' not in st.session_state:
    st.session_state['product_input'] = ""
if 'temp' not in st.session_state:
    st.session_state['temp'] = ""
if 'submitted' not in st.session_state:
    st.session_state['submitted'] = False

def prepare_reset():
    '''
    Prepares streamlit app for a reset by the user by setting the session state of 'reset' to True.
    Acts as a trigger to inform the application that a reset has been requested.
    Is called when the user clicks on 'Clear' button.
    '''
    st.session_state['reset'] = True

def reset_form():
    '''
    Checks if the session state of 'reset' has been set to True and, if so, resets the session states of potential user inputs to their default values.
    Is called automatically at the beginning of the script to ensure that inputs are clear before user interacts with them.
    '''
    if st.session_state.get('reset', False):
        st.session_state['reactant_input'] = ""
        st.session_state['product_input'] = ""
        st.session_state['temp'] = ""
        st.session_state['submitted'] = False
        st.session_state['reset'] = False  # Important: reset the flag!

'''
Streamlit application is launched and titled 'Chemical Equation Balancer'.

Labels and entry fields: 
- A label asks the user to choose between the simple name or SMILES format when inputting their reactants/products.
- Two entry fields are provided with labels indicating their use for the user to input their reactants and products in the desired format, using a comma to separate different compounds.
- A label and entry field indicate for the user to enter the temperature at which the reaction is taking place in Kelvins.

Submit button:
A button titled 'Balance Equation' is provided for the user to launch the equation balancer. 
Once clicked, the balance reaction is displayed in string format with the chemical formulas of the compounds and corresponding coefficients.
The reaction is also displayed with the structures of the molecules but without the stoichiometric coefficients.
The thermochemical data for the reaction is also displayed (enthalpy, entropy and Gibbs free energy).
Using the Gibbs free energy, the spontaneity of the reaction is also given.

Reset button:
A button titled 'Clear' allows to user to reinitialise the application, input new reactants/products and rerun the script.
'''

reset_form()

st.title('Chemical Equation Balancer')

input_type = st.radio("Choose the type of input for reactants and products:",
                      ['Names', 'SMILES'])

if input_type == 'Names':
    reactant_input = st.text_input("Enter reactant names (comma-separated):", value=st.session_state['reactant_input'], key='reactant_input')
    product_input = st.text_input("Enter product names (comma-separated):", value=st.session_state['product_input'], key='product_input')
    reactant_names = [name.strip() for name in reactant_input.split(',')]
    product_names = [name.strip() for name in product_input.split(',')]
    reactant_smiles = [get_smiles_from_name(name) for name in reactant_names]
    product_smiles = [get_smiles_from_name(name) for name in product_names]
else:
    reactant_input = st.text_input("Enter reactant SMILES (comma-separated):", value=st.session_state['reactant_input'], key='reactant_input')
    product_input = st.text_input("Enter product SMILES (comma-separated):", value=st.session_state['product_input'], key='product_input')
    reactant_smiles = [name.strip() for name in reactant_input.split(',')]
    product_smiles = [name.strip() for name in product_input.split(',')]

temp = st.number_input("Enter the temperature at which the reaction is taking place (in Kelvin):")

if st.button('Balance Equation'):
    st.session_state['submitted'] = True
    reactant_data, product_data = balance_chemical_equation(reactant_smiles, product_smiles)
    reaction_display = display_reaction(reactant_data, product_data)
    reaction_string = create_reaction_string(reactant_smiles, product_smiles)
    if reaction_display and st.session_state['submitted']:
        st.write("Balanced Chemical Reaction:", reaction_display)
        display_svg(draw_chemical_reaction(reaction_string))
    
    reactant_data_state = [(coeff, mol_formula, compound_state(mol_formula, temp)) for coeff, mol_formula in reactant_data]
    product_data_state = [(coeff, mol_formula, compound_state(mol_formula, temp)) for coeff, mol_formula in product_data]
    H_r = 0
    H_p = 0
    for coeff, reactant, state in reactant_data_state:
        H_r += enthalpy(coeff, reactant, state)
    for coeff, product, state in product_data_state:
        H_p += enthalpy(coeff, product, state)
    H = H_p - H_r
    st.write("Standard change in enthalpy:", H, "J/mol")
    S_r = 0
    S_p = 0
    for coeff, reactant, state in reactant_data_state:
        S_r += entropy(coeff, reactant, state)
    for coeff, product, state in product_data_state:
        S_p += entropy(coeff, product, state)
    S = S_p - S_r
    st.write("Standard change in entropy:", S, "J/mol*K")
    G = H - temp * S
    st.write("Standard change in Gibbs free energy:", G, "J/mol")
    if G < 0:
        st.write("The reaction is spontaneous at this temperature.")
    elif G >= 0:
        st.write("The reaction is not spontaneous at this temperature.")

if st.button("Clear", on_click=prepare_reset):
    st.experimental_rerun()

