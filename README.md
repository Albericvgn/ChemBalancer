![project logo](assets/prog.png)

# Chemical Equation Balancer

The Chemical Equation Equilibrator is a Python program that balances and equilibrates chemical equations based on reactants and products provided as SMILES strings.


## Features

- SMILES or name Input: Accepts reactants and products as SMILES or chemical compound name strings. This provides a convenient way to input reactants and products.
- Chemical Equation Balancing: Automatically balances chemical equations by solving an integer linear programming problem to find stoichiometric coefficients.
- Display Reaction: Formats and displays the balanced chemical reaction in a human-readable format.
- For a specific temperature, it calculates standard enthropy, enthalpy & free gibbs energy in order to display them and tell wether the reaction is spotaneous or not.

## Requirements

- [NumPy](https://github.com/numpy/numpy): Required for handling arrays and matrices.
- [RDKit](https://github.com/rdkit/rdkit): A cheminformatics toolkit used for parsing SMILES strings and generating 2D molecular structures.
- [PuLP](https://coin-or.github.io/pulp/): A linear programming library used for solving integer linear programming problems.
- [rxnmapper](https://github.com/rxn4chemistry/rxnmapper): A package for mapping reactions.
- [IPython](https://ipython.org/): Interactive computing in Python.
- [requests](https://requests.readthedocs.io/en/latest/): HTTP library for making requests.
- [base64](https://docs.python.org/3/library/base64.html#): Library for base64 encoding and decoding.
- [chemicals](https://github.com/CalebBell/chemicals): A package for chemical properties.
- [streamlit](https://streamlit.io/): Used for creating interactive web apps.
- [Rxn-INSIGHT](https://github.com/schwallergroup/Rxn-INSIGHT): Open-source algorithm, written in python, to classify and name chemical reactions, and suggest reaction conditions based on similarity and popularity.

All of the test environments can be run using the command tox from the top directory.
Alternatively, individual test environments can be run using the -e flag as 
in tox -e env-name. To run the tests, tests with coverage report, style checks, and
docs build, respectively:
```
tox -e py3
tox -e py3-coverage
tox -e style
tox -e docs
```
## Installation 

Ensure you have the required dependencies installed: NumPy, RDKit, and PuLP, PIL, rxnmapper, IPython, request and Rxn-INSIGHT.
```
pip install numpy
```
```
pip install -c conda-forge rdkit
```
```
pip install pulp
```
```
pip install rxnmapper
```
```
pip install ipython
```
```
pip install requests 
```
```
pip install pybase64 
```
```
pip install chemicals
```
```
pip install streamlit
```
A virtual environment can be installed with Anaconda as follows:
```
console
conda create -n rxn-insight python=3.10
conda activate rxn-insight
```

```
git clone https://github.com/schwallergroup/Rxn-INSIGHT.git
cd Rxn-INSIGHT
pip install .
```

Or, for developing with the optional dependencies, which are required to run the tests
and build the docs:
```
pip install -e ".[test,doc]"
```


## Usage

1. Clone this repository and navigate to it:
  ```
git clone https://github.com/Albericvgn/ChemBalancer
cd ChemBalancer
  ```

2. To open the application, run the following command:
  ```
streamlit run ChemBalancerInterface.py
  ```

3. Once the streamlit application is open and running, follow these steps to balance your reaction:
   - Select the desired format that you would like to input the reactants/products in (name or SMILES).
   - Input these compounds, separating them with a comma.
   - Input the temperature at which you wish to study this reaction (in Kelvins).
   - Click on the 'Balance Equation' button to generate the balanced reaction as well as the thermochemical data.


##  Reference

`M. R. Dobbelaere, I. Lengyel, C. V. Stevens, and K. M. Van Geem, 
‘Rxn-INSIGHT: fast chemical reaction analysis using bond-electron matrices’, J. Cheminform., vol. 16, no. 1, Mar. 2024.`

```
@ARTICLE{Dobbelaere2024-es,
  title     = "{Rxn-INSIGHT}: fast chemical reaction analysis using
               bond-electron matrices",
  author    = "Dobbelaere, Maarten R and Lengyel, Istv{\'a}n and Stevens,
               Christian V and Van Geem, Kevin M",
  journal   = "J. Cheminform.",
  publisher = "Springer Science and Business Media LLC",
  volume    =  16,
  number    =  1,
  month     =  mar,
  year      =  2024,
  copyright = "https://creativecommons.org/licenses/by/4.0",
  language  = "en"
}
```
