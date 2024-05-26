![project logo](assets/input.jpg)

# Chemical Equation Balancer

The Chemical Equation Equilibrator is a Python program that balances and equilibrates chemical equations based on reactants and products provided as SMILES strings.


## Features

- SMILES or name Input: Accepts reactants and products as SMILES or chemical compound name strings. This provides a convenient way to input reactants and products.
- Chemical Equation Balancing: Automatically balances chemical equations by solving an integer linear programming problem to find stoichiometric coefficients.
- Display Reaction: Formats and displays the balanced chemical reaction in a human-readable format.
- For a specific temperature, it calculates standard enthropy, enthalpy & free gibbs energy in order to display them and tell wether the reaction is spotaneous or not.

## Requirements

- NumPy: Required for handling arrays and matrices.

```
pip install numpy
```

- RDKit: A cheminformatics toolkit used for parsing SMILES strings and generating 2D molecular structures.

```
pip install -c conda-forge rdkit
```

- PuLP: A linear programming library used for solving integer linear programming problems.

```
pip install pulp
```

- PIL: Python Imaging Library used for image manipulation.

```
pip install pillow
```

- rxnmapper: A package for mapping reactions.

```
pip install rxnmapper
```

- IPython:Interactive computing in Python.

```
pip install ipython
```

- request:HTTP library for making requests.

```
pip install request 
```

- base64 : Library for base64 encoding and decoding.

```
pip install pybase64 
```

- chemicals: A package for chemical properties.

```
pip install chemicals
```

- streamlit: For creating interactive web apps

```
pip install streamlit
```


- Rxn-INSIGHT relies on NumPy, Pandas, RDKit, RDChiral, and RXNMapper.

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
## Installing ChemBalencer 
you can use pip:

```
pip install ChemBalancer
```

## Usage

1. Ensure you have the required dependencies installed: NumPy, RDKit, and PuLP, PIL, rxnmapper, IPython, request and Rxn-INSIGHT.
2. Run the Python script chemical_equation_equilibrator.py.
3. Follow the prompts to input the number of reactants and products, as well as the SMILES or names strings for each.
4. The program will automatically balance the chemical equation.
5. It will display the balanced reaction.
6. It will calculate standard enthropy, enthalpy and free Gibbs energy at a given temperature and display them.


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
