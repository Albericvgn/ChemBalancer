![project logo](assets/banner_project_prog.jpg)

# Chemical Equation Balancer

The Chemical Equation Balancer is a Python program that balances chemical reaction equations based on reactants and products provided by the user in the desired format (name or SMILES).

## Features

- SMILES or name Input: Accepts reactants and products as SMILES or chemical compound name strings. This provides a convenient way to input reactants and products.
- Chemical Equation Balancing: Automatically balances chemical equations by solving an integer linear programming problem to find stoichiometric coefficients.
- Display Reaction: Formats and displays the balanced chemical reaction in a human-readable format.
- For a specific temperature, it calculates standard enthropy, enthalpy & free gibbs energy in order to display them and tell wether the reaction is spotaneous or not.

## Requirements

- [NumPy](https://github.com/numpy/numpy): Required for handling arrays and matrices.
- [RDKit](https://github.com/rdkit/rdkit): A cheminformatics toolkit used for parsing SMILES strings and generating 2D molecular structures.
- [PuLP](https://coin-or.github.io/pulp/): A linear programming library used for solving integer linear programming problems.
- [RXNMapper](https://github.com/rxn4chemistry/rxnmapper): A package for mapping reactions.
- [IPython](https://ipython.org/): Interactive computing in Python.
- [requests](https://requests.readthedocs.io/en/latest/): HTTP library for making requests.
- [base64](https://docs.python.org/3/library/base64.html#): Library for base64 encoding and decoding.
- [chemicals](https://github.com/CalebBell/chemicals): A package for chemical properties.
- [streamlit](https://streamlit.io/): Used for creating interactive web apps.
- [Rxn-INSIGHT](https://github.com/schwallergroup/Rxn-INSIGHT): Open-source algorithm, written in python, to classify and name chemical reactions, and suggest reaction conditions based on similarity and popularity.

## Installation 

Before installing any of the required dependencies, it is recommended that you create and activate an environment to do this in as to avoid any potential conflicts.
1. Clone this repository and navigate to it:
  ```
git clone https://github.com/Albericvgn/ChemBalancer
cd ChemBalancer
  ```
2. Create a conda environment.
```
conda create -n your_env_name python=3.10
```
3. Activate this environment.
```
conda activate your_env_name
```
Now that the environment is activated, the dependencies required for the code must be installed.
1. Installing NumPy:
```
pip install numpy
```
2. Installing RDKit:
```
conda install -c conda-forge rdkit
```
3. Installing PuLP:
```
pip install pulp
```
4. Installing RXNMapper:
```
pip install rxnmapper
```
5. Installing IPython:
```
pip install ipython
```
6. Installing requests:
```
pip install requests 
```
7. Installing base64:
```
pip install pybase64 
```
8. Installing chemicals:
```
pip install chemicals
```
9. Installing streamlit:
```
pip install streamlit
```
10. Installing Rxn-INSIGHT:
It should be noted that Rxn-INSIGHT requires the following: NumPy, Pandas, RDKit, RDChiral, and RXNMapper. While some of these are required for our chemical equation balancer and should already be installed, Pandas and RDChiral may need to also be installed:
```
pip install pandas
pip install rdchiral
```
Once all the requirements are met, the following command can be used for the installation:
```
git clone https://github.com/schwallergroup/Rxn-INSIGHT.git
cd Rxn-INSIGHT
pip install .
```

For more details on the installation of this package, you can click on the following link: https://github.com/schwallergroup/Rxn-INSIGHT

Once all the dependencies are installed, you can run your Streamlit application. Make sure you are in the directory containing your script and then run:

```
streamlit run your_script_name.py
```

## Usage

1. Activate the environment and navigate to ChemBalancer:
  ```
conda activate your_env_name
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
   - Use the 'Clear' button if you wish to balance a new equation with different reactants/products.


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

## Developpers

- Fane SHALA, student in chemistry at EPFL. [<img src="./assets/GitHubSymb.png" alt="Profile Picture" width="60">](https://github.com/faneshala)
- John STEWART, student in chemical engineering at EPFL. [<img src="./assets/GitHubSymb.png" alt="Profile Picture" width="60">](https://github.com/johnstewartepfl)
- Albéric VIGNE, student in chemistry at EPFL. [<img src="./assets/GitHubSymb.png" alt="Profile Picture" width="60">](https://github.com/albericvgn)
