# test_get_smiles_from_name.py

import pytest
import requests
from requests.exceptions import HTTPError
from requests_mock import Mocker

from chembalancer import get_smiles_from_name  # Replace 'your_module' with the actual module name

def test_get_smiles_from_name_success():
    mock_response = {
        "PropertyTable": {
            "Properties": [{"CanonicalSMILES": "C[C@@H](C(=O)O)N"}]
        }
    }

    with Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetaminophen/property/CanonicalSMILES/JSON',
            json=mock_response,
            status_code=200
        )
        smiles = get_smiles_from_name('acetaminophen')
        assert smiles == "C[C@@H](C(=O)O)N"

def test_get_smiles_from_name_error():
    with Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/nonexistent/property/CanonicalSMILES/JSON',
            status_code=404
        )
        smiles = get_smiles_from_name('nonexistent')
        assert smiles == "No data found or error occurred."

def test_get_smiles_from_name_exception():
    with Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/error/property/CanonicalSMILES/JSON',
            exc=HTTPError("HTTP Error")
        )
        smiles = get_smiles_from_name('error')
        assert smiles == "No data found or error occurred."

