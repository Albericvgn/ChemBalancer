# test_get_smiles_from_name.py

import pytest
import requests
import requests_mock
from requests.exceptions import HTTPError
from requests_mock import Mocker
import sys
import os

try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
src_path = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_path)

from chembalancer.chembalancer import get_smiles_from_name  

# Test for successful response
def test_get_smiles_from_name_success():
    mock_response = {
        "PropertyTable": {
            "Properties": [{"CanonicalSMILES": "C[C@@H](C(=O)O)N"}]
        }
    }

    with requests_mock.Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetaminophen/property/CanonicalSMILES/JSON',
            json=mock_response,
            status_code=200
        )
        smiles = get_smiles_from_name('acetaminophen')
        assert smiles == "C[C@@H](C(=O)O)N"

# Test for 404 response
def test_get_smiles_from_name_error():
    with requests_mock.Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/nonexistent/property/CanonicalSMILES/JSON',
            status_code=404
        )
        smiles = get_smiles_from_name('nonexistent')
        assert "HTTP Error:" in smiles  # Adjusted assertion

# Test for HTTP error response
def test_get_smiles_from_name_exception():
    with requests_mock.Mocker() as m:
        m.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/error/property/CanonicalSMILES/JSON',
            exc=HTTPError("HTTP Error")
        )
        smiles = get_smiles_from_name('error')
        assert "HTTP Error:" in smiles  # Adjusted assertion

if __name__ == "__main__":
    pytest.main()
