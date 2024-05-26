
import pytest
import requests
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

from chembalancer import test_hello_smiles  


# Test the function
def test_hello_smiles():
    assert hello_smiles("C(=O)O") == "Hello, C(=O)O!", "Test failed: SMILES input"
