import pytest
from chemicals import CAS_from_any
from chembalancer.chembalancer import compound_state

@pytest.mark.parametrize(
    "compound, temp_kelvin, expected_state",
    [
        ("water", 298.15, "liquid"),    
        ("water", 383.15, "gas"),       
        ("ethanol", 298.15, "liquid"),  
        ("ethanol", 351.15, "liquid"),     
        ("iron", 298.15, "solid"),      
        ("oxygen", 73.15, "liquid"),       
    ]
)
def test_compound_state(compound, temp_kelvin, expected_state):
    assert compound_state(compound, temp_kelvin) == expected_state
