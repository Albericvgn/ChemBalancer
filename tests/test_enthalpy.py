import pytest
from unittest.mock import patch
from chembalancer.chembalancer import enthalpy
from chemicals import CAS_from_any, Tb, Tm, Tc, Hfs, Hfl, Hfg, S0s, S0l, S0g

@pytest.mark.parametrize(
    "coeff, compound, state, expected_enthalpy",
    [
        (1.0, "water", "solid", 0.0),      
        (2.0, "water", "liquid", -571650.0),  
        (3.0, "water", "gas", -725466.0),    
        (1.0, "ethanol", "solid", 0.0),  
        (2.0, "ethanol", "liquid", -554060),
        (3.0, "ethanol", "gas", -703710.0),  
        (1.0, "iron", "solid", 0.0),     
        (2.0, "iron", "liquid", 24800.0),     
        (3.0, "iron", "gas", 1248900.0),          
    ]
)
def test_enthalpy(coeff, compound, state, expected_enthalpy):
    with patch("chemicals.CAS_from_any") as mock_CAS_from_any:
        mock_CAS_from_any.return_value = "dummy_cas"
        
        with patch("chemicals.Hfs") as mock_Hfs, \
             patch("chemicals.Hfl") as mock_Hfl, \
             patch("chemicals.Hfg") as mock_Hfg:
            
            mock_Hfs.return_value = 0.0  
            mock_Hfl.return_value = 333.55  
            mock_Hfg.return_value = 891.8  
            
            result = enthalpy(coeff, compound, state)
            
            assert result == pytest.approx(expected_enthalpy, rel=1e-2)
