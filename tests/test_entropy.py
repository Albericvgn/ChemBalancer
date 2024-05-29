import pytest
from unittest.mock import patch
from chemicals import CAS_from_any, S0s, S0l, S0g
from chembalancer.chembalancer import entropy

@pytest.mark.parametrize(
    "coeff, compound, state, expected_entropy",
    [
        (1.0, "water", "solid", 0.0),      
        (2.0, "water", "liquid", 141),  
        (3.0, "water", "gas", 566.4),     
        (1.0, "ethanol", "solid", 0.0),  
        (2.0, "ethanol", "liquid", 321),   
        (3.0, "ethanol", "gas", 845),     
        (1.0, "iron", "solid", 27.2),      
        (2.0, "iron", "liquid", 69),      
        (3.0, "iron", "gas", 540),  
    ]
)

def test_entropy(coeff, compound, state, expected_entropy):
    with patch("chemicals.CAS_from_any") as mock_CAS_from_any:
        mock_CAS_from_any.return_value = "dummy_cas"
        
        with patch("chemicals.S0s") as mock_S0s, \
             patch("chemicals.S0l") as mock_S0l, \
             patch("chemicals.S0g") as mock_S0g:
            
            if state == "solid":
                mock_S0s.return_value = expected_entropy
                mock_S0l.return_value = None 
                mock_S0g.return_value = None
            elif state == "liquid":
                mock_S0s.return_value = None
                mock_S0l.return_value = expected_entropy
                mock_S0g.return_value = None
            else:
                mock_S0s.return_value = None 
                mock_S0l.return_value = None 
                mock_S0g.return_value = expected_entropy  
            
            result = entropy(coeff, compound, state)
            
            assert result == pytest.approx(expected_entropy, rel=1e-2)
