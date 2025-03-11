from vrs_anvil.phenotype import get_child_terms


def test_get_child_terms():
    result = get_child_terms("HP:0030850")
    assert result == {"HP:0030850", "HP:0030852", "HP:0030851"}
