import pytest
from coys import *


def test_calculate_properties():
    number_of_properties = 6
    properties_dict = calculate_properties("ACDEFGHKLCMNPQRSTVWY")
    assert len(properties_dict) == number_of_properties


def test_generate_sequence():
    l = 10
    aa_list = ["A", "D", "C", "G", "K"]
    best_seq, best_fitness, seq_properties, ga_obj = generate_seq(
        generations=50,
        length=l,
        aa_list=aa_list,
        charge=1,
        hydrophobic_ratio=0.6,
        hydrophobic_moment=0.6,
    )
    assert len(best_seq) == l
    assert set(best_seq) <= set(aa_list)
