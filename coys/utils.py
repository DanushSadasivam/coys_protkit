import numpy as np
import pandas as pd
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor
from modlamp.plot import plot_aa_distr
import pygad


def calculate_properties(seq):
    desc = GlobalDescriptor(seq)
    desc.length()
    desc.calculate_MW(append=True)
    desc.calculate_charge(append=True)
    desc.hydrophobic_ratio(append=True)
    desc.isoelectric_point(append=True)
    properties_vals = desc.descriptor.flatten()
    desc = PeptideDescriptor(seq)
    desc.calculate_moment()
    properties_vals = np.append(properties_vals, desc.descriptor.flatten()[0])
    properties_names = [
        "length",
        "molecular_weight",
        "charge",
        "hydrophobic_ratio",
        "isoelectric_point",
        "hydrophobic_moment",
    ]
    properties_dict = dict(zip(properties_names, properties_vals))
    return properties_dict


def sequence_summary(seq):
    properties_dict = calculate_properties(seq)
    properties_df = pd.DataFrame()
    properties_df["Property"] = list(properties_dict.keys())
    properties_df["Value"] = list(properties_dict.values())
    print(properties_df)
    plot_aa_distr(seq)


def generate_seq(generations=5, aa_list=None, length=5, **kwargs):
    target_dict = kwargs.copy()

    aa_keys = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    int_labels = range(1, 21)
    num2aa_table = dict(zip(int_labels, aa_keys))
    aa2num_table = dict(zip(aa_keys, int_labels))

    # AASeq to Num array func
    def aaseq2num(seq):
        arr = []
        for char in seq:
            arr.append(aa2num_table[char])

        return arr

    # Num array to AA seq func
    def num2aaseq(seq):
        aa_seq = ""
        for j, char in enumerate(seq):
            aa_seq = aa_seq + num2aa_table[char]
        return aa_seq

    # Fitness Function
    def prop_fitness_func(ga_aaconst, x, x_idx):
        seq = num2aaseq(x)
        properties_dict = calculate_properties(seq)
        total_dist = 0
        for prop in list(target_dict.keys()):
            prop_dist = (target_dict[prop] - properties_dict[prop]) ** 2
            total_dist += prop_dist

        total_dist = -1 * total_dist
        return total_dist

    l = length
    pop_size = l
    gene_arr = []
    if aa_list is None:
        aa_list = aa_keys.copy()

    for char in aa_list:
        gene_arr.append(aa2num_table[char])
    gene_space = gene_arr

    ga_aaconst = pygad.GA(
        num_generations=generations,
        num_parents_mating=2,
        sol_per_pop=pop_size,
        num_genes=l,
        fitness_func=prop_fitness_func,
        parent_selection_type="rank",
        keep_elitism=5,
        suppress_warnings=True,
        stop_criteria=["reach_0"],
        gene_type=int,
        gene_space=gene_space,
    )
    ga_aaconst.run()
    best_seq, best_fitness, _ = ga_aaconst.best_solution()
    best_seq = num2aaseq(best_seq)
    seq_properties = calculate_properties(best_seq)
    print("The best sequence is:", best_seq)
    return best_seq, best_fitness, seq_properties, ga_aaconst
