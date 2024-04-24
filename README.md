COYS: Compute protein properties and Optimize Your Sequences
===============================================================
# Usage

## Functions

### calculate_properties(seq)

Calculates various physicochemical properties of a protein sequence. The properties include length, molecular_weight, charge, hydrophobic_ratio, isoelectric_point, hydrophobic_moment.

#### Input:
- `seq`: str, the protein sequence

#### Returns:
- `properties_dict`: dict, a dictionary containing the calculated properties

### sequence_summary(seq)

Summarizes the properties of a protein sequence and visualizes the amino acid distribution.

#### Input:
- `seq`: str, the protein sequence

### generate_seq(generations=5, aa_list=None, length=5, **kwargs)

Generates a protein sequence optimized for user-defined properties using a genetic algorithm.

#### Inputs:
- `generations`: int, the number of generations for the genetic algorithm (default: 5)
- `aa_list`: list or None, a list of amino acids to consider during sequence generation (default: None, uses all 20 standard amino acids)
- `length`: int, the length of the generated sequence (default: 5)
- `**kwargs`: dict, keyword arguments specifying the target properties and their desired values. The key words arguments include:
  - `molecular_weight`
  - `charge`
  - `hydrophobic_ratio`
  - `isoelectric_point`
  - `hydrophobic_moment`

#### Returns:
- `best_seq`: str, the generated peptide sequence
- `best_fitness`: float, the fitness score of the generated sequence
- `seq_properties`: dict, a dictionary containing the properties of the generated sequence
- `ga_aaconst`: object, the genetic algorithm object used for sequence generation
