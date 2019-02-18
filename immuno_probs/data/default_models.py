# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""Contains function for loading and using build-in V(D)J models."""


import os


def is_valid_model(species, chain):
    """Checks if a given species and chain combination exists.

    species : str
        A string value representing a species name.
    chain : str
        A string value representing a chain name.

    Returns
    -------
    bool
        Specifying if the species and chain combination are supported.

    """
    species_chain_combination = {
        'human': ['alpha', 'beta', 'heavy'],
        'mouse': ['beta'],
    }
    for key, value in species_chain_combination.items():
        if key == species:
            if chain in value:
                return True
    return False


def get_default_model_file_paths(species, chain):
    """Returns a directory with file paths for a given species and chain.

    species : str
        A string value representing a species name.
    chain : str
        A string value representing a chain name.

    Returns
    -------
    dict
        Containing model marginals, model parameters, v_anchors and j_anchors
        file paths. If species and/or chain does not exist, return None.

    """
    # Set the file paths for the models
    file_path = os.path.dirname(__file__)
    default_models = {
        'human': {
            'alpha': {
                'marginals': os.path.join(file_path, 'human_T_alpha', 'model_marginals.txt'),
                'parameters': os.path.join(file_path, 'human_T_alpha', 'model_params.txt'),
                'v_anchors': os.path.join(file_path, 'human_T_alpha', 'V_gene_CDR3_anchors.csv'),
                'j_anchors': os.path.join(file_path, 'human_T_alpha', 'J_gene_CDR3_anchors.csv'),
            },
            'beta': {
                'marginals': os.path.join(file_path, 'human_T_beta', 'model_marginals.txt'),
                'parameters': os.path.join(file_path, 'human_T_beta', 'model_params.txt'),
                'v_anchors': os.path.join(file_path, 'human_T_beta', 'V_gene_CDR3_anchors.csv'),
                'j_anchors': os.path.join(file_path, 'human_T_beta', 'J_gene_CDR3_anchors.csv'),
            },
            'heavy': {
                'marginals': os.path.join(file_path, 'human_B_heavy', 'model_marginals.txt'),
                'parameters': os.path.join(file_path, 'human_B_heavy', 'model_params.txt'),
                'v_anchors': os.path.join(file_path, 'human_B_heavy', 'V_gene_CDR3_anchors.csv'),
                'j_anchors': os.path.join(file_path, 'human_B_heavy', 'J_gene_CDR3_anchors.csv'),
            },
        },
        'mouse': {
            'beta': {
                'marginals': os.path.join(file_path, 'mouse_T_beta', 'model_marginals.txt'),
                'parameters': os.path.join(file_path, 'mouse_T_beta', 'model_params.txt'),
                'v_anchors': os.path.join(file_path, 'mouse_T_beta', 'V_gene_CDR3_anchors.csv'),
                'j_anchors': os.path.join(file_path, 'mouse_T_beta', 'J_gene_CDR3_anchors.csv'),
            },
        },
    }

    # For the species and chain return dict with file paths
    for key1, value1 in default_models.items():
        if key1 == species:
            for key2, value2 in value1.items():
                if key2 == chain:
                    return value2
    return None
