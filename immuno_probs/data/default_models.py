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


def get_default_model_file_paths(model_name):
    """Returns a directory with file paths for a given model name.

    model_name : str
        A string value representing a model name in the dictionary.

    Returns
    -------
    dict
        Containing model marginals, model parameters, v_anchors and j_anchors
        file paths. If model name does not exist in the dictionary, return None.

    """
    # Set the file paths for the models
    file_path = os.path.dirname(__file__)
    default_models = {
        'human-t-alpha': {
            'marginals': os.path.join(file_path, 'human_T_alpha', 'model_marginals.txt'),
            'parameters': os.path.join(file_path, 'human_T_alpha', 'model_params.txt'),
            'v_anchors': os.path.join(file_path, 'human_T_alpha', 'V_gene_CDR3_anchors.csv'),
            'j_anchors': os.path.join(file_path, 'human_T_alpha', 'J_gene_CDR3_anchors.csv'),
        },
        'human-t-beta': {
            'marginals': os.path.join(file_path, 'human_T_beta', 'model_marginals.txt'),
            'parameters': os.path.join(file_path, 'human_T_beta', 'model_params.txt'),
            'v_anchors': os.path.join(file_path, 'human_T_beta', 'V_gene_CDR3_anchors.csv'),
            'j_anchors': os.path.join(file_path, 'human_T_beta', 'J_gene_CDR3_anchors.csv'),
        },
        'human-b-heavy': {
            'marginals': os.path.join(file_path, 'human_B_heavy', 'model_marginals.txt'),
            'parameters': os.path.join(file_path, 'human_B_heavy', 'model_params.txt'),
            'v_anchors': os.path.join(file_path, 'human_B_heavy', 'V_gene_CDR3_anchors.csv'),
            'j_anchors': os.path.join(file_path, 'human_B_heavy', 'J_gene_CDR3_anchors.csv'),
        },
        'mouse-t-beta': {
            'marginals': os.path.join(file_path, 'mouse_T_beta', 'model_marginals.txt'),
            'parameters': os.path.join(file_path, 'mouse_T_beta', 'model_params.txt'),
            'v_anchors': os.path.join(file_path, 'mouse_T_beta', 'V_gene_CDR3_anchors.csv'),
            'j_anchors': os.path.join(file_path, 'mouse_T_beta', 'J_gene_CDR3_anchors.csv'),
        },
    }

    # For the species and chain return dict with file paths
    for name, file_paths in default_models.items():
        if name == model_name:
            return file_paths
    return None
