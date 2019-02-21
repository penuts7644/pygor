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
from pkg_resources import resource_filename


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
    pkg_name = __name__.split('.')[0]
    default_models = {
        'human-t-alpha': {
            'marginals': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'model_marginals.txt')),
            'parameters': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'model_params.txt')),
            'v_anchors': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'V_gene_CDR3_anchors.csv')),
            'j_anchors': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'J_gene_CDR3_anchors.csv')),
            'reference': {
                'V': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'genomic_V.fasta')),
                'J': resource_filename(pkg_name, os.path.join('data', 'human_t_alpha', 'genomic_J.fasta')),
            },
        },
        'human-t-beta': {
            'marginals': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'model_marginals.txt')),
            'parameters': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'model_params.txt')),
            'v_anchors': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'V_gene_CDR3_anchors.csv')),
            'j_anchors': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'J_gene_CDR3_anchors.csv')),
            'reference': {
                'V': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'genomic_V.fasta')),
                'D': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'genomic_D.fasta')),
                'J': resource_filename(pkg_name, os.path.join('data', 'human_t_beta', 'genomic_J.fasta')),
            },
        },
        'human-b-heavy': {
            'marginals': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'model_marginals.txt')),
            'parameters': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'model_params.txt')),
            'v_anchors': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'V_gene_CDR3_anchors.csv')),
            'j_anchors': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'J_gene_CDR3_anchors.csv')),
            'reference': {
                'V': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'genomic_V.fasta')),
                'D': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'genomic_D.fasta')),
                'J': resource_filename(pkg_name, os.path.join('data', 'human_b_heavy', 'genomic_J.fasta')),
            },
        },
        'mouse-t-beta': {
            'marginals': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'model_marginals.txt')),
            'parameters': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'model_params.txt')),
            'v_anchors': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'V_gene_CDR3_anchors.csv')),
            'j_anchors': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'J_gene_CDR3_anchors.csv')),
            'reference': {
                'V': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'genomic_V.fasta')),
                'D': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'genomic_D.fasta')),
                'J': resource_filename(pkg_name, os.path.join('data', 'mouse_t_beta', 'genomic_J.fasta')),
            },
        },
    }
    # For the species and chain return dict with file paths
    for name, file_paths in default_models.items():
        if name == model_name:
            return file_paths
    return None
