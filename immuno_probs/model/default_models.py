# Create IGoR models and calculate the generation probability of V(D)J and
# CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt

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


"""Contains function for loading in and using pre-trained V(D)J models."""


import os
from pkg_resources import resource_filename


def get_default_model_file_paths(name=None):
    """Returns a directory with file paths for a given model identifier name.

    Parameters
    ----------
    name : str, optional
        A string value representing a model identifier name in the dictionary. If name is not specified, returns a list
        containing all the available model options.

    Returns
    -------
    dict
        Containing model marginals, model parameters, anchors and reference genome file paths. If model name does not exist
        in the dictionary, returns None instead.

    """
    # Set the file paths for the models
    pkg_name = __name__.split('.')[0]
    default_models = {
        'tutorial-model': {
            'type': 'beta',
            'marginals': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'model_marginals.txt')),
            'parameters': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'model_params.txt')),
            'v_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'V_gene_CDR3_anchors.tsv')),
            'j_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'J_gene_CDR3_anchors.tsv')),
            'reference': {
                'V': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBV.fasta')),
                'D': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBD.fasta')),
                'J': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBJ.fasta')),
            },
            'seqs': resource_filename(
                pkg_name, os.path.join('data', 'tutorial_model', 'VDJ_sequences.tsv')),
            'cdr3': resource_filename(
                pkg_name, os.path.join('data', 'tutorial_model', 'CDR3_sequences.tsv')),
        },
        'human-t-alpha': {
            'type': 'alpha',
            'marginals': resource_filename(
                pkg_name, os.path.join('data', 'human_t_alpha', 'model_marginals.txt')),
            'parameters': resource_filename(
                pkg_name, os.path.join('data', 'human_t_alpha', 'model_params.txt')),
            'v_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_alpha', 'V_gene_CDR3_anchors.tsv')),
            'j_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_alpha', 'J_gene_CDR3_anchors.tsv')),
            'reference': {
                'V': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_alpha', 'TRAV.fasta')),
                'J': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_alpha', 'TRAJ.fasta')),
            },
        },
        'human-t-beta': {
            'type': 'beta',
            'marginals': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'model_marginals.txt')),
            'parameters': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'model_params.txt')),
            'v_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'V_gene_CDR3_anchors.tsv')),
            'j_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_t_beta', 'J_gene_CDR3_anchors.tsv')),
            'reference': {
                'V': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBV.fasta')),
                'D': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBD.fasta')),
                'J': resource_filename(
                    pkg_name, os.path.join('data', 'human_t_beta', 'TRBJ.fasta')),
            },
        },
        'human-b-heavy': {
            'type': 'heavy',
            'marginals': resource_filename(
                pkg_name, os.path.join('data', 'human_b_heavy', 'model_marginals.txt')),
            'parameters': resource_filename(
                pkg_name, os.path.join('data', 'human_b_heavy', 'model_params.txt')),
            'v_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_b_heavy', 'V_gene_CDR3_anchors.tsv')),
            'j_anchors': resource_filename(
                pkg_name, os.path.join('data', 'human_b_heavy', 'J_gene_CDR3_anchors.tsv')),
            'reference': {
                'V': resource_filename(
                    pkg_name, os.path.join('data', 'human_b_heavy', 'IGHV.fasta')),
                'D': resource_filename(
                    pkg_name, os.path.join('data', 'human_b_heavy', 'IGHD.fasta')),
                'J': resource_filename(
                    pkg_name, os.path.join('data', 'human_b_heavy', 'IGHJ.fasta')),
            },
        },
        'mouse-t-beta': {
            'type': 'beta',
            'marginals': resource_filename(
                pkg_name, os.path.join('data', 'mouse_t_beta', 'model_marginals.txt')),
            'parameters': resource_filename(
                pkg_name, os.path.join('data', 'mouse_t_beta', 'model_params.txt')),
            'v_anchors': resource_filename(
                pkg_name, os.path.join('data', 'mouse_t_beta', 'V_gene_CDR3_anchors.tsv')),
            'j_anchors': resource_filename(
                pkg_name, os.path.join('data', 'mouse_t_beta', 'J_gene_CDR3_anchors.tsv')),
            'reference': {
                'V': resource_filename(
                    pkg_name, os.path.join('data', 'mouse_t_beta', 'TRBV.fasta')),
                'D': resource_filename(
                    pkg_name, os.path.join('data', 'mouse_t_beta', 'TRBD.fasta')),
                'J': resource_filename(
                    pkg_name, os.path.join('data', 'mouse_t_beta', 'TRBJ.fasta')),
            },
        },
    }

    # If name given, for the species and chain return dict with file paths
    if name:
        for model, file_paths in default_models.items():
            if model == name:
                return file_paths
    else:
        return list(default_models)
