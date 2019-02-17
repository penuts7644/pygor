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


"""IgorLoader class for loading a IGoR model including the CDR3 anchor files."""


import olga.load_model as olga_load_model

from immuno_probs.util.exception import ModelLoaderException


class IgorLoader(object):
    """Class loading in an IGoR model with CDR3 anchors.

    Parameters
    ----------
    model_params : string
        A file path location for the IGoR parameters model file.
    model_marginals : string
        A file path location for the IGoR marginals model file.

    Methods
    -------
    load_anchors(model_params, v_anchors, j_anchors)
        Load in the CDR3 V and J gene anchor files.
    is_vj()
        Return a boolean spcifying if the model is VJ or not.
    is_vdj()
        Return a boolean spcifying if the model is VDJ or not.
    get_genomic_data()
        Return the GenomicData OLGA object.
    get_generative_model()
        Return the GenerativeModel OLGA object.

    """
    def __init__(self, model_params, model_marginals):
        super(IgorLoader, self).__init__()
        self.type = self._check_type(model_marginals)
        self.data = self._load_params(model_params)
        self.model = self._load_model(model_marginals)

    @staticmethod
    def _check_type(model_marginals):
        """Private function to check the marginals file for D gene attributes.

        Parameters
        ----------
        model_marginals : string
            A file path location for the IGoR marginals model file.

        Returns
        -------
        boolean
            Specifying if the model marginals are VDJ or not.

        """
        # Parse the marginals file and search for VDJ classifiers.
        v_choice = False
        d_gene = False
        j_choice = False
        with open(model_marginals, 'r') as infile:
            for line in infile:
                if line == '@v_choice\n':
                    v_choice = True
                elif line == '@d_gene\n':
                    d_gene = True
                elif line == '@j_choice\n':
                    j_choice = True
        return (v_choice, d_gene, j_choice)

    def _load_params(self, model_params):
        """Private function for loading genomic parameter data for the IGoR model.

        Parameters
        ----------
        model_params : string
            A file path location for the IGoR parameters model file.

        Returns
        -------
        GenomicDataVJ or GenomicDataVDJ OLGA object
            The genomic data object class for a VJ or V(D)J model.

        Raises
        ------
        ModelLoaderException
            When the model input data cannot be loaded in as either a VJ or
            V(D)J model.

        """
        # Try to load the genomic data model for VDJ.
        try:
            genomic_data = None
            if self.is_vdj():
                genomic_data = olga_load_model.GenomicDataVDJ()
                genomic_data.genD = olga_load_model.read_igor_D_gene_parameters(model_params)
            elif self.is_vj():
                genomic_data = olga_load_model.GenomicDataVJ()
            else:
                raise ModelLoaderException("Model is not VJ or VDJ compliant")

            # Load the remainder of the data for the VJ model and return.
            genomic_data.genV = olga_load_model.read_igor_V_gene_parameters(model_params)
            genomic_data.genJ = olga_load_model.read_igor_J_gene_parameters(model_params)
            return genomic_data

        # If both VDJ, VJ loaders gave an exception, raise custom exception.
        except Exception as err:
            raise ModelLoaderException(err)

    def _load_model(self, model_marginals):
        """Private function for loading the IGoR model marginals.

        Parameters
        ----------
        model_marginals : string
            A file path location for the IGoR marginals model file.

        Returns
        -------
        GenerativeModelVJ or GenerativeModelVDJ OLGA object
            The IGoR generative model object class for a VJ or V(D)J model.

        Raises
        ------
        ModelLoaderException
            When the model input data cannot be loaded in as either a VJ or
            V(D)J model.

        """
        # Try to create the GenerativeModel object for VDJ or VJ.
        try:
            generative_model = None
            if self.is_vdj():
                generative_model = olga_load_model.GenerativeModelVDJ()
            elif self.is_vj():
                generative_model = olga_load_model.GenerativeModelVJ()
            else:
                raise ModelLoaderException("Model is not VJ or VDJ compliant")

            # Load the generative VDJ or VJ model marginals and return.
            generative_model.load_and_process_igor_model(model_marginals)
            return generative_model

        # If both VDJ, VJ loaders gave an exception, raise custom exception.
        except Exception as err:
            raise ModelLoaderException(err)

    def load_anchors(self, model_params, v_anchors, j_anchors):
        """Function for loading the CDR3 acnhor data for the IGoR model.

        Parameters
        ----------
        model_params : string
            A file path location for the IGoR parameters model file.
        v_anchors : string
            File path to the file containing the CDR3 anchors for the V gene.
        j_anchors : string
            File path to the file containing the CDR3 anchors for the J gene.

        Raises
        ------
        ModelLoaderException
            When the model input data cannot be loaded in as either a VJ or
            V(D)J model.

        Notes
        -----
        This function uses the model that has already been loaded in and updates
        the existing model object.

        """
        # Try to load the anchor files into the data model for VDJ.
        try:
            self.data.anchor_and_curate_genV_and_genJ(v_anchors, j_anchors)
            if self.is_vdj():
                self.data.read_VDJ_palindrome_parameters(model_params)
                self.data.generate_cutD_genomic_CDR3_segs()
            elif self.is_vj():
                self.data.read_igor_VJ_palindrome_parameters(model_params)
            else:
                raise ModelLoaderException("Model is not VJ or VDJ compliant")

            # Load the remainder of the data model
            self.data.generate_cutV_genomic_CDR3_segs()
            self.data.generate_cutJ_genomic_CDR3_segs()

        # If both VDJ, VJ loaders gave an exception, raise custom exception.
        except Exception as err:
            raise ModelLoaderException(err)

    def is_vj(self):
        """Function to specify if the loaded model is a VJ model.

        Returns
        -------
        boolean
            Indicating if the model is VJ compliant.

        """
        return bool(self.type[0] and self.type[2] and not self.type[1])

    def is_vdj(self):
        """Function to specify if the loaded model is a VDJ model.

        Returns
        -------
        boolean
            Indicating if the model is VDJ compliant.

        """
        return bool(self.type[0] and self.type[1] and self.type[2])

    def get_genomic_data(self):
        """Getter function for collecting the GenomicData OLGA object.

        Returns
        -------
        GenomicDataVJ or GenomicDataVDJ OLGA object
            The genomic data object class for a VJ or V(D)J model.

        """
        return self.data

    def get_generative_model(self):
        """Getter function for collecting the GenerativeModel OLGA object.

        Returns
        -------
        GenerativeModelVJ or GenerativeModelVDJ OLGA object
            The IGoR generative model object class for a VJ or V(D)J model.

        """
        return self.model


def main():
    """Function to be called when file executed via terminal."""
    print(__doc__)


if __name__ == "__main__":
    main()
