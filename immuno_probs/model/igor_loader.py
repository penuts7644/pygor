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


"""Contains IgorLoader class for loading in a IGoR model files."""


import olga.load_model as olga_load_model


class IgorLoader(object):
    """Loads in an IGoR model as well as corresponding with CDR3 anchor files.

    Parameters
    ----------
    model_type : str
        The type of the input model, either 'VJ' or 'VDJ'.
    model_params : str
        A file path location for the IGoR parameters model file.
    model_marginals : str
        A file path location for the IGoR marginals model file.

    Methods
    -------
    set_anchor(gene, file)
        Set the model's CDR3 V or J gene anchors.
    initialize_model()
        Initializes the model including the anchor files.
    get_type()
        Returns the type of the model ('VDJ' or 'VJ').
    get_genomic_data()
        Return the OLGA's GenomicData object.
    get_generative_model()
        Return the OLGA's GenerativeModel object.

    """
    def __init__(self, model_type, model_params, model_marginals):
        super(IgorLoader, self).__init__()
        self.type = self._check_type(model_type, model_marginals)
        self.data = self._load_params(model_params)
        self.model = self._load_model(model_marginals)
        self.params = model_params
        self.v_anchors = None
        self.j_anchors = None

    @staticmethod
    def _check_type(model_type, model_marginals):
        """Private function to check the model marginals file for possible D gene attributes. If these are found, the model
        classifies as VDJ else, VJ.

        Parameters
        ----------
        model_type : str
            The type of the input model: alpha, beta, light or heavy.
        model_marginals : str
            A file path location for the IGoR model marginals file.

        Returns
        -------
        str
            Specifying the model type, either 'VJ' or 'VDJ'.

        Raises
        ------
        TypeError
            When the model marginals are not compliant to the given model type.

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
        if (v_choice and d_gene and j_choice) and (model_type in ['beta', 'heavy']):
            return 'VDJ'
        if (v_choice and j_choice and not d_gene) and (model_type in ['alpha', 'light']):
            return 'VJ'
        raise TypeError("Model is not compliant to the given type: '{}''".format(model_type))

    def _load_params(self, model_params):
        """Private function for loading in the genomic parameter data for the IGoR model.

        Parameters
        ----------
        model_params : str
            A file path location for the IGoR parameters model file.

        Returns
        -------
        GenomicDataVJ or GenomicDataVDJ OLGA object
            The genomic data object class for a VJ or VDJ model.

        Raises
        ------
        TypeError
            When the model input data cannot be loaded in as either a VJ or VDJ model.
        OSError
            When OLGA produces and system error with the input data.

        """
        # Try to load the genomic data model for VDJ.
        try:
            genomic_data = None
            if self.type == 'VDJ':
                genomic_data = olga_load_model.GenomicDataVDJ()
                genomic_data.genD = olga_load_model.read_igor_D_gene_parameters(model_params)
            elif self.type == 'VJ':
                genomic_data = olga_load_model.GenomicDataVJ()
            else:
                raise TypeError("Model genomic data could not be loaded as 'VDJ' or 'VJ' type")

            # Load the remainder of the data for the VJ model and return.
            genomic_data.genV = olga_load_model.read_igor_V_gene_parameters(model_params)
            genomic_data.genJ = olga_load_model.read_igor_J_gene_parameters(model_params)
            return genomic_data

        except Exception as err:
            raise OSError(err)

    def _load_model(self, model_marginals):
        """Private function for loading in the IGoR model marginals.

        Parameters
        ----------
        model_marginals : str
            A file path location for the IGoR marginals model file.

        Returns
        -------
        GenerativeModelVJ or GenerativeModelVDJ OLGA object
            The IGoR generative model object class for a VJ or VDJ model.

        Raises
        ------
        TypeError
            When the model input data cannot be loaded in as either a VJ or VDJ model.
        OSError
            When OLGA produces and system error with the input data.

        """
        # Try to create the GenerativeModel object for VDJ or VJ.
        try:
            generative_model = None
            if self.type == "VDJ":
                generative_model = olga_load_model.GenerativeModelVDJ()
            elif self.type == "VJ":
                generative_model = olga_load_model.GenerativeModelVJ()
            else:
                raise TypeError("Generative model could not be loaded as 'VDJ' or 'VJ' type")

            # Load the generative VDJ or VJ model marginals and return.
            generative_model.load_and_process_igor_model(model_marginals)
            return generative_model

        except Exception as err:
            raise OSError(err)

    def set_anchor(self, gene, file):
        """Sets the CDR3 anchor file path for a given gene.

        Parameters
        ----------
        gene : str
            A gene identifier, either 'V' or 'J', specifying the alignement's origin gene.
        file : str
            File path location for the CDR3 anchor positions for the given gene.

        Returns
        -------
        str
            The input gene character value after passing the validation test.

        Raises
        ------
        ValueError
            When the given gene character does not equal 'V' or 'J'.

        """
        gene = gene.upper()
        if gene == "V":
            self.v_anchors = file
        elif gene == "J":
            self.j_anchors = file
        else:
            raise ValueError("Gene identifier should be either 'V' or 'J'", gene)

    def initialize_model(self):
        """Initializes the model data with CDR3 anchor position data files.

        Raises
        ------
        TypeError
            When the model input data cannot be loaded in as either a VJ or VDJ model.
        OSError
            When OLGA produces and system error with the input data.

        Notes
        -----
        This function uses the model marginals and parameters data that has already been loaded through trhe class constructor
        and updates this model object with the CDR3 anchor files.

        """

        # Try to load the anchor files into the data model for VDJ.
        try:
            self.data.anchor_and_curate_genV_and_genJ(self.v_anchors, self.j_anchors)
            if self.type == "VDJ":
                self.data.read_VDJ_palindrome_parameters(self.params)
                self.data.generate_cutD_genomic_CDR3_segs()
            elif self.type == "VJ":
                self.data.read_igor_VJ_palindrome_parameters(self.params)
            else:
                raise TypeError("Model type could not be identified as 'VDJ' or 'VJ'")

            # Load the remainder of the data model
            self.data.generate_cutV_genomic_CDR3_segs()
            self.data.generate_cutJ_genomic_CDR3_segs()

        except Exception as err:
            raise OSError(err)

    def get_type(self):
        """Collects and returns the type of the model.

        Returns
        -------
        str
            Indicates the model type, either 'VDJ' or 'VJ'.

        """
        return self.type

    def get_genomic_data(self):
        """Collects and returns the GenomicData OLGA object.

        Returns
        -------
        GenomicDataVJ or GenomicDataVDJ OLGA object
            The genomic data object class for a VJ or VDJ model.

        """
        return self.data

    def get_generative_model(self):
        """Collects and returns the GenerativeModel OLGA object.

        Returns
        -------
        GenerativeModelVJ or GenerativeModelVDJ OLGA object
            The IGoR generative model object class for a VJ or VDJ model.

        """
        return self.model
