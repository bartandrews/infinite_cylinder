from tenpy.networks.mps import MPS

from models.hofstadter_old import BosonicHofstadterModel, FermionicHofstadterModel
from models.hex_1 import BosonicHex1Model, FermionicHex1Model
from models.hex_1_hex_5 import BosonicHex1Hex5Model, FermionicHex1Hex5Model
from models.hex_1_hex_5_orbital import BosonicHex1Hex5OrbitalModel, FermionicHex1Hex5OrbitalModel

from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH

import sys
import pickle


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):

    if model not in ['BosonicHofstadter', 'FermionicHofstadter',
                     'BosonicHex1', 'FermionicHex1',
                     'BosonicHex1Hex5', 'FermionicHex1Hex5',
                     'BosonicHex1Hex5Orbital', 'FermionicHex1Hex5Orbital']:
        sys.exit('Error: Unknown model.')

    stem = ("%s_%s_chi_%s_"
            % (tool, model, chi_max))

    return stem


################################################
# select_initial_psi (selects the initial psi) #
################################################


def select_initial_psi(model, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly):

    if model == 'BosonicHofstadter':
        if int((nnvalue / ndvalue) / (pvalue / qvalue)) == 1:
            if (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 1, 0, 0, 0,
                                 1, 0, 0, 0, 1, 0, 0, 0]
        else:
            if (nnvalue, ndvalue) == (1, 6) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (4, 9) and (pvalue, qvalue) == (4, 3) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 1, 0, 1, 0, 1, 0, 0,
                                 1, 0, 1, 0, 1, 0, 1, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
    elif model == 'FermionicHofstadter':
        if (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 9):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 9):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 9):
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif model == 'BosonicHex1' or model == 'BosonicHex1Hex5':
        if int((nnvalue/ndvalue)/(pvalue/qvalue))==1:
            if (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 6) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
        else:
            if (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                    product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 4):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
    elif model == 'FermionicHex1' or model == 'FermionicHex1Hex5':
        if int((nnvalue / ndvalue) / (pvalue / qvalue)) == 1:
            if (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
        else:
            if (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 5):
                # polarized2 (gives same results as unpolarized2 but faster)
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 10):
                # polarized2 (gives same results as unpolarized2 but faster)
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 20) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 20) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 25) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 25) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (2, 21) and (pvalue, qvalue) == (2, 7) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (2, 21) and (pvalue, qvalue) == (2, 7) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (8, 27) and (pvalue, qvalue) == (8, 9) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (8, 27) and (pvalue, qvalue) == (8, 9) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (10, 33) and (pvalue, qvalue) == (10, 11) and (Lx_MUC, Ly) == (1, 6):
                product_state = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (10, 33) and (pvalue, qvalue) == (10, 11) and (Lx_MUC, Ly) == (1, 9):
                product_state = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
    elif model == 'BosonicHex1Hex5Orbital':
        if int((nnvalue / ndvalue) / (pvalue / qvalue)) == 1:
            if (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 6) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
        else:
            if (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 4):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
    elif model == 'FermionicHex1Hex5Orbital':
        if int((nnvalue / ndvalue) / (pvalue / qvalue)) == 1:
            if (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 5):  # working IQH
                # polarized2 (expected to give same results as unpolarized2 but faster)
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0,
                                 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0]
                # unpolarized2
                # product_state = ['full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0,
                #                  'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0]
                # full1
                # product_state = ['full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0,
                #                  0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0]
                # polarized2d
                # product_state = ['full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0,
                #                  'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0, 'full_x empty_y', 'full_x empty_y', 0]
                # unpolarized2d
                # product_state = ['full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 0]
                # full1d
                # product_state = ['full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0,
                #                  'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0]
            elif (nnvalue, ndvalue) == (6, 7) and (pvalue, qvalue) == (6, 7) and (Lx_MUC, Ly) == (1, 6):
                # fill2
                # product_state = ['full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0,
                #                  'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0, 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 'full_x empty_y', 'empty_x full_y', 0]
                # fill3
                product_state = ['full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0,
                                 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0,
                                 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0,
                                 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0,
                                 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0,
                                 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0, 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 'full_x full_y', 0]
            elif (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                # polarized1 (original)
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0,
                #                  'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0]
                # polarized2 (double_fill)
                product_state = ['full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0,
                                 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0, 'full_x empty_y', 0, 0]
                # unpolarized1
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 'empty_x full_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0,
                #                  'empty_x full_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'empty_x full_y', 0, 0, 0, 0, 0]
                # unpolarized2
                # product_state = ['full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0,
                #                  'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0, 'full_x empty_y', 0, 0, 'empty_x full_y', 0, 0]
                # full1 (double_fill2)
                # product_state = ['full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0,
                #                  'full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0]
                # full2 (double_fill3)
                # product_state = ['full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0,
                #                  'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0, 'full_x full_y', 0, 0]
            elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
        else:
            if (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 5):  # working FQH
                # polarized2 (gives same results as unpolarized2 but faster)
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                # unpolarized2
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'empty_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                # full1
                # product_state = ['full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 10):
                # polarized2 (gives same results as unpolarized2 but faster)
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 20) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 25) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 5):
                # polarized2
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 6):
                # polarized1 (original)
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                # polarized2
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0]
                # unpolarized1
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'empty_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                # unpolarized2
                # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'empty_x full_y', 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'empty_x full_y', 0, 0, 0, 0, 0, 0, 0, 0]
                # full1
                product_state = ['full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                # full2
                # product_state = ['full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0,
                #                  'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx_MUC, Ly) == (1, 9):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx_MUC, Ly) == (1, 9):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 6):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx_MUC, Ly) == (1, 9):
                product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            else:
                sys.exit('Error: Unknown initial_state configuration.')
    else:
        sys.exit('Error: Unknown initial_state.')

    return product_state


#############################################################################
# define_iDMRG_model (defines the changable parameters for the iDMRG model) #
#############################################################################

def define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):

    if model == 'BosonicHofstadter':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHofstadterModel(model_params)

    elif model == 'FermionicHofstadter':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, V=10,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHofstadterModel(model_params)

    elif model == 'BosonicHex1':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Model(model_params)

    elif model == 'FermionicHex1':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1Model(model_params)

    elif model == 'BosonicHex1Hex5':
        model_params = dict(conserve='N', t1=t1, t2=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Hex5Model(model_params)

    elif model == 'FermionicHex1Hex5':
        model_params = dict(conserve='N', t1=t1, t2=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1Hex5Model(model_params)

    elif model == 'BosonicHex1Hex5Orbital':
        model_params = dict(conserve='N', t1=t1, t2=t2, t2dash=t2dash, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Hex5OrbitalModel(model_params)

    elif model == 'FermionicHex1Hex5Orbital':
        model_params = dict(conserve='N', t1=t1, t2=t2, t2dash=t2dash, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                            Lx_MUC=Lx_MUC, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1Hex5OrbitalModel(model_params)

    return M


#######################################################
# define iDMRG (used when we want to reuse the state) #
#######################################################


def define_iDMRG_engine_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("engine", model, chi_max)
        pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_MUC_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:
        with open(pickle_file, 'rb') as file1:
            engine = pickle.load(file1)
    else:
        engine = define_iDMRG_engine(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump(engine, file2)

    return engine


def define_iDMRG_engine(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):

    M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
    product_state = select_initial_psi(model, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly)
    print(product_state)
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-6,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        'max_sweeps': 1000,
        'verbose': 1,
        'N_sweeps_check': 10,
        # 'diag_method': 'lanczos'
    }

    # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
    engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)

    return engine


##########################################################
# run iDMRG (used when we want to recalculate the state) #
##########################################################


def run_iDMRG_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("E_psi_M", model, chi_max)
        pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_MUC_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:
        with open(pickle_file, 'rb') as file1:
            [E, psi, M] = pickle.load(file1)
    else:
        (E, psi, M) = run_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump([E, psi, M], file2)

    return E, psi, M


def run_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):

    M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
    product_state = select_initial_psi(model, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly)
    print(product_state)
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-6,
        'max_S_err': 1.e-6,
        'max_sweeps': 1000,
        'verbose': 1,
        'N_sweeps_check': 10,
        # 'diag_method': 'lanczos'
    }

    info = dmrg.run(psi, M, dmrg_params)
    E = info['E']

    return E, psi, M
