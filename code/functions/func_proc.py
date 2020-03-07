# --- python imports
import sys
import os
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):

    if model not in ['BosHofSqu1', 'FerHofSqu1',
                     'BosHofHex1', 'FerHofHex1',
                     'BosHofHex1Hex5', 'FerHofHex1Hex5',
                     'BosHofHex1Hex5Orbital', 'FerHofHex1Hex5Orbital']:
        raise ValueError('Unknown model for the file_name_stem function.')

    stem = f"{tool}_{model}_chi_{chi_max}_"

    return stem


##################################################################
# prepare_output_files (creates the output directories and files) #
##################################################################


def prepare_output_files(tools, model, chi_max, leaf, chi_max_K=0):

    stem, file, data = [dict()]*3
    for tool in tools:
        stem.update({tool: file_name_stem(tool, model, chi_max)})
        os.makedirs(f"data/{tool}/{model}/", exist_ok=True)
        if tool == 'ent_spec_mom':
            file.update({tool: f"data/{tool}/{model}/" + stem[tool].replace(" ", "_") + f"chi_K_{chi_max_K}_" + leaf})
        else:
            file.update({tool: f"data/{tool}/{model}/" + stem[tool].replace(" ", "_") + leaf})
        open(file[tool], "w")
        data[tool] = open(file[tool], "a", buffering=1)

    return data


######################################################################
# Logger (writes stdout and stderr both to the screen and to a file) #
######################################################################


class Logger(object):
    def __init__(self, flow, model, leaf):
        self.terminal = sys.stdout or sys.stderr
        os.makedirs(f"logs/{flow}/{model}/", exist_ok=True)
        self.log = open(f"logs/{flow}/{model}/log_{flow}_{model}_{leaf}", 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


######################################################################################
# print_LylB_headings (prints the Ly/lB ratio to file if there is a range of values) #
######################################################################################


def print_LylB_headings(model, Ly, ndvalue, nd_min, pvalue, qvalue, q_min, Ly_min):
    if "Squ" in model:
        LylB = Ly * np.sqrt(2 * np.pi * (pvalue / qvalue))
    elif "Hex" in model:
        LylB = Ly * np.sqrt((4 * np.pi * (pvalue / qvalue)) / np.sqrt(3))
    else:
        raise ValueError("Unknown model for the print_LylB_headings function.")

    if ndvalue == nd_min and qvalue == q_min and Ly == Ly_min:
        data_line = f"LylB={LylB:.15f}\n"
    else:
        data_line = f"\n\nLylB={LylB:.15f}\n"

    return data_line, LylB


if __name__ == '__main__':

    models = []
    for m in CouplingMPOModel.__subclasses__():
        models.append(m.__name__)
    print(models)
