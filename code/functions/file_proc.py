# --- python imports
import sys
import numpy as np


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):
    if model not in ['BosonicHofstadter', 'FermionicHofstadter',
                     'BosonicHex1', 'FermionicHex1',
                     'BosonicHex1Hex5', 'FermionicHex1Hex5',
                     'BosonicHex1Hex5Orbital', 'FermionicHex1Hex5Orbital']:
        sys.exit('Error: Unknown model.')

    stem = f"{tool}_{model}_chi_{chi_max}_"

    return stem


######################################################################
# Logger (writes stdout and stderr both to the screen and to a file) #
######################################################################


class Logger(object):
    def __init__(self, flow, model, leaf):
        self.terminal = sys.stdout or sys.stderr
        self.log = open(f"data/output/{flow}/{model}/output_{flow}_{model}_{leaf}", 'w', buffering=1)

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
    if "Hofstadter" in model:
        LylB = Ly * np.sqrt(2 * np.pi * (pvalue / qvalue))
    elif "Hex" in model:
        LylB = Ly * np.sqrt((4 * np.pi * (pvalue / qvalue)) / np.sqrt(3))
    else:
        sys.exit("Error: Unknown model for the print_LylB_headings function.")

    if ndvalue == nd_min and qvalue == q_min and Ly == Ly_min:
        data_line = f"LylB={LylB:.15f}\n"
    else:
        data_line = f"\n\nLylB={LylB:.15f}\n"

    return data_line, LylB
