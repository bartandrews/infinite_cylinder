# --- python imports
import argparse


#############################################################
# input_arguments (parse the input arguments for a program) #
#############################################################


def parse_input_arguments(program):

    parser = argparse.ArgumentParser(prog=program)
    prog = parser.add_argument_group("program sub-arguments")
    stem = parser.add_argument_group("stem sub-arguments")
    leaf = parser.add_argument_group("leaf sub-arguments")

    models = ["BosHofSqu1", "FerHofSqu1", "BosHofHex1", "FerHofHex1",
              "BosHofHex1Hex5", "FerHofHex1Hex5", "BosHofHex1Hex5Orbital", "FerHofHex1Hex5Orbital"]
    Vtypes = ["Coulomb", "Yukawa"]

    prog.add_argument("-thr", "--threads", type=int, default=1, help="number of threads")
    stem.add_argument("-mod", "--model", type=str, default="BosHofSqu1", choices=models, required=True,
                        help="name of model")
    stem.add_argument("-chi", "--chi_max", type=int, default=50, required=True, help="maximum MPS bond dimension")
    prog.add_argument("-u_p", "--use_pickle", type=bool, default=False, help="use a pickled state")
    prog.add_argument("-m_p", "--make_pickle", type=bool, default=False, help="make a pickled state")

    for i in range(1, 11):  # search up to 10th-NN hoppings for both t and tdash
        tdefault = 1 if i == 1 else 0
        leaf.add_argument(f"-t{i}", type=float, default=tdefault, help=f"{i}-NN hopping parameter")
        leaf.add_argument(f"-t{i}dash", type=float, default=0, help=f"{i}-NN subsidiary hopping parameter")

    if program == "kappa_flow":
        leaf.add_argument("-kappa_min", type=float, default=0, required=True, help="minimum kappa")
        leaf.add_argument("-kappa_max", type=float, default=1, required=True, help="maximum kappa")
        leaf.add_argument("-kappa_samp", type=int, default=11, required=True, help="number of kappa samples")

    if program == "U_flow":
        leaf.add_argument("-U_min", type=float, default=0, required=True, help="minimum onsite interaction strength")
        leaf.add_argument("-U_max", type=float, default=10, required=True, help="maximum onsite interaction strength")
        leaf.add_argument("-U_samp", type=int, default=11, required=True,
                            help="number of onsite interaction strength samples")
    else:
        leaf.add_argument("-U", type=float, default=0, help="onsite interaction strength")

    leaf.add_argument("-mu", type=float, default=0, help="chemical potential")

    if program == "V_flow":
        leaf.add_argument("-V_min", type=float, default=0, required=True, help="minimum offsite interaction strength")
        leaf.add_argument("-V_max", type=float, default=10, required=True,
                            help="maximum offsite interaction strength")
        leaf.add_argument("-V_samp", type=int, default=11, required=True,
                            help="number of offsite interaction strength samples")
    else:
        leaf.add_argument("-V", type=float, default=0, help="offsite interaction strength")
    leaf.add_argument("-Vtype", type=str, default="Coulomb", choices=Vtypes,
                        help="offsite interaction type")
    leaf.add_argument("-Vrange", type=float, default=0, choices=range(11),
                      help="offsite interaction range (in units of nearest neighbors)")

    leaf.add_argument("-n", nargs=2, type=int, default=[1, 8], required=True, help="filling of the MPS unit cell")
    leaf.add_argument("-nphi", nargs=2, type=int, default=[1, 4], required=True, help="flux density")
    leaf.add_argument("-LxMUC", type=int, default=1, required=True,
                        help="width of MPS unit cell (in units of magnetic unit cell)")
    leaf.add_argument("-Ly", type=int, default=4, required=True,
                        help="height of MPS unit cell (in units of lattice unit cell)")

    if program == "phi_flow":
        leaf.add_argument("-phi_min", type=float, default=0, required=True,
                            help="minimum value of external flux (in units of 2*pi)")
        leaf.add_argument("-phi_max", type=float, default=2, required=True,
                            help="maximum value of external flux (in units of 2*pi)")
        leaf.add_argument("-phi_samp", type=int, default=21, required=True, help="number of external flux samples")
    else:
        leaf.add_argument("-phi", type=float, default=0, help="external flux (in units of 2*pi)")

    leaf.add_argument("-tag", type=str, default="", help="tag to append to filename")

    args = vars(parser.parse_args())
    __check_input_arguments(program, args)

    prog_args, stem_args, leaf_args = dict(), dict(), dict()
    for prog_key in ['threads', 'use_pickle', 'make_pickle']:
        prog_args.update({prog_key: args.pop(prog_key, None)})
    for stem_key in ['model', 'chi_max']:
        stem_args.update({stem_key: args.pop(stem_key, None)})
    leaf_args = args

    return prog_args, stem_args, leaf_args


##########################################################################
# check_input_arguments (check the requirements for the input arguments) #
##########################################################################


def __check_input_arguments(program, args):

    if args['threads'] < 0:
        raise ValueError("threads needs to be positive.")

    if args['chi_max'] < 0:
        raise ValueError("chi_max needs to be positive.")

    if "flow" in program:
        if args[f"{program.replace('flow', 'min')}"] > args[f"{program.replace('flow', 'max')}"]:
            raise ValueError(f"{program.replace('flow', 'max')} has to be greater "
                             f"than {program.replace('flow', 'min')}.")

    if (args['V'] == 0 and args['Vrange'] != 0) or (args['V'] != 0 and args['Vrange'] == 0):
        raise ValueError("Cannot have zero interaction over a finite range, or a finite interaction over zero range.")

    if args['n'][0] < 0 or args['n'][1] < 0:
        raise ValueError("n needs to have positive entries.")

    if args['nphi'][0] < 0 or args['nphi'][1] < 0:
        raise ValueError("nphi needs to have positive entries.")

    if args['LxMUC'] < 0 or args['Ly'] < 0:
        raise ValueError("LxMUX and Ly need to be positive.")

    return
