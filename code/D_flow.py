# --- python imports
import time
import sys
import numpy as np
# from scipy.linalg import expm
# from tenpy.linalg.np_conserved import expm
import tenpy.linalg.np_conserved as npc
# --- TeNPy imports
import tenpy.tools.process as prc
from tenpy.networks.mps import TransferMatrix
from tenpy.networks.site import SpinSite
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_args as fa


def spatial_inversion(self):
    """Perform a spatial inversion along the MPS.
    Exchanges the first with the last tensor and so on,
    i.e., exchange site `i` with site ``L-1 - i``.
    This is equivalent to a mirror/reflection with the bond left of L/2 (even L) or the site
    (L-1)/2 (odd L) as a fixpoint.
    For infinite MPS, the bond between MPS unit cells is another fix point.
    """
    self.sites = self.sites[::-1]
    self.form = [(f if f is None else (f[1], f[0])) for f in self.form[::-1]]
    self._B = [
        B.replace_labels(['vL', 'vR'], ['vR', 'vL']).transpose(self._B_labels)
        for B in self._B[::-1]
    ]
    self._S = self._S[::-1]
    self.test_sanity()
    return self


def my_D_flow(path_flag, threads, model, chi_max, ham_params):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("D_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("D_flow", path, model, chi_max, leaf)

    tools = ["ent_D_flow", "O_sym_flow"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ####################################################################################################################

    for t2 in np.linspace(ham_params['D_min'], ham_params['D_max'], ham_params['D_samp']):

        ham_params.update(t2=t2)
        state_data = fd.my_iDMRG_pickle("D_flow", path, model, chi_max, ham_params, run=True)
        psi = state_data['psi']

        ##############
        # ent_D_flow #
        ##############

        SvN = psi.entanglement_entropy()[0]

        data_line = f"{t2:.15f}\t{SvN:.15f}"
        print(data_line)
        data['ent_D_flow'].write(data_line + "\n")

        ##############
        # O_sym_flow #
        ##############

        site = SpinSite(S=1, conserve=None)

        # O_Z2xZ2 ######################################################################################################

        # compute Ux
        sx = site.get_op('Sx')
        scaled_sx = sx.copy()
        scaled_sx.iscale_prefactor(1j * np.pi)
        Rx = npc.expm(scaled_sx)
        psi_t = psi.copy()
        for i in range(psi.L):
            psi_t.apply_local_op(i, Rx, renormalize=True)
        TM = TransferMatrix(psi, psi_t, transpose=False, charge_sector=None)
        eta_x, w_x = TM.eigenvectors()
        Ux = w_x[0].split_legs().to_ndarray().conj().T
        print("eta_x, Ux.shape = ", eta_x, Ux.shape)

        # compute Uz
        sz = site.get_op('Sz')
        scaled_sz = sz.copy()
        scaled_sz.iscale_prefactor(1j * np.pi)
        Rz = npc.expm(scaled_sz)
        psi_t = psi.copy()
        for i in range(psi.L):
            psi_t.apply_local_op(i, Rz, renormalize=True)
        TM = TransferMatrix(psi, psi_t, transpose=False, charge_sector=None)
        eta_z, w_z = TM.eigenvectors()
        Uz = w_z[0].split_legs().to_ndarray().conj().T
        print("eta_z, Uz.shape = ", eta_z, Uz.shape)

        if Ux.shape[0] != Uz.shape[0] or Ux.shape[1] != Uz.shape[1] or Ux.shape[0] != Uz.shape[1]:
            raise ValueError("Ux and Uz bond dimension mismatch.")

        if abs(eta_x) > 0.99 and abs(eta_z) > 0.99:
            O_Z2xZ2 = Ux.shape[0] * np.real(np.trace(Ux @ Uz @ Ux.conj().T @ Uz.conj().T))  # spurious chi^2 factor
        else:
            O_Z2xZ2 = 0

        # O_I ##########################################################################################################

        # compute UI
        psi_t = psi.copy()
        psi_t = spatial_inversion(psi_t)
        TM = TransferMatrix(psi, psi_t, transpose=False, charge_sector=None)
        eta_I, w_I = TM.eigenvectors()
        UI = w_I[0].split_legs().to_ndarray().conj().T
        print("eta_I, UI.shape = ", eta_I, UI.shape)

        if UI.shape[0] != UI.shape[1]:
            raise ValueError("UI bond dimension mismatch.")

        if abs(eta_I) > 0.99:
            O_I = np.real(np.trace(UI @ UI.conj()))  # spurious chi factor
        else:
            O_I = 0

        # O_TR #########################################################################################################

        # compute Uy
        sy = site.get_op('Sy')
        scaled_sy = sy.copy()
        scaled_sy.iscale_prefactor(1j * np.pi)
        Ry = npc.expm(scaled_sy)
        psi_t = psi.copy()
        for i in range(psi.L):
            psi_t.apply_local_op(i, Ry, renormalize=True)
        TM = TransferMatrix(psi, psi_t, transpose=False, charge_sector=None)
        eta_y, w_y = TM.eigenvectors()
        Uy = w_y[0].split_legs().to_ndarray().conj().T
        print("eta_y, Uy.shape = ", eta_y, Uy.shape)

        if Uy.shape[0] != Uy.shape[1]:
            raise ValueError("Uy bond dimension mismatch.")

        if abs(eta_y) > 0.99:
            O_TR = np.real(np.trace(Uy @ Uy.conj()))  # spurious chi factor
        else:
            O_TR = 0

        ################################################################################################################

        data_line = f"{t2:.15f}\t{O_Z2xZ2:.15f}\t{O_I:.15f}\t{O_TR:.15f}"
        print(data_line)
        data['O_sym_flow'].write(data_line + "\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("D_flow")

    my_D_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args)
