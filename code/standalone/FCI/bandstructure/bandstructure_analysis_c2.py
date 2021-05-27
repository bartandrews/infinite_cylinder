# --- python imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import csv
import matplotlib.gridspec as gridspec

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

a0 = 1  # lattice constant
t = 1  # hopping amplitude


def define_unit_cell(model_val, q_val=4, C_val=1):

    if 'Squ' in model_val:
        if model_val is 'HofSqu1':
            num_bands_val = q_val
            # lattice vectors
            a1 = a0 * np.array([num_bands_val, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / num_bands_val, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]
        elif model_val in ['HalSquCN']:
            if C_val == 1:
                num_bands_val = 2
                lx = 2  # rectangular unit cell
            else:
                num_bands_val = C_val
                lx = 1  # square unit cell
            # lattice vectors
            a1 = a0 * np.array([lx, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / lx, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            if C_val == 1:  # rectangular BZ
                GA = np.array([0, 0])
                Y = np.array([0, 0.5])
                S = np.array([0.5, 0.5])
                X = np.array([0.5, 0])
                sym_points_val = [GA, Y, S, X]
            else:  # square BZ
                GA = np.array([0, 0])
                M = np.array([0.5, 0.5])
                X = np.array([0.5, 0])
                sym_points_val = [GA, M, X]
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif 'Tri' in model_val:
        if model_val is 'HalTriC3':
            num_bands_val = 2
            # lattice vectors
            a1 = a0 * np.array([1, 0])
            a2 = a0 * np.array([1/2, np.sqrt(3)/2])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1, -1/np.sqrt(3)])
            b2 = (2. * np.pi) / a0 * np.array([0, 2/np.sqrt(3)])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            sym_points_val = []  # to implement
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif 'Hex' in model_val or model_val is 'graphene':
        num_bands_val = 2
        # lattice vectors
        a1 = (a0 / 2) * np.array([3, np.sqrt(3)])
        a2 = (a0 / 2) * np.array([3, -np.sqrt(3)])
        avec_val = np.vstack((a1, a2))
        # reciprocal lattice vectors
        b1 = (2. * np.pi) / (3 * a0) * np.array([1, np.sqrt(3)])
        b2 = (2. * np.pi) / (3 * a0) * np.array([1, -np.sqrt(3)])
        bvec_val = np.vstack((b1, b2))
        # symmetry points
        K1 = np.array([2/3, 1/3])
        GA = np.array([0., 0.])
        MM = np.array([0.5, 0.5])
        K2 = np.array([1/3, 2/3])
        sym_points_val = [K1, GA, MM, K2]
    else:
        return ValueError("Requested lattice cannot be read from model name.")

    return num_bands_val, avec_val, bvec_val, sym_points_val


def hamiltonian(model_val, k_val, num_bands_val, avec_val, p_val=1, C_val=1, tx_factor_val=1):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

    if 'Squ' in model_val:

        # nearest neighbors
        delta = np.zeros((2, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([0, 1])

        if model_val is 'HofSqu1':

            q_val = num_bands_val
            nphi = p_val/q_val

            def h(k_val_val, m_val):
                return 2 * np.cos(2 * np.pi * nphi * m_val + k_val_val[1] * a0)

            for n in range(q_val):
                Hamiltonian[n][n] = t * h(k_val, n)

            for n in range(q_val-1):
                Hamiltonian[n][n+1] = tx_factor_val*t * np.exp(+1j*k_val[0]*a0)
                Hamiltonian[n+1][n] = tx_factor_val*t * np.exp(-1j*k_val[0]*a0)

            Hamiltonian[0][q_val-1] = tx_factor_val*t * np.exp(-1j*k_val[0]*a0)
            Hamiltonian[q_val-1][0] = tx_factor_val*t * np.exp(+1j*k_val[0]*a0)

        elif model_val in ['HalSquCN']:
            if C_val in [1, 2]:
                tphi = tx_factor_val * t * np.exp(1j * np.pi / 4)
                tdash = t / (2 + np.sqrt(2))
                tddash = t / (2 + 2 * np.sqrt(2))

                # first-nearest neighbors
                firstNN = np.zeros((4, 2))
                # positive direction for BdA
                firstNN[0, :] = delta[1, :]  # up
                firstNN[1, :] = -delta[1, :]  # down
                # positive direction for AdB
                firstNN[2, :] = -delta[0, :]  # left
                firstNN[3, :] = delta[0, :]  # right

                # second-nearest neighbors
                secondNN = np.zeros((2, 2))
                secondNN[0, :] = (delta[0, :] + delta[1, :])
                secondNN[1, :] = (delta[0, :] - delta[1, :])

                # third-nearest neighbors
                thirdNN = np.zeros((2, 2))
                thirdNN[0, :] = 2 * delta[0, :]
                thirdNN[1, :] = 2 * delta[1, :]

                f_AdB = 0
                for m in range(0, 2):
                    f_AdB += tphi * np.exp(1j * k_val.dot(firstNN[m, :]))
                f_BdA = 0
                for m in range(2, 4):
                    f_BdA += tphi * np.exp(1j * k_val.dot(firstNN[m, :]))

                fdash_A = tdash * np.exp(1j * k_val.dot(secondNN[0, :])) - tdash * np.exp(
                    1j * k_val.dot(secondNN[1, :]))
                fdash_B = -tdash * np.exp(1j * k_val.dot(secondNN[0, :])) + tdash * np.exp(
                    1j * k_val.dot(secondNN[1, :]))

                fddash = 0
                for m in range(0, 2):
                    fddash += tddash * np.exp(1j * k_val.dot(thirdNN[m, :]))

                Hamiltonian[0][0] = fdash_A + fddash
                Hamiltonian[0][1] = f_AdB
                Hamiltonian[1][0] = f_BdA
                Hamiltonian[1][1] = fdash_B + fddash

                Hamiltonian += Hamiltonian.conj().transpose()
            else:
                N = C_val
                t1 = tx_factor_val * t
                if C_val < 100:  # need to include the t2 term otherwise cannot tune gap-to-width ratio
                    t2 = - t / np.sqrt(N)
                else:
                    t2 = 0
                phi = np.pi / N

                # first-nearest neighbors
                firstNN = np.zeros((2, 2))
                # positive phase direction for inter-orbital creation
                firstNN[0, :] = -delta[1, :]  # down
                # neutral direction
                firstNN[1, :] = delta[0, :]  # right

                # second-nearest neighbors
                secondNN = np.zeros((3, 2))
                # positive phase direction for inter-orbital creation
                secondNN[0, :] = (delta[0, :] - delta[1, :])  # bottom right
                # positive phase direction for same-orbital creation
                secondNN[1, :] = (-delta[0, :] - delta[1, :])  # bottom left
                # negative phase direction for same orbital creation
                secondNN[2, :] = (delta[0, :] + delta[1, :])  # top right

                for a in range(N):
                    Hamiltonian[(a+1) % N][a] = t1 * (np.exp(1j * k_val.dot(firstNN[1, :])) + np.exp(1j * (k_val.dot(firstNN[0, :]) + 2*(a+1) * phi)))
                    Hamiltonian[a][a] = t2 * np.exp(1j * (k_val.dot(secondNN[2, :]) - (2*(a+1)-1) * phi))
                    Hamiltonian[(a+2) % N][a] = t2 * np.exp(1j * (k_val.dot(secondNN[0, :]) + (2*(a+1)+1) * phi))

                Hamiltonian += Hamiltonian.conj().transpose()
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    elif 'Tri' in model_val:

        # nearest neighbors
        delta = np.zeros((3, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([-1/2, np.sqrt(3)/2])
        delta[2, :] = a0 * np.array([-1/2, -np.sqrt(3)/2])

        if model_val is 'HalTriC3':
            t1 = t
            t2 = 0.39 * t * np.exp(1j * np.pi/2)
            t3 = -0.34 * t

            # first-nearest neighbors
            firstNN = np.zeros((3, 2))
            firstNN[0, :] = delta[0, :]
            firstNN[1, :] = delta[1, :]
            firstNN[2, :] = delta[2, :]

            # second-nearest neighbors
            secondNN = np.zeros((6, 2))
            # positive direction for A
            secondNN[0, :] = a0 * np.array([0, np.sqrt(3)])
            secondNN[1, :] = a0 * np.array([1.5, -np.sqrt(3)/2])
            secondNN[2, :] = a0 * np.array([-1.5, -np.sqrt(3)/2])
            # positive direction for B
            secondNN[3, :] = a0 * np.array([1.5, np.sqrt(3)/2])
            secondNN[4, :] = a0 * np.array([-1.5, np.sqrt(3)/2])
            secondNN[5, :] = a0 * np.array([0, -np.sqrt(3)])

            # third-nearest neighbors
            thirdNN = np.zeros((3, 2))
            thirdNN[0, :] = 2 * delta[0, :]
            thirdNN[1, :] = 2 * delta[1, :]
            thirdNN[2, :] = 2 * delta[2, :]

            f1 = 0
            for m in range(0, 3):
                f1 += -t1 * np.exp(1j * k_val.dot(firstNN[m, :]))

            f2_A = 0
            for m in range(0, 3):
                f2_A += -t2 * np.exp(1j * k_val.dot(secondNN[m, :]))
            f2_B = 0
            for m in range(3, 6):
                f2_B += -t2 * np.exp(1j * k_val.dot(secondNN[m, :]))

            f3 = 0
            for m in range(0, 3):
                f3 += -t3 * np.exp(1j * k_val.dot(thirdNN[m, :]))

            Hamiltonian[0][0] = f2_A
            Hamiltonian[0][1] = f1
            Hamiltonian[1][0] = f3
            Hamiltonian[1][1] = f2_B

            Hamiltonian += Hamiltonian.conj().transpose()
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    elif 'Hex' in model_val or model_val is 'graphene':

        # nearest neighbors
        delta = np.zeros((3, 2))
        delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
        delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
        delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

        if model_val is 'graphene':
            f = 0
            for m in range(3):
                f += -t * np.exp(1j * k_val.dot(delta[m, :]))

            Hamiltonian[0][0] = 0
            Hamiltonian[0][1] = f
            Hamiltonian[1][0] = np.conj(f)
            Hamiltonian[1][1] = 0
        elif model_val is 'HalHexC1':
            t2 = t * np.sqrt(129) / 36
            phi = np.arccos(3*np.sqrt(3/43))

            # next-nearest neighbors
            secondNN = np.zeros((6, 2))
            # positive direction for A sites / negative direction for B sites
            secondNN[0, :] = avec_val[0, :]
            secondNN[1, :] = -avec_val[0, :] + avec_val[1, :]
            secondNN[2, :] = -avec_val[1, :]
            # negative direction for A sites / positive direction for B sites
            secondNN[3, :] = avec_val[1, :]
            secondNN[4, :] = -avec_val[0, :]
            secondNN[5, :] = avec_val[0, :] - avec_val[1, :]

            f = 0
            for m in range(0, 3):
                f += t * np.exp(1j * k_val.dot(delta[m, :]))
            f1 = 0
            for m in range(0, 3):
                f1 += t2 * np.exp(1j * k_val.dot(secondNN[m, :]))
            f2 = 0
            for m in range(3, 6):
                f2 += t2 * np.exp(1j * k_val.dot(secondNN[m, :]))

            Hamiltonian[0][0] = np.exp(1j * phi) * f1 + np.exp(-1j * phi) * f2
            Hamiltonian[0][1] = f
            Hamiltonian[1][0] = np.conj(f)
            Hamiltonian[1][1] = np.exp(1j * phi) * f2 + np.exp(-1j * phi) * f1
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")
    else:
        return ValueError("Requested lattice cannot be read from model name.")

    return Hamiltonian


def berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    Berry_curv = - np.imag(np.log(np.conj(eigvec).dot(eigvec_alpha) * np.conj(eigvec_alpha).dot(eigvec_alpha_beta)
                           * np.conj(eigvec_alpha_beta).dot(eigvec_beta) * np.conj(eigvec_beta).dot(eigvec)))

    return Berry_curv


if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 2.5))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.5, width_ratios=[1.3, 1])
    left_cell = outer_grid[0, 0]
    right_cell = outer_grid[0, 1]
    right_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, right_cell, hspace=0)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # initialization ###################################################################################################
    num_samples = 101
    model = 'HofSqu1'  # (HofSqu1, HalSquCN), (HalTriC3), (graphene, HalHexC1)
    mining = False  # data mining mode for 2D band structures
    if not mining:
        flag_3D = True  # choose between 3D or 2D band structure
        p, q = 3, 5  # for Hofstadter model only
        C = 1  # for HalSquCN model only
    else:
        tx_min, tx_max, tx_samp = 0, 1000, 1001
        flag_3D = False  # only works in 2D mode
        p = [1, 4, 4, 4, 4]
        q = [4, 7, 11, 15, 19]  # lists need to be the same length
        C = [1, 2, 3, 4, 5]
        if 'Hof' in model:
            path_to_file = "code/models/hofstadter/squ_1_hoppings_ratio.dat"
        elif 'Hal' in model:
            path_to_file = "code/models/haldane/squ_hoppings_ratio.dat"
        else:
            raise ValueError("Chosen model is not implemented for mining.")
    ####################################################################################################################

    if flag_3D:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q, C_val=C)

        # construct bands
        eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
        eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
        for band in range(num_bands):
            for idx_x in range(num_samples):
                frac_kx = idx_x / (num_samples-1)
                for idx_y in range(num_samples):
                    frac_ky = idx_y / (num_samples-1)
                    k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                    eigvals, eigvecs = np.linalg.eig(hamiltonian(model, k, num_bands, avec, p_val=p, C_val=C, tx_factor_val=1))
                    idx = np.argsort(eigvals)
                    eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                    eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

        # compute Chern numbers
        berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
        for band in range(num_bands):
            for idx_x in range(num_samples - 1):
                for idx_y in range(num_samples - 1):
                    berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y + 1])
        chern_numbers = np.zeros(num_bands)
        for band in range(num_bands):
            chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
            print(f"Chern number ({band}) = {chern_numbers[band]}")

        # construct figure
        ax = plt.subplot(left_cell, projection='3d')
        idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            # ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([0.5, 0.5, 1, 1]))
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky])
        ax.xaxis.labelpad = -1
        ax.yaxis.labelpad = -1
        ax.zaxis.labelpad = -5
        ax.set_xlabel('$k_x / |\mathbf{b}_x|$')
        ax.set_ylabel('$k_y / |\mathbf{b}_y|$')


        def custom(value, tick_number):

            if value == 0:
                return "$0$"
            elif value == num_samples - 1:
                return "$1$"
            else:
                return "${}$".format(round(value / (num_samples - 1), 2))


        ax.xaxis.set_major_formatter(plt.FuncFormatter(custom))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(custom))

        ax.set_xticks([0, 50, 100])
        ax.set_yticks([0, 50, 100])

        ax.set_zlabel('$E$')
        ax.set_title('$p=3$, $C=2$ $\Rightarrow\;n_\phi=\\frac{3}{5}$', pad=-5)
        ax.tick_params(axis='x', which='major', pad=-1.5)
        ax.tick_params(axis='y', which='major', pad=-1.5)
        ax.tick_params(axis='z', which='major', pad=-1.5)
        ax.view_init(elev=8, azim=135)
    elif mining:
        total_num_chern = len(C)
        ratio = np.zeros((total_num_chern, tx_samp))
        # loop through all chern numbers
        for num_chern in range(total_num_chern):
            for num_sample, tx_factor in enumerate(np.linspace(tx_min, tx_max, tx_samp)):
                # define unit cell
                num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q[num_chern], C_val=C[num_chern])
                # construct bands
                num_paths = len(sym_points)
                points_per_path = int(num_samples / num_paths)
                num_points = num_paths * points_per_path
                eigenvalues = np.zeros((num_bands, num_points))  # real
                count = 0
                for i in range(num_paths):
                    for j in range(points_per_path):
                        k = sym_points[i] + (sym_points[(i + 1) % num_paths] - sym_points[i]) * float(j) / float(points_per_path - 1)
                        k = np.matmul(k, bvec)
                        eigvals = np.linalg.eigvals(hamiltonian(model, k, num_bands, avec, p_val=p[num_chern], C_val=C[num_chern], tx_factor_val=int(tx_factor)))
                        idx = np.argsort(eigvals)
                        for band in range(num_bands):
                            eigenvalues[band, count] = np.real(eigvals[idx[band]])
                        count += 1
                band_gap = np.min(eigenvalues[1]) - np.max(eigenvalues[0])
                band_width = np.max(eigenvalues[0]) - np.min(eigenvalues[0])
                ratio[num_chern][num_sample] = band_gap/band_width

        with open(path_to_file, 'w') as file:
            for num_sample, tx_factor in enumerate(np.linspace(tx_min, tx_max, tx_samp)):
                file.write(f"{tx_factor}\t{ratio[0][num_sample]}\t{ratio[1][num_sample]}\t{ratio[2][num_sample]}\t{ratio[3][num_sample]}\t{ratio[4][num_sample]}\n")
    else:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q, C_val=C)

        # construct bands
        num_paths = len(sym_points)
        points_per_path = int(num_samples/num_paths)
        num_points = num_paths*points_per_path
        eigenvalues = np.zeros((num_bands, num_points))  # real
        count = 0
        for i in range(num_paths):
            for j in range(points_per_path):
                k = sym_points[i] + (sym_points[(i+1) % num_paths] - sym_points[i]) * float(j) / float(points_per_path - 1)
                k = np.matmul(k, bvec)
                eigvals = np.linalg.eigvals(hamiltonian(model, k, num_bands, avec, p_val=p, C_val=C, tx_factor_val=1))
                idx = np.argsort(eigvals)
                for band in range(num_bands):
                    eigenvalues[band, count] = np.real(eigvals[idx[band]])
                count += 1

        # construct figure
        fig = plt.figure()
        ax = plt.subplot(111)
        for i in range(num_bands):
            ax.plot(eigenvalues[i])
        for i in range(1, num_paths):
            ax.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
        ax.set_xlim([0, num_points])
        ax.set_xlabel('symmetry path')
        ax.set_ylabel('$E$')

    ####################################################################################################################

    # extract data from file
    with open('/home/bart/PycharmProjects/infinite_cylinder/code/standalone/FCI/bandstructure/gaps.txt', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        chern_numbers = []
        p_values = []
        widths = []
        gaps = []
        flatness_ratios = []
        for row in plots:
            chern_numbers.append(float(row[0]))
            p_values.append(float(row[1]))
            widths.append(float(row[2]))
            gaps.append(float(row[3]))
            flatness_ratios.append(float(row[3])/float(row[2]))

    ax1 = plt.subplot(right_inner_grid[0])
    ax2 = plt.subplot(right_inner_grid[1])
    ax3 = plt.subplot(right_inner_grid[2])

    for i in range(5):
        if i == 0:
            p_min = 0
            p_max = 6
        else:
            p_min = 6 + 7*(i-1)
            p_max = 6 + 7*(i-1) + 7
        ax1.plot(p_values[p_min:p_max], np.log(gaps[p_min:p_max]), '.', marker=markers[i], markersize=5, fillstyle='none',
                 c=f"C{i}", label=f"${int(chern_numbers[p_min])}$")
        ax2.plot(p_values[p_min:p_max], np.log(widths[p_min:p_max]), '.', marker=markers[i], markersize=5, fillstyle='none',
                 c=f"C{i}", label=f"${int(chern_numbers[p_min])}$")
        ax3.plot(p_values[p_min:p_max], np.log(flatness_ratios[p_min:p_max]), '.', marker=markers[i], markersize=5, fillstyle='none',
                 c=f"C{i}", label=f"${int(chern_numbers[p_min])}$")

    ax1.set_ylabel('$\ln(\Delta)$')
    ax2.set_ylabel('$\ln(W)$')
    ax3.set_ylabel('$\ln(\Delta / W)$')
    ax3.set_xlabel('$p$')
    ax3.set_xticks([3, 4, 5, 6, 7, 8, 9])

    leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=6, columnspacing=0.5,
                     bbox_to_anchor=(0.5, 1.8), title='$C$')
    leg.get_frame().set_linewidth(0.5)

    ####################################################################################################################

    ax.annotate('', xy=(0.17, 0.28), xycoords='axes fraction', xytext=(0.17, 0.3),
                arrowprops=dict(arrowstyle="<->", color='k'))
    ax.annotate('', xy=(0.82, 0.24), xycoords='axes fraction', xytext=(0.82, 0.27),
                arrowprops=dict(arrowstyle="<->", color='k'))
    fig.text(0.16, 0.31, "$\Delta$")
    fig.text(0.42, 0.3, "$W$")

    fig.text(0.25, 0.285, "$2$", c='w')
    fig.text(0.3, 0.335, "$-3$", c='w')
    fig.text(0.25, 0.4875, "$2$", c='w')
    fig.text(0.3, 0.64, "$-3$", c='w')
    fig.text(0.25, 0.684, "$2$", c='w')

    fig.text(0.13, 0.95, "(a)", fontsize=12)
    fig.text(0.525, 0.95, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/bandstructure_analysis_c2.png", bbox_inches='tight', dpi=300)
    plt.show()
