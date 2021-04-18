import pickle 
import numpy as np 
import scipy.fft as fft
import scipy.signal as sgn

from utils import build_neighbors_matrix

with open('sinusoidal_simulation_frames.pkl', 'rb') as file:
    T = pickle.load(file)

t1, t2, t3 = T.frames.shape
B = T.frames.copy().reshape((t1 * t2, t3))

# Find connected points - identify regions
W = t1

L = build_neighbors_matrix(W, W)

def SpatioTemporalDetectionOFRecurrence(signals, 
                                        L,
                                        window_length=10,
                                        overlap=5,
                                        threshold_spectral_concentration=0.1,
                                        min_numberof_points_in_a_region=20,
                                        min_prominence=0.6):

    # Signals: matrix collecting the signal/time series from all points in the
    #          hyperplane under analysis (rows are points, and columns are samples)
    # L: sparse adjacency matrix collecting information about proximity among
    #    the points in the hyperplane (each row is a point, and the entries are
    #    the labels/positions of the points near it)

    # SpatioTemporalDetectionOFRecurrence(B, L, 40, 20, 0.2, 100)

    step = window_length - overlap

    count_clusters = 0
    count_iteration = 1
    total_clusters = []
    interval_cluster = []
    number_of_clusters = []

    signals_n_rows, signals_n_cols = signals.shape

    for index_window in range(0, signals_n_cols, step):

        if index_window + window_length <= signals_n_cols:
            B = signals[:, indexwindow : indexwindow + WindowLength]
        else:
            B = signals[:, indexwindow : signals_n_cols]

        B0 = B.copy()
        U, S, V = np.linalg.svd(B0)
        A = U * S / np.sqrt(signals_n_rows-1)
        Z = np.sqrt(signals_n_rows-1) * V.T.conj()

        Yz = np.abs(fft.fft(Z.T.conj()))
        count_PC = 0
        pc_index = []

        for i in range(0, Yz.shape[1]):
            sy = Yz[:, i] / np.max(Yz[:, i])

            # Review this part (how to translate findpeaks to python)
            ltemp, _ = sgn.find_peaks(a)
            wtemp = np.diff(locs) / 2

            indmaxp = np.max(sy[ltemp])
            l = ltemp[indmaxp]
            w = wtemp[indmaxp]
            if l - round(w/2) > 1 and np.sum( Yz[l - round(w/2) : l + round(w/2), i]**2) / np.sum(Yz[:, i]**2) > threshold_spectral_concentration:
                countPC += 1
                pc_index.append(i)
            elif l - round(w/2) > 1 and np.sum( Yz[l - round(w/2) : l + round(w/2), i]**2) / np.sum(Yz[:, i]**2) < threshold_spectral_concentration and np.sum( Yz[l - round(w/2) : l + round(w/2), i]**2) / np.sum(Yz[:, i]**2) > 0:
                break

        J = countPC
        positions = []
        for counter in range(J):
            pcs = np.abs(A[:, pc_index[counter]])
            max_pc = np.max(pcs)
            max_normalization = pcs / max_pc
            std_range = 3 * np.std(pcs)
            std_normalization = std_range / max_normalization
            is_outlier_pcs = max_normalization > std_normalization
            # is_outlier_pc = np.abs(A[:, pc_index(counter)]) / np.max(np.abs(A[:, pc_index(counter)])) > 3 * np.std(np.abs(A[:, pc_index(counter)]) / np.max(np.abs(A[:, pc_index(counter)])))
            positions += is_outlier_pcs

            temppos = list(set(positions))
            positions = temppos

            positions[np.diff(positions) > 3] = []