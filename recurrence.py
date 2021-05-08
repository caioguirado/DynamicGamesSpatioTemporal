import numpy as np 
import scipy.io
from scipy.fftpack import fft
import scipy.signal as sgn
from tqdm import tqdm

from utils import ismember, get_neighbors, Graph

def spatio_emporal_detection_of_recurrence(signals,
                                           L,
                                           window_length=40,
                                           overlap=20,
                                           threshold_spectral_concentration=0.2,
                                           MinNumberofPointsInaRegion=100,
                                           min_prominence=0.6):
    step = window_length - overlap
    
    count_clusters = 0
    count_iteration = 0
    total_clusters = []
    interval_cluster = None
    numberofclusters = []
    
    signals_n_rows, signals_n_cols = signals.shape
    
    for index_window in tqdm(range(0, signals_n_cols, step)):

        
        if index_window + window_length <= signals_n_cols:
            B = signals[:, index_window : index_window + window_length]
        else:
            B = signals[:, index_window : signals_n_cols]
        
        B0 = B.copy()
        U, S, V = scipy.linalg.svd(B0, full_matrices=False)
        V = V.T # Make matrices compatible with matlab code
        A = U * S / np.sqrt(signals_n_rows-1)
        Z = np.sqrt(signals_n_rows-1) * V.T.conj()
        
        Yz = np.abs(fft(Z))
        Yz = Yz.T # Make matrices compatible with matlab code
        count_PC = 0
        pc_index = []
        for i in range(0, Yz.shape[1]):
            
            sy = Yz[:, i] / np.max(Yz[:, i])
        
            # This part outputs same code between matlab and python, only for sy it's 
            # mismatching (middle numbers are doubled)
            ltemp, _ = sgn.find_peaks(sy, prominence=min_prominence)
            
            # If no peaks, go to next element in loop
            if len(ltemp) == 0:
                continue
            
            prominences, left_bases, right_bases = sgn.peak_prominences(sy, ltemp)
            wtemp, _, _, _ = sgn.peak_widths(
                sy, 
                ltemp, 
                rel_height=0.5,
                prominence_data=(sy[ltemp], left_bases, right_bases)
            )
            indmaxp = int(np.max(sy[ltemp]))
            l = ltemp[indmaxp-1] + 1
            w = wtemp[indmaxp-1] 
            height_width_proportion = l - round(w/2)
            spectral_concentration = np.sum( Yz[int(l - round(w/2)) - 1 : int(l + round(w/2)), i]**2) / np.sum(Yz[:, i]**2) # review name
            if  height_width_proportion > 1 and spectral_concentration > threshold_spectral_concentration:
                count_PC += 1
                pc_index.append(i)
            elif height_width_proportion > 1 and spectral_concentration < threshold_spectral_concentration and spectral_concentration > 0:
                break
            else:
                pass
        
        
        J = count_PC
        positions = []
        for counter in range(J):
            pcs = np.abs(A[:, pc_index[counter]])
            max_pc = np.max(pcs)
            max_normalization = pcs / max_pc
            std_range = 3 * np.std(max_normalization)
            is_outlier_pcs = max_normalization > std_range
            positions += np.where(is_outlier_pcs)[0].tolist()
        
        temppos = list(set(positions))
        positions = np.array(sorted(temppos.copy()))
        
        diffs = (np.diff(positions) > 3).tolist()
        idx_greater_than_3_diff = [i for i, x in enumerate(diffs) if x]
        keep_indexes = [i for i, x in enumerate(positions) if i not in idx_greater_than_3_diff]
        positions = positions[keep_indexes].tolist()
        
        
        W = 246
        g = Graph(W * W) 
        
        clusters = []
        for pos in positions:
            neighbors = get_neighbors(matrix_shape=(W, W), i=pos//W, j=pos%W, complete=False)
            if len(set(neighbors).intersection(set(positions))) == 0:
                clusters.append([pos])
            else:
                for n in neighbors:
                    if n in positions:
                        g.addEdge(pos, n)
                    
        cc = g.connectedComponents()
        cc = [x for x in cc if len(x) > 1]
        clusters += cc
        
        positions4 = [c for c in clusters if len(c) >= MinNumberofPointsInaRegion]
        
        numberofclusters.append(len(positions4))     

        newclusters = []
        if numberofclusters[count_iteration] > 0 and count_clusters == 0:
            for i, cluster in enumerate(positions4):
                count_clusters += 1
                total_clusters.append(cluster)
                if interval_cluster is None:
                    interval_cluster = np.array([[index_window, 0]])
                else:
                    interval_cluster = np.vstack((interval_cluster, [index_window, 0]))
        elif numberofclusters[count_iteration] > 0 and count_clusters > 0:
            newclusters = positions4.copy()
            # check in the current window for clusters already identified in the previous window
            for pc_idx, prev_cluster in enumerate(total_clusters):
                flagcluster = 0
                for nc_idx, new_cluster in enumerate(newclusters):
                    # i1, i2 = ismember(prev_cluster, new_cluster)
                    intersection = set(prev_cluster).intersection(new_cluster)
                    # if sum(i1) > 0:
                    if len(intersection) > 0:
                        total_clusters[pc_idx] = total_clusters[pc_idx] + newclusters[nc_idx]
                        newclusters.pop(nc_idx)
                        flagcluster = 1
                        break

                if flagcluster == 0 and interval_cluster[pc_idx, 1] == 0:
                    interval_cluster[pc_idx, 1] = index_window
    
        if numberofclusters[count_iteration] == 0 and count_clusters > 0:
            interval_cluster[interval_cluster[:, 1] == 0, 1] = index_window
    
        # check for new clusters in the current window
        for nc_idx, new_cluster in enumerate(newclusters):
            if len(newclusters[nc_idx]) > 0:
                count_clusters += 1
                total_clusters.append(newclusters[nc_idx])
                interval_cluster = np.vstack((interval_cluster, [index_window, 0]))
                    
        count_iteration += 1
        
    ## End index_window loop
    
    if len(interval_cluster) > 0:
        interval_cluster[interval_cluster[:, 1] == 0, 1] = signals.shape[1]
    
    for index1 in range(len(total_clusters)):
        poscluster = []
        for index2 in range(index1+1, len(total_clusters)):
            i11, _ = ismember(total_clusters[index2], total_clusters[index1])
            if sum(i11) > 0:
                poscluster.append(index2)
        for index3 in range(len(poscluster)):
            total_clusters[index1] = list(set(total_clusters[index1] + total_clusters[poscluster[index3]]))
    
        if len(poscluster) != 0:
            interval_cluster[index1, :] = [interval_cluster[index1, 1], max(interval_cluster[[index1] + poscluster, 1])]
    
        for index3 in range(len(poscluster)):
            total_clusters[poscluster[index3]] = []
             
    tempTcluster = total_clusters
    tempinterval = interval_cluster
    
    countT = 0
    for index in range(len(tempTcluster)):
        if len(tempTcluster[index]) != 0:
            total_clusters[countT] = tempTcluster[index]
            interval_cluster[countT, :] = tempinterval[index, :]
            countT += 1

    return total_clusters, interval_cluster, numberofclusters



