import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from Frames import Frames
from utils import create_wave_pattern, animate_frames, build_neighbors_matrix
from recurrence import spatio_emporal_detection_of_recurrence

# First wave parameters
W1 = 32
H1 = 32
freqRGB = 1./12
angleRGB = 2 * np.pi * 1./7
phaseRGB = 0

# Second wave parameters
W2 = 20 
H2 = 20
freqRGB2 = 1./5
angleRGB2 = 2 * np.pi * 1./3
phaseRGB2 = 0

# General parameters
W = 246 # Picture dimension (W x W)
n_frames = 150 # Number of frames

# Creating the two waves
imgr = create_wave_pattern(height=H1, width=W1, freq=freqRGB, angle=angleRGB, phase=phaseRGB)
imgr2 = create_wave_pattern(height=H2, width=W2, freq=freqRGB2, angle=angleRGB2,  phase=phaseRGB2)

# Unpacking images shapes to variables
nx, ny = (W, W)
nx1, ny1 = imgr.shape
nx2, ny2 = imgr2.shape

# Place at random location
posx = np.random.choice(range(nx1, nx-nx1+1))
posy = np.random.choice(range(ny1, ny-ny1+1))
posx2 = np.random.choice(range(nx2, nx-nx2+1))
posy2 = np.random.choice(range(nx2, ny-nx2+1))

# Creating frames object
T = Frames(frame_size=W, n_frames=n_frames)

for counter in range(40, 140):
    if counter < 80:
        # Only first wave alone
        phaseRGB += np.pi/180 * 20
        imgr = create_wave_pattern(height=H1, width=W1, freq=freqRGB, angle=angleRGB, phase=phaseRGB)
        T.edit_frame(frame_n=counter, frame_edit_range=(posx, posy, nx1, ny1), frame_edit_new_content=imgr)

    elif counter >= 80 and counter < 100:
        # Both first and second wave together
        phaseRGB += np.pi/180 * 20
        imgr = create_wave_pattern(height=H1, width=W1, freq=freqRGB, angle=angleRGB, phase=phaseRGB)
        T.edit_frame(frame_n=counter, frame_edit_range=(posx, posy, nx1, ny1), frame_edit_new_content=imgr)

        phaseRGB2 += np.pi/180 * 20
        imgr2 = create_wave_pattern(height=H2, width=W2, freq=freqRGB2, angle=angleRGB2, phase=phaseRGB2)
        T.edit_frame(frame_n=counter, frame_edit_range=(posx2, posy2, nx2, ny2), frame_edit_new_content=imgr2)

    elif counter >= 100 and counter <= 140:
        # Only second wave alone
        phaseRGB2 += np.pi/180 * 20
        imgr2 = create_wave_pattern(height=H2, width=W2, freq=freqRGB2, angle=angleRGB2, phase=phaseRGB2)
        T.edit_frame(frame_n=counter, frame_edit_range=(posx2, posy2, nx2, ny2), frame_edit_new_content=imgr2)


if __name__ == "__main__":
    
    with open('sinusoidal_simulation_frames.pkl', 'wb') as file:
        pickle.dump(T, file)

    B = T.get_multivariate_matrix()
    
    ani = animate_frames(T)

    L = build_neighbors_matrix(W, W)
    
    total_clusters, interval_cluster, numberofclusters = spatio_emporal_detection_of_recurrence(
                                                                   signals=B,
                                                                   L=L,
                                                                   window_length=40,
                                                                   overlap=20,
                                                                   threshold_spectral_concentration=0.2,
                                                                   MinNumberofPointsInaRegion=100,
                                                                   min_prominence=0.6
                                                        )
    
    for cluster in total_clusters:
        img = np.zeros((W * W))
        img[cluster] = 1
        img = img.reshape((W, W))
        plt.imshow(img)

    