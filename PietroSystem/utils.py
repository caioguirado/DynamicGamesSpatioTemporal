import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def create_wave_pattern(height, width, freq, angle, phase=0):
    X, Y = np.meshgrid(range(width), range(height))
    wave_pattern = np.sin(2 * np.pi * freq * (X * np.cos(angle) + Y * np.sin(angle)) + phase)
    return wave_pattern

def animate_frames(frames_obj):
    fig, ax = plt.subplots()
    ims = []
    for i in range(frames_obj.frames.shape[-1]):
        im = ax.imshow(frames_obj.get_frame(frame_n=i), animated=True, cmap='Greys',  interpolation='nearest')
        ims.append([im])

    ani = animation.ArtistAnimation(fig, 
                                    ims, 
                                    interval=50, 
                                    repeat_delay=1000)
    plt.show()

def get_neighbors(matrix_shape, i, j):
    total_neighbors = 8
    n_rows, n_cols = matrix_shape
    neighbors = []
    for row in range(i-1, i+2):
        for col in range(j-1, j+2):
            if row != i or col != j:
                if row >= 0 and row < n_rows and col >= 0 and col < n_cols:
                    element_order_1d = row*n_cols + col
                    neighbors.append(element_order_1d)

    # Complete with zeros to keep same format as L
    neighbors_size = len(neighbors)
    if neighbors_size != total_neighbors:
        for _ in range(total_neighbors - neighbors_size):
            neighbors.append(0)

    return neighbors

def build_neighbors_matrix(n_rows, n_cols):
    all_neighbors = []
    for i in range(n_rows):
        for j in range(n_cols):
            neighbors = get_neighbors((n_rows, n_cols), i, j)
            all_neighbors.append(neighbors)
    L = np.array(all_neighbors)
    return L


# L = np.zeros((W**2, 8)) 

# L[0, :] = [2, 1+W, 2+W, 0, 0, 0, 0, 0]
# L[W-1, :] = [W-1, 2*W, 2*W-1, 0, 0, 0, 0, 0]
# L[W**2-W+1-1, :] = [W**2-W+2, W**2-2*W+1, W**2-2*W+2, 0, 0, 0, 0, 0]
# L[W**2-1, :] = [W**2-1, W**2-W, W**2-W-1, 0, 0, 0, 0, 0]
    
# for ind in range(1, W):
#     L[ind, :] = [ind, ind+2, ind+W, ind+W+1, ind+W+2, 0, 0, 0]

# for ind in range(W**2-W+1, W**2):
#     L[ind, :] = [ind, ind+2, ind-W, ind-W+1, ind-W+2, 0, 0, 0]

# for ind in range(W, W**2-W+1):
#     if (ind+1) % W != 0 and (ind+1) % W != 1:
#         L[ind, :] = [ind, ind+2, ind-W, ind-W+1, ind-W+2, ind+W, ind+W+1, ind+W+2]
#     elif (ind+1) % W == 0:
#         L[ind, :] = [ind-W, ind-W+1, ind, ind+W, ind+W+1, 0, 0, 0]
#     elif (ind+1) % W == 1:
#         L[ind, :] = [ind-W+1, ind-W+2, ind+2, ind+W+1, ind+W+2, 0, 0, 0]