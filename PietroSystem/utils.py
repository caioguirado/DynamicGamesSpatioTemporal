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

def get_neighbors(matrix_shape, i, j, complete=True):
    total_neighbors = 8
    n_rows, n_cols = matrix_shape
    neighbors = []
    for row in range(i-1, i+2):
        for col in range(j-1, j+2):
            if row != i or col != j:
                if row >= 0 and row < n_rows and col >= 0 and col < n_cols:
                    element_order_1d = row*n_cols + col
                    neighbors.append(element_order_1d)

    if complete:
        # Complete with zeros to keep same format as L
        neighbors_size = len(neighbors)
        if neighbors_size != total_neighbors:
            for _ in range(total_neighbors - neighbors_size):
                neighbors.append(0)
    else:
        return neighbors

    return neighbors

def build_neighbors_matrix(n_rows, n_cols):
    all_neighbors = []
    for i in range(n_rows):
        for j in range(n_cols):
            neighbors = get_neighbors((n_rows, n_cols), i, j)
            all_neighbors.append(neighbors)
    L = np.array(all_neighbors)
    return L

def ismember(a, b):
    if type(a) == int:
        a = [a]
    position_list = []
    # is_member = [1 if x in b else 0 for x in a]
    is_member = np.in1d(a, b).tolist()
    for i, x in enumerate(a):
        try:
            index = b.index(x)
        except Exception as e:
            index = -1 # review            
        position_list.append(index)
    
    return is_member, position_list

class Graph:
 
    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]
 
    def DFSUtil(self, temp, v, visited):
 
        # Mark the current vertex as visited
        visited[v] = True
 
        # Store the vertex to list
        temp.append(v)
 
        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
 
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp
 
    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)
 
    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc
    
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