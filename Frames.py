import numpy as np

class Frames:

    def __init__(self, frame_size=10, n_frames=10, frames=None):
        if frames is None:
            self.frame_size = frame_size
            self.n_frames = n_frames
            self.frames = np.dstack([np.random.normal(size=(frame_size, frame_size)) for _ in range(n_frames)])
        else:
            self.frame_size = frame_size
            self.n_frames = n_frames
            self.frames = frames
            
    def edit_frame(self, frame_n, frame_edit_range, frame_edit_new_content):
        posx, posy, nx1, ny1 = frame_edit_range
        self.frames[:, :, frame_n][posx:posx + nx1, posy:posy + ny1] = frame_edit_new_content

    def get_frame(self, frame_n):
        return self.frames[:, :, frame_n]
    
    def get_multivariate_matrix(self):
        return self.frames.reshape((self.frames.shape[0] * self.frames.shape[1], self.frames.shape[2]))