import numpy as npfrom tqdm import tqdmclass GameOfLife:        # TODO: check iniFrame        def __init__(self, initial_config: np.ndarray, size: tuple, n_frames: int, position: tuple) -> None:        self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()            def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0],               self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config                return cells        def generate_simulation(self) -> None:                for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame                            def step(self):                nxt = np.zeros(self.grid.shape)        for r, c in np.ndindex(self.grid.shape):            num_alive = np.sum(self.grid[r-1:r+2, c-1:c+2]) - self.grid[r, c]                if (self.grid[r, c] == 1 and 2 <= num_alive <= 3) or (self.grid[r, c] == 0 and num_alive == 3):                nxt[r, c] = 1            return nxt                # =============================================================================# initial_pattern = np.array([ [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],#                                 [1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                                 [1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],#                                 [0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]])# # gol = GameOfLife(initial_pattern, (100, 100), 50, (3,3))# gol.generate_simulation()# # from utils import animate_frames# from Frames import Frames# # b = Frames(frames=gol.simulation)# # animate_frames(b)# # np.savez_compressed('initial_config', initial_pattern)# # a = np.load('output_compresed.npz')# a['arr_0']# =============================================================================