import numpy as npfrom tqdm import tqdmimport randomimport mathimport itertoolsimport copyclass EvolutionaryGame:        # TODO: check iniFrame        def __init__(self,                  initial_config: np.ndarray,                  size: tuple,                  n_frames: int,                  position: tuple,                 fitness_matrix: np.ndarray,                 death_rate: float,                 inter_radius: int,                 offspring_radius: int                 ) -> None:                self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()        self.fitness_matrix = fitness_matrix        self.death_rate = death_rate        self.inter_radius = inter_radius        self.offspring_radius = offspring_radius                # colors        self.col_about_to_die = (200, 200, 225)        self.col_type1 = (255, 255, 0)        self.col_type2 = (255, 0, 255)        self.col_type3 = (67, 216, 255)        self.col_background = (10, 10, 40)        self.col_grid = (30, 30, 60)            def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0],               self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config                return cells        def generate_simulation(self) -> None:                for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame        frames_3d = []        for i in range(self.simulation.shape[-1]):            frame = self.simulation[:, :, i]            frame_3d = np.zeros((*self.grid.shape, 3))            mask1 = frame == 1            mask2 = frame == 2            mask3 = frame == 3            for color_channel in range(3): # 3 color channels                frame_3d[:, :, color_channel][mask1] = self.col_type1[color_channel]                frame_3d[:, :, color_channel][mask2] = self.col_type2[color_channel]                frame_3d[:, :, color_channel][mask3] = self.col_type3[color_channel]                            frames_3d.append(frame_3d)                    self.simulation_3d = np.array(frames_3d)    def step(self) -> np.ndarray:                nxt = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        for r, c in np.ndindex(self.grid.shape):            if nxt[r,c] == 0:                vocal = self.grid[r,c]                nxt[r,c] = vocal                if vocal != 0:                    # lets kill some of the old generation                    death = np.random.uniform(0.0, 1.0, 1)                        if death <= self.death_rate:                        nxt[r,c]= 0                        # pick a random individual from inetarction_radius to mate                    mate_pos, mate_value = self.radius_rand(r, c, 1)                    if vocal > 3 or vocal < 0:                        raise Exception("The type of the player is unknown")                        # Placing the offspring in offspring_radius according to A = [[1,2],[2,1]]                    if mate_value is not None:                            nr_offspring = self.fitness_matrix[int(vocal-1), int(mate_value-1)]                        for i in range(nr_offspring):                            offspring_pos, offspring_value = self.radius_rand(r, c, 0)                            if offspring_value is not None:                                nxt[offspring_pos[0], offspring_pos[1]] = vocal        return nxt        def radius_rand(self, r, c, value):        chosen_pos = None        chosen_value = None        R = np.arange(0, self.grid.shape[0], 1, dtype=int)        C = np.arange(0, self.grid.shape[1], 1, dtype=int)                # defining the begin and end of the neighberhood        beginR = r - self.offspring_radius        endR = r + self.offspring_radius         beginC = c - self.offspring_radius        endC = c + self.offspring_radius            for i in range(8 * self.offspring_radius):            rand_r = random.randint(beginR, endR)            rand_c = random.randint(beginC, endC)                        rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]            rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]                        if value == 0:                if self.grid[rand_r, rand_c] == 0:                    chosen_pos = [rand_r, rand_c]                    chosen_value = 0                    return chosen_pos, chosen_value            else:                 if self.grid[rand_r, rand_c] != 0:                    chosen_pos = [rand_r, rand_c]                    chosen_value = self.grid[rand_r, rand_c]                    return chosen_pos, chosen_value                        return chosen_pos, chosen_valueclass EvolutionaryGame2:    # TODO: check iniFrame    def __init__(self,                 initial_config: np.ndarray,                 size: tuple,                 n_frames: int,                 position: tuple,                 fitness_matrix: np.ndarray,                 death_rate: float,                 inter_radius: int,                 offspring_radius: int                 ) -> None:        self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()        self.fitness_matrix = fitness_matrix        self.death_rate = death_rate        self.inter_radius = inter_radius        self.offspring_radius = offspring_radius        # colors        self.col_about_to_die = (200, 200, 225)        self.col_type1 = (255, 255, 0)        self.col_type2 = (255, 0, 255)        self.col_type3 = (67, 216, 255)        self.col_background = (10, 10, 40)        self.col_grid = (30, 30, 60)    def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0],        self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config        return cells    def generate_simulation(self) -> None:        for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame        frames_3d = []        for i in range(self.simulation.shape[-1]):            frame = self.simulation[:, :, i]            frame_3d = np.zeros((*self.grid.shape, 3))            mask1 = frame == 1            mask2 = frame == 2            mask3 = frame == 3            for color_channel in range(3):  # 3 color channels                frame_3d[:, :, color_channel][mask1] = self.col_type1[color_channel]                frame_3d[:, :, color_channel][mask2] = self.col_type2[color_channel]                frame_3d[:, :, color_channel][mask3] = self.col_type3[color_channel]            frames_3d.append(frame_3d)        self.simulation_3d = np.array(frames_3d)    def step(self) -> np.ndarray:        nxt = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        nxt2 = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        for r, c in np.ndindex(self.grid.shape):            # if nxt[r, c] == 0:                vocal = self.grid[r, c]                # nxt[r, c] = vocal                if vocal != 0:                    # lets kill some of the old generation                    death = np.random.uniform(0.0, 1.0, 1)                    if death <= self.death_rate:                        nxt[r, c] = 0                    # pick a random individual from inetarction_radius to mate                    mate_pos, mate_value = self.radius_rand(r, c, 1)                    if vocal > 3 or vocal < 0:                        raise Exception("The type of the player is unknown")                    # Placing the offspring in offspring_radius according to A = [[1,2],[2,1]]                    if mate_value is not None:                        nr_offspring = self.fitness_matrix[int(vocal - 1), int(mate_value - 1)]                        for i in range(nr_offspring):                            offspring_pos, offspring_value = self.radius_rand(r, c, 0, mate_value)                            if offspring_value is not None:                                nxt[offspring_pos[0], offspring_pos[1]] = vocal        for r, c in np.ndindex(self.grid.shape):            nxt2[r, c] = nxt[r, c] if nxt[r, c] > 0 else self.grid[r, c]        return nxt2    def radius_rand(self, r, c, value, mate_value = 0):        chosen_pos = None        chosen_value = None        R = np.arange(0, self.grid.shape[0], 1, dtype=int)        C = np.arange(0, self.grid.shape[1], 1, dtype=int)        # defining the begin and end of the neighberhood        beginR = r - self.offspring_radius        endR = r + self.offspring_radius        beginC = c - self.offspring_radius        endC = c + self.offspring_radius        for i in range(8 * self.offspring_radius):            rand_r = random.randint(beginR, endR)            rand_c = random.randint(beginC, endC)            rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]            rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]            if value == 0:                if self.grid[rand_r, rand_c] == 0 or self.grid[rand_r, rand_c]  == mate_value:                    chosen_pos = [rand_r, rand_c]                    chosen_value = 0                    return chosen_pos, chosen_value            else:                if self.grid[rand_r, rand_c] != 0:                    chosen_pos = [rand_r, rand_c]                    chosen_value = self.grid[rand_r, rand_c]                    return chosen_pos, chosen_value        return chosen_pos, chosen_valueclass EvolutionaryGame3:    # TODO: check iniFrame    def __init__(self,                 initial_config: np.ndarray,                 size: tuple,                 n_frames: int,                 position: tuple,                 fitness_matrix: np.ndarray,                 death_rate: float,                 inter_radius: int,                 offspring_radius: int                 ) -> None:        self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()        self.fitness_matrix = fitness_matrix        self.death_rate = death_rate        self.inter_radius = inter_radius        self.offspring_radius = offspring_radius        # colors        self.col_about_to_die = (200, 200, 225)        self.col_type1 = (255, 255, 0)        self.col_type2 = (255, 0, 255)        self.col_type3 = (67, 216, 255)        self.col_background = (10, 10, 40)        self.col_grid = (30, 30, 60)    def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0],        self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config        return cells    def generate_simulation(self) -> None:        for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame        frames_3d = []        for i in range(self.simulation.shape[-1]):            frame = self.simulation[:, :, i]            frame_3d = np.zeros((*self.grid.shape, 3))            mask1 = frame == 1            mask2 = frame == 2            mask3 = frame == 3            for color_channel in range(3):  # 3 color channels                frame_3d[:, :, color_channel][mask1] = self.col_type1[color_channel]                frame_3d[:, :, color_channel][mask2] = self.col_type2[color_channel]                frame_3d[:, :, color_channel][mask3] = self.col_type3[color_channel]            frames_3d.append(frame_3d)        self.simulation_3d = np.array(frames_3d)    def step(self) -> np.ndarray:        nxt = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        cells = []        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        for r, c in np.ndindex(self.grid.shape):            vocal = int(self.grid[r, c])            if vocal != 0:                # NO DEATH RATE                # pick a random individual from interaction_radius to mate                mate_pos, mate_value = self.radius_rand(r, c, 1)                if vocal > 3 or vocal < 0:                    raise Exception("The type of the player is unknown")                # Placing the offspring in offspring_radius according to A = [[1,2],[2,1]]                if mate_value is not None:                    mate_value = int(mate_value)                    nr_offspring = self.fitness_matrix[int(vocal - 1), int(mate_value - 1)]                    negative = True if nr_offspring < 0 else False                    for i in range(abs(nr_offspring)):                        offspring_pos, offspring_value = self.radius_rand(r, c, 0, mate_value)                        if offspring_value is not None:                            if negative:                                cells[mate_value-1][offspring_pos[0], offspring_pos[1]] = mate_value                            else:                                cells[vocal - 1][offspring_pos[0], offspring_pos[1]] = vocal        for r, c in np.ndindex(self.grid.shape):            randomList = list(range(3))            random.shuffle(randomList)            nxt[r, c] = self.grid[r, c]            for choice in randomList:                if cells[choice][r, c] != 0:                    nxt[r, c] = cells[choice][r, c]                    break;        return nxt    def radius_rand(self, r, c, value, mate_value = 0):        chosen_pos = None        chosen_value = None        R = np.arange(0, self.grid.shape[0], 1, dtype=int)        C = np.arange(0, self.grid.shape[1], 1, dtype=int)        # defining the begin and end of the neighberhood        beginR = r - self.offspring_radius        endR = r + self.offspring_radius        beginC = c - self.offspring_radius        endC = c + self.offspring_radius        for i in range(8 * self.offspring_radius):            rand_r = random.randint(beginR, endR)            rand_c = random.randint(beginC, endC)            rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]            rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]            if value == 0:                if self.grid[rand_r, rand_c] == 0 or self.grid[rand_r, rand_c]  == mate_value:                    chosen_pos = [rand_r, rand_c]                    chosen_value = 0                    return chosen_pos, chosen_value            else:                if self.grid[rand_r, rand_c] != 0:                    chosen_pos = [rand_r, rand_c]                    chosen_value = self.grid[rand_r, rand_c]                    return chosen_pos, chosen_value        return chosen_pos, chosen_valueclass EvolutionaryGame4:    # TODO: check iniFrame    def __init__(self,                 initial_config: np.ndarray,                 size: tuple,                 n_frames: int,                 position: tuple,                 fitness_matrix: np.ndarray,                 death_rate: float,                 inter_radius: int,                 offspring_radius: int,                 d_when: tuple,                 d_where: tuple,                 d_pattern: np.ndarray                 ) -> None:        self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()        self.fitness_matrix = fitness_matrix        self.death_rate = death_rate        self.inter_radius = inter_radius        self.offspring_radius = offspring_radius        self.d_when = d_when        self.d_where = d_where        self.d_pattern = d_pattern    def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0], self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config        return cells    def generate_simulation(self) -> None:        for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            if i in self.d_when:                next_frame[self.d_where[0]:self.d_where[0] + self.d_pattern.shape[0], self.d_where[1]:self.d_where[1] + self.d_pattern.shape[1]] = self.d_pattern            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame    def step(self) -> np.ndarray:        nxt = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        cells = []        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        for r, c in np.ndindex(self.grid.shape):            vocal = int(self.grid[r, c])            if vocal != 0:                # Create a list with all the neighbors relative position                list_of_pos_inter = list(itertools.product(np.arange(-self.inter_radius, self.inter_radius + 1, 1), repeat=2))                list_of_pos_inter.remove((0, 0))                # NO DEATH RATE                # pick a random individual from interaction_radius to mate                mate_pos, mate_value, _, mate_rel = self.radius_rand(r, c, 1, list_of_pos_inter)                if vocal > 3 or vocal < 0:                    raise Exception("The type of the player is unknown")                # Placing the offspring in offspring_radius according to A = [[1,2],[2,1]]                if mate_value is not None:                    # Create a list with all the neighbors relative position                    list_of_pos_off=list(itertools.product(np.arange(-self.offspring_radius, self.offspring_radius + 1, 1),repeat=2))                    list_of_pos_off.remove(mate_rel)  ## Take a look here                    # list_of_pos_off.remove((0, 0)) ## Take a look here                    mate_value = int(mate_value)                    nr_offspring = self.fitness_matrix[int(vocal - 1), int(mate_value - 1)]                    negative = True if nr_offspring < 0 else False                    for i in range(abs(nr_offspring)):                        offspring_pos, offspring_value, list_of_pos_off, _ = self.radius_rand(r, c, 0, list_of_pos_off, mate_value)                        if offspring_value is not None:                            if negative:                                cells[mate_value-1][offspring_pos[0], offspring_pos[1]] = mate_value                            else:                                cells[vocal - 1][offspring_pos[0], offspring_pos[1]] = vocal                else:                    cells[vocal - 1][mate_pos[0], mate_pos[1]] = vocal        for r, c in np.ndindex(self.grid.shape):            randomList = list(range(3))            random.shuffle(randomList)            nxt[r, c] = self.grid[r, c]            for choice in randomList:                if cells[choice][r, c] != 0:                    nxt[r, c] = cells[choice][r, c]                    break;        return nxt    def radius_rand(self, r, c, value, list_of_pos, mate_value = 0):        chosen_pos = None        chosen_value = None        mate_rel = None        R = np.arange(0, self.grid.shape[0], 1, dtype=int)        C = np.arange(0, self.grid.shape[1], 1, dtype=int)        if list_of_pos:            choice = random.choice(list_of_pos)            rand_r = r + choice[0]            rand_c = c + choice[1]            list_of_pos.remove(choice)            rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]            rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]            chosen_pos = [rand_r, rand_c]            if value == 0:                if self.grid[rand_r, rand_c] == 0 or self.grid[rand_r, rand_c]  == mate_value:                    chosen_value = 0            else:                if self.grid[rand_r, rand_c] != 0:                    chosen_value = self.grid[rand_r, rand_c]                    mate_rel = (-choice[0],-choice[1])        return chosen_pos, chosen_value, list_of_pos, mate_relclass EvolutionaryGame5:    # TODO: check iniFrame    def __init__(self,                 initial_config: np.ndarray,                 size: tuple,                 n_frames: int,                 position: tuple,                 fitness_matrix: np.ndarray,                 death_rate: float,                 inter_radius: int,                 offspring_radius: int,                 d_when: tuple,                 d_where: tuple,                 d_pattern: np.ndarray                 ) -> None:        self.initial_config = initial_config        self.size = size        self.n_frames = n_frames        self.position = position        self.grid = self._init_grid()        self.simulation = self.grid.copy()        self.fitness_matrix = fitness_matrix        self.death_rate = death_rate        self.inter_radius = inter_radius        self.offspring_radius = offspring_radius        self.d_when = d_when        self.d_where = d_where        self.d_pattern = d_pattern    def _init_grid(self) -> np.ndarray:        cells = np.zeros(self.size)        cells[self.position[0]:self.position[0] + self.initial_config.shape[0], self.position[1]:self.position[1] + self.initial_config.shape[1]] = self.initial_config        return cells    def generate_simulation(self) -> None:        for i in tqdm(range(1, self.n_frames)):            next_frame = self.step()            if i in self.d_when:                next_frame[self.d_where[0]:self.d_where[0] + self.d_pattern.shape[0], self.d_where[1]:self.d_where[1] + self.d_pattern.shape[1]] = self.d_pattern            self.simulation = np.dstack((self.simulation, next_frame))            self.grid = next_frame    def step(self) -> np.ndarray:        nxt = np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8)        cells = []        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        cells.append(np.zeros((self.grid.shape[0], self.grid.shape[1]), dtype=np.int8))        for r, c in np.ndindex(self.grid.shape):            vocal = int(self.grid[r, c])            if vocal != 0:                # Create a list with all the neighbors relative position                list_of_pos_inter = list(itertools.product(np.arange(-self.inter_radius, self.inter_radius + 1, 1), repeat=2))                list_of_pos_inter.remove((0, 0))                # NO DEATH RATE                # pick a random individual from interaction_radius to mate                mate_pos, mate_value, _, mate_rel = self.radius_rand(r, c, 1, list_of_pos_inter)                if vocal > 3 or vocal < 0:                    raise Exception("The type of the player is unknown")                # Placing the offspring in offspring_radius according to A = [[1,2],[2,1]]                if mate_value is not None:                    # Create a list with all the neighbors relative position                    list_of_pos_off=list(itertools.product(np.arange(-self.offspring_radius, self.offspring_radius + 1, 1),repeat=2))                    if self.offspring_radius >= self.inter_radius:                        list_of_pos_off.remove(mate_rel)  ## Take a look here                    mate_value = int(mate_value)                    nr_offspring = self.fitness_matrix[int(vocal - 1), int(mate_value - 1)]                    negative = True if nr_offspring < 0 else False                    for i in range(abs(nr_offspring)):                        offspring_pos, offspring_value, list_of_pos_off, _ = self.radius_rand(mate_pos[0], mate_pos[1], 0, list_of_pos_off, mate_value)                        if offspring_value is not None:                            if negative:                                cells[mate_value-1][offspring_pos[0], offspring_pos[1]] = mate_value                            else:                                cells[vocal - 1][offspring_pos[0], offspring_pos[1]] = vocal                else:                    cells[vocal - 1][mate_pos[0], mate_pos[1]] = vocal        for r, c in np.ndindex(self.grid.shape):            randomList = list(range(3))            random.shuffle(randomList)            nxt[r, c] = self.grid[r, c]            for choice in randomList:                if cells[choice][r, c] != 0:                    nxt[r, c] = cells[choice][r, c]                    break;        return nxt    def radius_rand(self, r, c, value, list_of_pos, mate_value = 0):        vocal = self.grid[r, c]        chosen_pos = None        chosen_value = None        mate_rel = None        R = np.arange(0, self.grid.shape[0], 1, dtype=int)        C = np.arange(0, self.grid.shape[1], 1, dtype=int)        if list_of_pos:            choice = random.choice(list_of_pos)            rand_r = r + choice[0]            rand_c = c + choice[1]            list_of_pos.remove(choice)            rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]            rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]            chosen_pos = [rand_r, rand_c]            if value == 0:                if self.grid[rand_r, rand_c] == 0 or self.grid[rand_r, rand_c]  == mate_value:                    chosen_value = 0            else:                for a in list_of_pos:                    choice = random.choice(list_of_pos)                    rand_r = r + choice[0]                    rand_c = c + choice[1]                    list_of_pos.remove(choice)                    rand_r = R[rand_r] if rand_r < self.grid.shape[0] else R[rand_r - self.grid.shape[0]]                    rand_c = C[rand_c] if rand_c < self.grid.shape[1] else R[rand_c - self.grid.shape[1]]                    if self.grid[rand_r, rand_c] != 0 and self.grid[rand_r, rand_c] != vocal:                        chosen_value = self.grid[rand_r, rand_c]                        mate_rel = (-choice[0],-choice[1])                        break        return chosen_pos, chosen_value, list_of_pos, mate_rel# =============================================================================# initial = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0],#                     [1,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [2,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#                     [0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]], dtype=np.int8)# SIZE = (100, 100)# N_FRAMES = 50# POSITION = (5,10)# FITNESS_MATRIX = np.array([[1,8,1],#                            [1,1,8],#                            [8,1,1]])# DEATH_RATE = 0.1# INTER_RADIUS = 1# OFFSPRINT_RADIUS = 1# evo_game = EvolutionaryGame(initial_config=initial, #                             size=SIZE,#                             n_frames=N_FRAMES,#                             position=POSITION,#                             fitness_matrix=FITNESS_MATRIX,#                             death_rate=DEATH_RATE,#                             inter_radius=INTER_RADIUS,#                             offspring_radius=OFFSPRINT_RADIUS)# # evo_game.generate_simulation()# =============================================================================