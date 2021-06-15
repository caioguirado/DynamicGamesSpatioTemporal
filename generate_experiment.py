#!/usr/bin/pythonimport sysimport getoptfrom enum import Enumfrom game_of_life import GameOfLifefrom evolutionary_game import EvolutionaryGame, EvolutionaryGame2, EvolutionaryGame3, EvolutionaryGame4, EvolutionaryGame5from utils import create_wave_pattern, animate_frames, build_neighbors_matrixfrom Frames import Framesimport os.pathfrom recurrence import spatio_emporal_detection_of_recurrenceimport numpy as npimport matplotlib.pyplot as pltimport matplotlib.animation as animationfrom PIL import Imageimport PILimport resourceresource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))sys.setrecursionlimit(10**6)class GameType(Enum):    GAME_OF_LIFE=1    EVO_GAME=2    EVO_GAME2=3    EVO_GAME3=4    EVO_GAME4=5    EVO_GAME5=6    def run_recurrence(frames: Frames, exp_path: str) -> None:    B = frames.get_multivariate_matrix()    L = build_neighbors_matrix(*frames.frames.shape[:2])        variables = {}    exec(open(exp_path + '/config.py').read(), variables)        window_length = variables['WINDOW_LENGTH']    overlap = variables['OVERLAP']    threshold_spectral_concentration = variables['THRESHOLD_SPECTRAL_CONCENTRATION']    MinNumberofPointsInaRegion = variables['MinNumberofPointsInaRegion']    min_prominence = variables['MIN_PROMINENCE']        print('============ Running Pattern Recognition ============')    total_clusters, interval_cluster, numberofclusters = spatio_emporal_detection_of_recurrence(                                                                   signals=B,                                                                   L=L,                                                                   window_length=window_length,                                                                   overlap=overlap,                                                                   threshold_spectral_concentration=threshold_spectral_concentration,                                                                   MinNumberofPointsInaRegion=MinNumberofPointsInaRegion,                                                                   min_prominence=min_prominence                                                        )            patterns = []    for cluster in total_clusters:        w1, w2, w3 = frames.frames.shape        img = np.zeros((w1*w2))        img[cluster] = 1        img = img.reshape((w1, w2))        patterns.append(img)            combined_patterns = np.dstack(patterns)    combined_patterns_frames = Frames(frames=combined_patterns)    ani = animate_frames(combined_patterns_frames,                          show=False,                          titles=interval_cluster)        print('============ Saving Experiment Animation ============')    ani.save(exp_path + '/experiment.gif', writer='Pillow', fps=60)    ani.save(exp_path + '/experiment.mp4', writer=animation.FFMpegWriter(fps=30))    combined_clusters = np.any(combined_patterns, axis=2)    resize = 5    Image.fromarray(combined_clusters).convert('RGB').resize((combined_clusters.shape[0]*resize, combined_clusters.shape[1]*resize)).save(exp_path + '/patterns1.png')    plt.imsave(exp_path + '/patterns.png', combined_clusters)    def generate_experiment(argv):    try:       opts, args = getopt.getopt(argv,                                   shortopts="",                                   longopts=['game_type=',                                             'pattern_folder_name=',                                            'exp_folder_name='])           except getopt.GetoptError as e:       print('Error: ', e)       sys.exit(2)           for opt, arg in opts:        if opt not in ('--game_type', '--exp_folder_name', '--pattern_folder_name'):            print('Invalid parameter')            sys.exit()        elif opt == '--game_type':            game_type = arg        elif opt == '--exp_folder_name':            exp_folder = arg        elif opt == '--pattern_folder_name':            pattern_folder = arg        if GameType[game_type].value == 1:        # game of life        base_path = './simulations/game_of_life/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder                if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:            variables = {}             exec(open(base_path + '/config.py').read(), variables)                        initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            # generate output files to exp_path            gol = GameOfLife(initial_config=initial_config,                              size=size,                              n_frames=n_frames,                              position=position)                        print('============ Generating Game Simulation ============')            gol.generate_simulation()            frames = Frames(frames=gol.simulation)            ani = animate_frames(frames, show=False)                        print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)                run_recurrence(frames=frames, exp_path=exp_path)                    elif GameType[game_type].value == 2:         # evo game        base_path = './simulations/evo_game/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder                if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:                variables = {}             exec(open(base_path + '/config.py').read(), variables)                        initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            fitness_matrix = variables['FITNESS_MATRIX']            death_rate = variables['DEATH_RATE']            inter_radius = variables['INTER_RADIUS']            offspring_radius = variables['OFFSPRINT_RADIUS']                        eg = EvolutionaryGame(initial_config=initial_config,                                   size=size,                                   n_frames=n_frames,                                   position=position,                                   fitness_matrix=fitness_matrix,                                   death_rate=death_rate,                                   inter_radius=inter_radius,                                   offspring_radius=offspring_radius)                        print('============ Generating Game Simulation ============')            eg.generate_simulation()            frames = Frames(frames=eg.simulation)            ani = animate_frames(frames, show=False, cmap='plasma')                        print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)                    run_recurrence(frames=frames, exp_path=exp_path)    elif GameType[game_type].value == 3:        # evo game 2        base_path = './simulations/evo_game2/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder        if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:            variables = {}            exec(open(base_path + '/config.py').read(), variables)            initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            fitness_matrix = variables['FITNESS_MATRIX']            death_rate = variables['DEATH_RATE']            inter_radius = variables['INTER_RADIUS']            offspring_radius = variables['OFFSPRINT_RADIUS']            eg = EvolutionaryGame2(initial_config=initial_config,                                  size=size,                                  n_frames=n_frames,                                  position=position,                                  fitness_matrix=fitness_matrix,                                  death_rate=death_rate,                                  inter_radius=inter_radius,                                  offspring_radius=offspring_radius)            print('============ Generating Game Simulation ============')            eg.generate_simulation()            frames = Frames(frames=eg.simulation)            ani = animate_frames(frames, show=False, cmap='plasma')            print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)        run_recurrence(frames=frames, exp_path=exp_path)    elif GameType[game_type].value == 4:        # evo game 3        base_path = './simulations/evo_game3/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder        if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:            variables = {}            exec(open(base_path + '/config.py').read(), variables)            initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            fitness_matrix = variables['FITNESS_MATRIX']            death_rate = variables['DEATH_RATE']            inter_radius = variables['INTER_RADIUS']            offspring_radius = variables['OFFSPRINT_RADIUS']            eg = EvolutionaryGame3(initial_config=initial_config,                                  size=size,                                  n_frames=n_frames,                                  position=position,                                  fitness_matrix=fitness_matrix,                                  death_rate=death_rate,                                  inter_radius=inter_radius,                                  offspring_radius=offspring_radius)            print('============ Generating Game Simulation ============')            eg.generate_simulation()            frames = Frames(frames=eg.simulation)            ani = animate_frames(frames, show=False, cmap='RdYlBu')            print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)        run_recurrence(frames=frames, exp_path=exp_path)    elif GameType[game_type].value == 5:        # evo game 4        base_path = './simulations/evo_game4/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder        if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:            variables = {}            exec(open(base_path + '/config.py').read(), variables)            initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            fitness_matrix = variables['FITNESS_MATRIX']            death_rate = variables['DEATH_RATE']            inter_radius = variables['INTER_RADIUS']            offspring_radius = variables['OFFSPRINT_RADIUS']            d_when = variables['DISASTER_WHEN']            d_where = variables['DISASTER_WHERE']            d_pattern = variables['DISASTER_PATTERN']            eg = EvolutionaryGame4(initial_config=initial_config,                                  size=size,                                  n_frames=n_frames,                                  position=position,                                  fitness_matrix=fitness_matrix,                                  death_rate=death_rate,                                  inter_radius=inter_radius,                                  offspring_radius=offspring_radius,                                  d_when=d_when,                                  d_where=d_where,                                  d_pattern=d_pattern)            print('============ Generating Game Simulation ============')            eg.generate_simulation()            frames = Frames(frames=eg.simulation)            ani = animate_frames(frames, show=False, cmap='plasma')            print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)        run_recurrence(frames=frames, exp_path=exp_path)    elif GameType[game_type].value == 6:        # evo game 5        base_path = './simulations/evo_game5/' + pattern_folder        exp_path = base_path + '/experiments/' + exp_folder        if os.path.isfile(base_path + '/game_frames.npz'):            arr = np.load(base_path + '/game_frames.npz', allow_pickle=True)['arr_0']            frames = Frames(frames=arr)        else:            variables = {}            exec(open(base_path + '/config.py').read(), variables)            initial_config_file = variables['INITIAL_CONFIG_FILE']            initial_config = np.load(base_path + '/' + initial_config_file)['arr_0']            size = variables['SIZE']            n_frames = variables['N_FRAMES']            position = variables['POSITION']            fitness_matrix = variables['FITNESS_MATRIX']            death_rate = variables['DEATH_RATE']            inter_radius = variables['INTER_RADIUS']            offspring_radius = variables['OFFSPRINT_RADIUS']            d_when = variables['DISASTER_WHEN']            d_where = variables['DISASTER_WHERE']            d_pattern = variables['DISASTER_PATTERN']            eg = EvolutionaryGame5(initial_config=initial_config,                                  size=size,                                  n_frames=n_frames,                                  position=position,                                  fitness_matrix=fitness_matrix,                                  death_rate=death_rate,                                  inter_radius=inter_radius,                                  offspring_radius=offspring_radius,                                  d_when=d_when,                                  d_where=d_where,                                  d_pattern=d_pattern)            print('============ Generating Game Simulation ============')            eg.generate_simulation()            frames = Frames(frames=eg.simulation)            ani = animate_frames(frames, show=False, cmap='plasma')            print('============ Saving Game Animation ============')            ani.save(base_path + '/game.gif', writer='Pillow', fps=60)            ani.save(base_path + '/game.mp4', writer=animation.FFMpegWriter(fps=30))            np.savez_compressed(base_path + '/game_frames', frames.frames)        run_recurrence(frames=frames, exp_path=exp_path)if __name__ == "__main__":   generate_experiment(sys.argv[1:])   