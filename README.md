# Spatio-Temporal Recurrent Patterns in Dynamic Games

```
root/
  ├── .gitignore
  ├── utils.py
  ├── game_of_life.py
  ├── evolutionary_game.py
  ├── Frames.py
  ├── recurrence.py
  ├── sinusoidal_simulation.py
  ├── generate_experiment.py
  ├── README.md
  └── simulations/
          └── game_of_life/
                    └── pattern_1/
                            └── config.py
                            └── pattern_matrix.npz
                                  └── experiments/
                                          └── experiment_1/
                                                  └── config.py
          └── evo_game/
                    └── pattern_1/
                            └── config.py
                            └── pattern_matrix.npz
                                  └── experiments/
                                          └── experiment_1/
                                                  └── config.py
```

## Implement a new experiment

To generate a new analysis, add a directory under `game_of_life/` or `evo_game/` structures. Name this directory however you prefer (in the example above this directory is `pattern_1`).
Inside this directory, a `config.py` file is needed to input the parameters of the game, and a `.npz` matrix with the initial configuration of the game is also necessary.
Then, create a new experiment folder under `experiments/` and name this directory however you prefer (in the example above this directory is `experiment_1`).
Inside this directory, a `config.py` file is needed to input the parameters of the experiment.

To run the experiment, in the root directory, run:
```zsh
  python3 generate_experiment.py --game_type <game_type> --pattern_folder_name <pattern_folder_name> --exp_folder_name <exp_folder_name>
```

Where:
| Variable      | Description   |
| ------------- |:-------------:|
| <game_type>                | `GAME_OF_LIFE` or `EVO_GAME` |
| <pattern_folder_name>      | the name of the folder of the game |
| <exp_folder_name>          | the name of the folder of the experiment |

For the tree above, the command would be (for the evolutionary game):
```zsh
   python3 generate_experiment.py --game_type EVO_GAME --pattern_folder_name pattern_1 --exp_folder_name experiment_1
```

for the game of life:
```zsh
  python3 generate_experiment.py --game_type GAME_OF_LIFE --pattern_folder_name pattern_1 --exp_folder_name experiment_1
```
