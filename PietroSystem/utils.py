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