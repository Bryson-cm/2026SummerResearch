import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from matplotlib.animation import PillowWriter
import matplotlib.animation as animation
from matplotlib.animation import HTMLWriter
import sys
import importlib
#import shutil
#print(shutil.which("ffmpeg"))
##################################################
# Universal Parameters
me = 9.1093837e-31
c = 299792458.0
ne = 1e15 * 1e6
e0 = 8.85418782e-12
e = 1.60217663e-19
##################################################
# Define Parameters
wp = np.sqrt((ne*e**2)/(me*e0))
gamma = 110
k = 0.475 * me * wp**2
rb = 0.65 * c/wp
px = 110 * me * c
f0 = px**2/(gamma*me*k)
##################################################
#Define the function for the cylindrical wakefield
def f_cyl(i):
    return f0/(2*rb) * 1/np.sqrt(1-(i/rb)**2)

##################################################


input_module = importlib.import_module(sys.argv[1])

N = input_module.N                 #Number of Electrons
distance = input_module.distance         #Distance from Wakefield

##################################################


# t values to animate
t_values = np.linspace(-1, 0, 100)

fig, ax = plt.subplots()




def update(frame):

    t = t_values[frame]

    y0rb = np.linspace(t, 0.99, N)

    electron = (-(y0rb*rb/f_cyl(y0rb*rb))*distance+ y0rb*rb)

    arr = electron/rb

    ax.clear()
    ax.hist(arr, bins=60, edgecolor='black')

    ax.set_xlabel('y0/rb')
    ax.set_ylabel('# of electrons')
    ax.set_xlim(-1, 1)
    ax.set_title(f't = {t:.3f}')

    ax.text(
    0.02, 0.95,
    rf'$y_{{0,\min}}/r_b = {t:.3f}$',
    transform=ax.transAxes,
    fontsize=14,
    verticalalignment='top',
    bbox=dict(boxstyle='round', alpha=0.8) 
    )

ani = FuncAnimation(
    fig,
    update,
    frames=len(t_values),
    interval=100,
    repeat=True
)

plt.close()

html_content = ani.to_jshtml()

# Save to file
with open("electron_histogram.html", "w") as f:
    f.write(html_content)

# Display inline in notebook
HTML(html_content)
