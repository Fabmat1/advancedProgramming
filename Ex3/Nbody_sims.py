import os
import subprocess
import cv2
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from GravitySimulationClass import GravitySimulation
import numpy as np


def square(x, a, c):
    return a + c * x ** 2


def generate_video_from_images(directory, prefix, output_name, fps=30):
    images = [img for img in os.listdir(directory) if img.endswith(".png") and img.startswith(prefix)]
    images.sort(key=lambda x: int(x.split(".")[0].split("_")[-1]))

    # Determine the width and height from the first image
    frame = cv2.imread(os.path.join(directory, images[0]))

    height, width, layers = frame.shape

    video = cv2.VideoWriter(output_name, 0, fps, (width, height))

    # add images to video
    for image in images:
        video.write(cv2.imread(os.path.join(directory, image)))

    cv2.destroyAllWindows()
    video.release()

def makeplots_cpp(filename, dt, prefix='cpp_plot'):
    df = pd.read_csv(filename, delimiter=';')

    # Set the directory where you want to save the plot images
    output_dir = 'output/'

    # Create the directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each row of the dataframe
    for i in range(len(df)):
        # Extract particle positions
        positions = df.iloc[i]
        n = len(positions)
        xs = [positions[j] for j in range(0, n, 3)]
        ys = [positions[j] for j in range(1, n, 3)]
        zs = [positions[j] for j in range(2, n, 3)]

        # Create a plot
        fig = plt.figure(figsize=(10, 10))
        ax = plt.axes(projection="3d")
        ax.scatter(xs, ys, zs)

        # Set title and labels
        ax.set_title(f'Time: {round(dt*(i+1), 1)}', pad=20)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.set_xlim((-10, 10))
        ax.set_ylim((-10, 10))
        ax.set_zlim((-10, 10))

        plt.tight_layout()

        # Save the plot
        fig.savefig(os.path.join(output_dir, f'{prefix}_{i+1}.png'))

        # Close the plot to free up memory
        plt.close(fig)



os.makedirs("output/", exist_ok=True)

print("Analyzing C++ results...")

try:
    makeplots_cpp("../output/cpp_Nbody_coords.csv", 0.1)
    generate_video_from_images("output", "cpp_plot", "animation_cpp.mp4")
except:
    print("Please build and execute the C++ target 'quickNbody' before executing this script!")

print("Running the Simulation with dt = 0.01; T_max = 30; N = 512...")

sim = GravitySimulation(OutputPrefix="output/forvideo")
sim.InitializeHernquistHalo(Nparticles=512)
sim.RunSimulation(TimebetweenPlots=0.1)

generate_video_from_images("output", "forvideo", "animation.mp4")

Ns = [2 ** (i + 5) for i in range(7)]
avg_times = []

for N in Ns:
    print(f"Running the Simulation with dt = 0.01; T_max = 2; N = {N}...")
    sim = GravitySimulation(OutputPrefix=f"output/{N}")
    sim.InitializeHernquistHalo(Nparticles=N)
    sim.RunSimulation(tmax=2)
    avg_times.append(np.mean(sim.AccelerationTime))

# Fit a polynomial of grade 2 to the time increase, to show N^2 relationship

params, _ = curve_fit(square,
                      np.array(Ns),
                      np.array(avg_times),
                      bounds=[[0, 0],
                              [np.inf, np.inf]])

print(
    f"The time it takes to calculate a simulation step increases as {params[0]}+{params[1]}*N^2")

nspace = np.linspace(Ns[0], Ns[-1], 100)
plt.scatter(Ns, avg_times, color="darkred", zorder=5)

# Plot times from c++
if os.path.isfile("../output/cpp_Nbody_times.csv"):
    cpp_times = np.loadtxt("../output/cpp_Nbody_times.csv").flatten()
    plt.scatter(Ns, cpp_times, color="darkgreen", zorder=4)

plt.plot(nspace, square(nspace, *params), color="darkgrey", zorder=1)#

plt.loglog()
plt.ylabel("Time per simulation step [s]")
plt.xlabel("Number of particles")

try:
    plt.legend(["Measured Timedelta", "Measured Timedelta (C++)", "Fit"])
except:
    plt.legend(["Measured Timedelta", "Fit"])

plt.tight_layout()
plt.savefig("time_scaling.png", dpi=300)
plt.show()
