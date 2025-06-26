import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read, write
from ase.visualize.plot import plot_atoms
import argparse
from matplotlib.patches import FancyArrow
from ase.build import make_supercell
from ase import Atoms



def get_avg_pos(atoms_list):
    # Convert to NumPy array (shape: num_frames x num_atoms x 3)
    positions = np.array([atoms.get_positions() for atoms in atoms_list])
    # Compute the average positions over all frames
    return np.mean(positions, axis=0)  # Shape: (num_atoms, 3)


def lammps2AseAtoms(lammps_atoms, atom_type_symbol_pair):
    symbols = [atom_type_symbol_pair[key] for key in lammps_atoms.get_atomic_numbers()]
    return Atoms(symbols=symbols, positions=lammps_atoms.positions, cell=lammps_atoms.cell)


# Load trajectory from EXTXYZ file
def get_guest_pos(atoms_list, n_frame):
    positions = np.vstack([atoms[n_frame:].get_positions() for atoms in atoms_list])  # Stack all positions
    return positions


def get_guest_pos_supercell(atoms_list, n_frame, P):
    positions = np.vstack([make_supercell(atoms[n_frame:], P).get_positions() for atoms in atoms_list])  # Stack all positions
    return positions


# Compute 3D histogram (probability distribution)
def compute_distribution3D(data, bins=50, range_=None):
    hist, edges = np.histogramdd(data, bins=bins, range=range_, density=True)
    return hist, edges


# Visualize the probability distribution
def plot_distribution3D(hist, edges):
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111, projection='3d')
    # Extract bin centers
    x_edges, y_edges, z_edges = edges
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    z_centers = (z_edges[:-1] + z_edges[1:]) / 2

    # Flatten the histogram for plotting
    X, Y, Z = np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')
    values = hist.flatten()
    mask = values > 0  # Plot only occupied bins
    ax.scatter(X.flatten()[mask], Y.flatten()[mask], Z.flatten()[mask], c=values[mask], cmap='viridis', marker='o')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.title("Gas Molecule Probability Distribution")
    plt.savefig("prob_dist.png")
    #  plt.show()


# Compute 2D histogram (probability distribution)
def compute_distribution2D(data, directions, bins=50, range_=None):
    hist, x_edges, y_edges = np.histogram2d(data[:, directions[0]],
                                            data[:, directions[1]],
                                            bins=bins, range=range_, density=True)
    return hist, x_edges, y_edges

# Visualize the probability distribution in 2D
def plot_distribution2D(hist, x_edges, y_edges):
    plt.figure(figsize=(8, 4))
    plt.imshow(hist.T, origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], aspect='auto', cmap='viridis')
    plt.colorbar(label='Probability Density')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Gas Molecule Probability Distribution (2D)")
    plt.savefig("prob_dist2D.png")
    #  plt.show()

# Load MOF structure from XYZ or CIF file
def get_frame_positons(frame0):
    return frame0.get_positions()[:, :2]  # Extract only X, Y coordinates


# Visualize the probability distribution in 2D with MOF structure
def plot_distribution2DwithMOF(hist, x_edges, y_edges, frame):
    #  plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots(figsize=(8, 4))

    # Plot probability distribution
    img = ax.imshow(hist.T, origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], aspect='auto', cmap='viridis')
    plt.colorbar(img, ax=ax, label='Probability Density')  # Pass 'img' to colorbar


        # Overlay MOF structure
    #  plot_atoms(frame, ax, radii=0.3, rotation=('0x,0y,0z'))

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.title("Gas Molecule Probability Distribution (2D) with MOF Structure")
    plt.savefig("prob_dist2DwithMOF.png")
    #  plt.show()



def counter_plot_distribution2D(hist, x_edges, y_edges, frame):
    fig, ax = plt.subplots(figsize=(12, 8))
    #  fig, ax = plt.subplots(figsize=(12, 6))

    # Generate meshgrid for contour plot
    X, Y = np.meshgrid((x_edges[:-1] + x_edges[1:]) / 2, (y_edges[:-1] + y_edges[1:]) / 2)

    # Plot probability distribution using contour plot
    #  contour = ax.contourf(X, Y, hist.T, levels=20, cmap='coolwarm')
    contour = ax.contourf(X, Y, hist.T, levels=20, cmap='afmhot')
    #  plt.colorbar(contour, ax=ax, label='Probability Density')
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label(label='Increasing Probability Density', fontsize=16)
    cbar.ax.set_yticklabels([])  # Remove numbers from colorbar
    cbar.ax.tick_params(size=0)  # Remove ticks from colorbar


    # Overlay MOF structure

    #  frame_positions = frame.get_positions()[:, :2]  # Extract only X and Y coordinates
    #  ax.scatter(frame_positions[:, 0], frame_positions[:, 1], c='red', s=10, label='MOF Atoms')
    #  ax.legend()

    # as an good alternatif using ASE 
    # NOTE: there is some scaling issue
    #  plot_atoms(frame, ax, radii=0.5, offset=(-28, -1.5), rotation=('0x,0y,0z'))
    #  plot_atoms(frame, ax, radii=0.5, offset=(0, -1.5), rotation=('0x,90y,0z'))
        # Get colorbar position
    cbar_pos = cbar.ax.get_position()
    arrow_x = cbar_pos.x1 - 0.072
    arrow_y_start = cbar_pos.y0 -0.15
    arrow_y_end = cbar_pos.y1 -0.15


    # Add arrow on top of the colorbar
    fig.patches.append(FancyArrow(
        x=arrow_x, y=arrow_y_end, dx=0, dy=0.05,
        width=0.012, head_width=0.025, head_length=0.04,
        color='firebrick', transform=fig.transFigure, clip_on=False
    ))

     # Add an increasing direction arrow to the colorbar
    #  cbar.ax.annotate(
    #      "↑", xy=(0.5, 0.9), xycoords="axes fraction", fontsize=24, ha="center", va="bottom"
    #  )

    #  ax.set_xlim(0, 25)
    #  ax.set_ylim(10, 35)
    #  ax.set_xlim(-30, -5)
    #  ax.set_ylim(-45, -20)
    ax.set_xlim(-25, 5)
    ax.set_ylim(-30, -8)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.savefig("prob_dist2DwithMOF_test.png")



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-traj_file", type=str, required=True, help="")
parser.add_argument("-frame_file", type=str, required=True, help="")
args = parser.parse_args()

# Example usage
traj_file = args.traj_file
frame_file = args.frame_file

atoms_list = read(traj_file, index=slice(-750, -1, 1))  # Read all frames

replica = [2, 2, 2]
P = [[0, 0, -replica[0]], [0, -replica[1], 0], [-replica[2], 0, 0]]

frame0 = read(frame_file)
#  print(atoms_list[0].get_positions())
n_frame = len(frame0)
frame = atoms_list[-1][:n_frame]
#  write("test.extxyz", make_supercell(atoms_list[-1][n_frame:], P))
#  quit()

#  pos_data = get_guest_pos(atoms_list, n_frame)
pos_data = get_guest_pos_supercell(atoms_list, n_frame, P)
#  frame = atoms_list[-1][:n_frame]

# get avg structure
#  frame_atoms_list = [atoms[:n_frame] for atoms in atoms_list]
#  frame.set_positions(get_avg_pos(frame_atoms_list))
#  print(frame.get_positions())
#  write("frame.cif", frame)
#  quit()

# NOTE before lammps2AseAtoms due to wrap position problem
frame = make_supercell(frame, P)

if traj_file.endswith(".lammpstrj"):
    atom_type_symbol_pair = {1:"Mg", 2:"O", 3:"C", 4:"H"}
    frame = lammps2AseAtoms(frame, atom_type_symbol_pair)

write("frame.cif", frame)
#  frame_positions = get_frame_positons(frame0)
#  histogram, bin_edges = compute_distribution3D(pos_data, bins=50)
#  plot_distribution3D(histogram, bin_edges)
histogram, x_edges, y_edges = compute_distribution2D(pos_data, directions=(1,2), bins=50)
#  histogram, x_edges, y_edges = compute_distribution2D(pos_data, directions=(0,1), bins=50)
#  plot_distribution2D(histogram, x_edges, y_edges)
#  plot_distribution2DwithMOF(histogram, x_edges, y_edges, frame_positions)
counter_plot_distribution2D(histogram, x_edges, y_edges, frame)



