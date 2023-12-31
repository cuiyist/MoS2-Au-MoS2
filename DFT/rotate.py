import numpy as np
import os
from ase.io import read, write

"""
Rotate the gold layer by a certain angle relative to MoS2 lattice
"""

def rotation_matrix(axis, theta):

    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.)
    b, c, d = -axis*np.sin(theta/2.)

    R = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    return R


def rotate(coords, axis, theta):

	R = rotation_matrix(axis, theta)

	coords_new = np.dot(R, coords.T).T

	return coords_new


path = "/Users/yicui/Desktop/Stanford PhD Project/MoS2-Au-MoS2 Moire/energy calculation/Nov. 5/MoS2-3DAu"
os.chdir(path)


# bottom_layer = read("moire_supercell_2-3_bottom.cif")
# top_layer = read("moire_supercell_2-3_top.cif")
# gold_layer = read("Au_111_3x3x3.xyz")

bilayer = read("moire_supercell_MoS2-3DAu_30.cif")


symbols = bilayer.get_chemical_symbols()
Au_indices = [i for i, symbol in enumerate(symbols) if symbol == "Au"]

coords = bilayer.get_positions()
Au_coords = coords[Au_indices]

#print(Au_coords)

####################################################

shift_vec = np.array([11.05160, 11.16610, 0.])
Au_coords -= shift_vec

Au_coords_new = rotate(Au_coords, axis=np.array([0., 0., 1.]), theta=(38/180.)*np.pi)
Au_coords_new += shift_vec

for i, i_ind in enumerate(Au_indices):
    coords[i_ind] = Au_coords_new[i]
    
bilayer.set_positions(coords)

write("moire_supercell_MoS2-3DAu_.cif", bilayer)
