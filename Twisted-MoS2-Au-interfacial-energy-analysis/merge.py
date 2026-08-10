import numpy as np
import os

from ase.io import read, write

#-------------------------------------------------------------------

path = "energy calculation/Nov. 5"
os.chdir(path)

bottom_layer = read("moire_supercell_8-4_bottom.cif")
cell = bottom_layer.get_cell()[:]
cell[2][2] = 28.
bottom_layer.set_cell(cell)

top_layer = read("3DAu.xyz")
'''
for atom in top_layer:
    bottom_layer.append(atom)
'''
bottom_layer_frac_coords = bottom_layer.get_scaled_positions()
top_layer_frac_coords = top_layer.get_scaled_positions()

top_layer_frac_coords += np.array([-2., 3., 0.])

bottom_layer.set_scaled_positions(bottom_layer_frac_coords)
top_layer.set_scaled_positions(top_layer_frac_coords)
'''
#center Au layer
center_Au_scaled_pos = top_layer.get_scaled_positions()[-31]
shift = np.array([0.5, 0.5, center_Au_scaled_pos[-1]]) - center_Au_scaled_pos

scaled_pos = top_layer.get_scaled_positions()
symbols = top_layer.get_chemical_symbols()

for i, symbol in enumerate(symbols):
    if symbol == "Au":
        scaled_pos[i] += shift

top_layer.set_scaled_positions(scaled_pos)
'''
bottom_layer.translate([0.9211, 0.,0.])
top_layer.translate([-0.2213, 0.5029, -14.7852])

bottom_layer.extend(top_layer)

write("moire_supercell.cif", bottom_layer)
