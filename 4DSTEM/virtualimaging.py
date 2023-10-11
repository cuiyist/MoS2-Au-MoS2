#Initialization - import the needed packages.
import py4DSTEM
from py4DSTEM.visualize import show
import numpy as np
from numpy import savetxt
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load the data
file_data = 'sample2_ss=1nm_C2=10um_alpha=0p4_spot11_exposure50ms_CL=195mm_bin2_300kV.dm4'
py4DSTEM.io.import_file(file_data)

datacube = py4DSTEM.io.import_file(
    file_data,
    data_id = 'datacube_0'
)

# Calculate maximum diffraction pattern:
datacube.get_dp_max()
datacube.get_dp_mean()

# Estimate the radius of the BF disk, and the center coordinates
probe_semiangle, probe_qx0, probe_qy0 = py4DSTEM.process.calibration.get_probe_size(
    datacube.tree['dp_mean'].data,
)

# Overlay the estimated probe position and radius on the maximum diffraction pattern
py4DSTEM.visualize.show(
    datacube.tree['dp_max'].data, 
    scaling='log',
    circle = {
      'center':(probe_qx0, probe_qy0),
      'R': probe_semiangle,
      'alpha':0.3,
      'fill':True
    },
    intensity_range='absolute',
    vmin=3,
    vmax=10
)

# Print the estimate probe radius
print('Estimated probe radius =', '%.2f' % probe_semiangle, 'pixels')
print('Position of the central peak is','x = %.5f' % probe_qx0 )
print('Position of the central peak is','y = %.5f' % probe_qy0 )

# Create a bright field (BF) virtual detector using the the center beam position, and expanding the radius slightly.
expand_BF = 5

# Overlay the estimated probe position and radius on the maximum diffraction pattern
py4DSTEM.visualize.show(
    datacube.tree['dp_max'], 
    scaling='log',
    circle = {
      'center':(probe_qx0, probe_qy0),
      'R': probe_semiangle + expand_BF,
      'alpha':0.3,
      'fill':True
    },
    intensity_range='absolute',
    vmin=6,
    vmax=10
)

# virtual ADF
radii_DF = (50/bin_factor,80/bin_factor)
radii_Au = (116/bin_factor,126/bin_factor)

center = (probe_qx0, probe_qy0)

py4DSTEM.visualize.show(
    dataset.tree('dp_max'),
    annulus = {
        'center':[center, center],
        'radii':[radii_DF, radii_Au],
        'alpha':0.3,
        'fill':True,
        'color':['c','y'],
    },
    vmax = 0.998,
    power = 0.5,
)

# Generate virtual ADF images
dataset.get_virtual_image(
    mode = 'annulus',
    geometry = (center,radii_DF),
    name = 'dark_field',
);
dataset.get_virtual_image(
    mode = 'annulus',
    geometry = (center,radii_Au),
    name = 'dark_field_Au',
);

# Create an off-axis dark-field (DF) virtual detector using the the center beam position, and expanding the radius slightly.
qx0_DF,qy0_DF = 396,485
r_DF = 15

#show the DF detector
py4DSTEM.visualize.show(
    datacube.tree['dp_max'], 
    scaling='log',
    circle = {
      'center':(qx0_DF,qy0_DF),
      'R': r_DF,
      'alpha':0.3,
      'fill':True,
      'color':'c'
    },
    intensity_range='absolute',
    vmin=2,
    vmax=10
)
datacube.get_virtual_image(
    mode = 'circle',
    geometry = ((qx0_DF, qy0_DF), r_DF),
    name = 'virt_dark_field_01'
)
# Generate virtual DF image
py4DSTEM.visualize.show(datacube.tree['virt_dark_field_01'])
