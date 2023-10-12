import py4DSTEM
from py4DSTEM.visualize import show
import numpy as np
from numpy import savetxt
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file_data = 'sample2_ss=1nm_C2=10um_alpha=0p4_spot11_exposure50ms_CL=195mm_bin2_300kV.dm4'
bin_factor = 2

dataset = py4DSTEM.import_file(
    file_data,
)

# Apply binning
dataset.bin_Q(bin_factor)

# HAADF detector mask
dataset.get_beamstop_mask(
    threshold = 0.1,
    distance_edge = 4,
    include_edges = False,
    sigma = 2,
);
fitting_mask = dataset.tree["mask_beamstop"].data == False
py4DSTEM.show(
    fitting_mask,
    intensity_range='absolute',
    vmin=0,
    vmax=1,
)

# stain mapping using wholepatternfit
bound = 4.0 / bin_factor
radius_init = 4.56 / bin_factor
width_init = 2.5 / bin_factor

center_guess = (115.64, 143.43)

WPF = py4DSTEM.process.wholepatternfit.WholePatternFit(
    dataset,
    x0=(center_guess[0], 1.0/bin_factor),
    y0=(center_guess[1], 1.0/bin_factor),
    mask = fitting_mask,
    
    meanCBED = dataset.tree('dp_mean').data * fitting_mask,
)

model_list = []

model_list.append(
    py4DSTEM.process.wholepatternfit.DCBackground(
        background_value=(100, 0, None)
    )
)
model_list.append(py4DSTEM.process.wholepatternfit.GaussianBackground(
    sigma=(16.0 / bin_factor,1,None),         
    intensity=(100,0,None)
))

# parameters for MoS2 lattice 1
model_list.append(
    py4DSTEM.process.wholepatternfit.SyntheticDiskLattice(
        WPF,
        ux=(63.79/bin_factor, bound),
        uy=(3.81/bin_factor, bound),
        vx=(27.15/bin_factor, bound),
        vy=(56.34/bin_factor, bound),
        disk_radius=(radius_init, 1e-4),
        disk_width=(width_init, 1e-4),
        u_max = 3,
        v_max = 3,
        intensity_0=(100, 0, None),
        refine_radius=True,
        refine_width=True,
        global_center=True,
        include_indices = [
            
            [0,-4],
            [1,-4],
            [2,-4],
            [3,-4],

            [-1,-3],
            [0,-3],
            [1,-3],
            [2,-3],
            [3,-3],
                        
            [-2,-2],
            [-1,-2],
            [0,-2],
            [1,-2],
            [2,-2],
            [3,-2],

            [-3,-1],
            [-2,-1],
            [-1,-1],
            [0,-1],
            [1,-1],
            [2,-1],
            [3,-1],
            
            [-3,0],
            [-2,0],
            [-1,0],
            [0,0],
            [1,0],
            [2,0],
            [3,0],
                        
            [-4,1],
            [-3,1],
            [-2,1],
            [-1,1],
            [0,1],
            [1,1],
            [2,1],
            
            [-4,2],
            [-3,2],
            [-2,2],
            [-1,2],
            [0,2],
            [1,2],
            
            [-3,3],
            [-2,3],
            [-1,3],
            [0,3],
        ],
        name="MoS2 Lattice 1",
    ),
)

# parameters for MoS2 lattice 2
model_list.append(
    py4DSTEM.process.wholepatternfit.SyntheticDiskLattice(
        WPF,
        ux=(64.04/bin_factor, bound),
        uy=(-0.66/bin_factor, bound),
        vx=(31.07/bin_factor, bound),
        vy=(54.05/bin_factor, bound),
        disk_radius=(radius_init, 1e-4),
        disk_width=(width_init, 1e-4),
        u_max = 3,
        v_max = 3,
        intensity_0=(100, 0, None),
        refine_radius=True,
        refine_width=True,
        global_center=True,
        include_indices = [
            [0,-4],
            [1,-4],
            [2,-4],
            [3,-4],

            [-1,-3],
            [0,-3],
            [1,-3],
            [2,-3],
            [3,-3],
                        
            [-2,-2],
            [-1,-2],
            [0,-2],
            [1,-2],
            [2,-2],
            [3,-2],

            [-3,-1],
            [-2,-1],
            [-1,-1],
            [0,-1],
            [1,-1],
            [2,-1],
            [3,-1],
            
            [-3,0],
            [-2,0],
            [-1,0],
            [0,0],
            [1,0],
            [2,0],
            [3,0],
                        
            [-4,1],
            [-3,1],
            [-2,1],
            [-1,1],
            [0,1],
            [1,1],
            [2,1],
            
            [-4,2],
            [-3,2],
            [-2,2],
            [-1,2],
            [0,2],
            [1,2],
            
            [-3,3],
            [-2,3],
            [-1,3],
            [0,3],
        ],
        name="MoS2 Lattice 2",
    ),
)

# parameters for Au lattice
model_list.append(
    py4DSTEM.process.wholepatternfit.SyntheticDiskLattice(
        WPF,
        ux=(109.33/bin_factor, bound),
        uy=(-57.35/bin_factor, bound),
        vx=(102.94/bin_factor, bound),
        vy=(62.17/bin_factor, bound),
        disk_radius=(radius_init, 1e-4),
        disk_width=(width_init, 1e-4),
        u_max = 2,
        v_max = 2,
        intensity_0=(100, 0, None),
        refine_radius=True,
        refine_width=True,
        global_center=True,
        include_indices = [
            
            [0,-2],
            [1,-2],
            [2,-2],
            
            [-1,-1],
            [0,-1],
            [1,-1],
            [2,-1],
            
            [-1,0],
            [0,0],
            [1,0],
            
            [-2,1],
            [-1,1],
            [0,1],
        ],
        name="Au Lattice",
    ),
)

WPF.add_model_list(model_list)

opt = WPF.fit_to_mean_CBED(xtol=1e-8)
WPF.accept_mean_CBED_fit()
WPF.show_lattice_points(figsize=(10, 10))
WPF.show_model_grid(figsize=(14, 9))

fit_all, fit_metrics = WPF.fit_all_patterns(
    max_nfev = 2, 
    xtol = 1e-4,
)

# generate strain maps for all lattices
g_maps = WPF.get_lattice_maps()
for gm in g_maps:
    g_ref = py4DSTEM.process.latticevectors.get_reference_g1g2(
        gm, np.ones(gm.data.shape[:2], dtype=np.bool_)
    )
    strain = py4DSTEM.process.latticevectors.get_strain_from_reference_g1g2(gm, *g_ref)
    print(gm.name, flush=True)
    py4DSTEM.visualize.show_strain(
        strain,
        [-0.7, 0.7], 
        [-0.2, 0.2], 
        layout=1, 
        axes_x0=5, 
        axes_y0=5, 
        figsize=(20, 8),
    )
