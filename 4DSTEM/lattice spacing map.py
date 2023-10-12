import py4DSTEM
from py4DSTEM.visualize import show
import numpy as np
from numpy import savetxt
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file_data = 'sample2_ss=1nm_C2=10um_alpha=0p4_spot11_exposure50ms_CL=195mm_bin2_300kV.dm4'
py4DSTEM.io.import_file(file_data) #sample 2

# Load the data
datacube = py4DSTEM.io.import_file(
    file_data,
    data_id = 'datacube_0'
)

#initialize the size of the data
shape=datacube.shape
Rxdim=shape[0]
Rydim=shape[1]
Qxdim=shape[2]
Qydim=shape[3]

# Calculate maximum diffraction pattern:
datacube.get_dp_max()
datacube.get_dp_mean()

#Fit lattice spacing
dc_1=np.zeros((Qxdim,Qydim))
for i in range(0,Qxdim):
    for j in range(0,Qydim):
        for a in range(28,31):
            for b in range(24,27):
                dc_1[i,j]+=datacube[a,b,i,j]  #average 16 pixels

# mask the diffraction spot and fit its position using Gaussian function
image_size=230

qx0_DF,qy0_DF = image_size+98,image_size+183 # max value = 70, zero strain = 0.00490

r_DF = 10
center_x0=round(probe_qx0)
center_y0=round(probe_qy0)

py4DSTEM.visualize.show(
    dc_1[center_x0-image_size:center_x0+image_size,center_y0-image_size:center_y0+image_size],
    scaling='log',
    circle = {
      'center':(qx0_DF,qy0_DF),
      'R': r_DF,
      'alpha':0.3,
      'fill':True,
      'color':'c'
    },
    intensity_range='absolute',
    vmin=4,
    vmax=10
)

center_x1=center_x0-image_size+qx0_DF
center_y1=center_y0-image_size+qy0_DF

dc_dfrt=dc_1[center_x1-r_DF:center_x1+r_DF,center_y1-r_DF:center_y1+r_DF]
py4DSTEM.visualize.show(
    dc_dfrt,
    scaling='log',
    circle = {
      'center':(qx0_DF,qy0_DF),
      'R': r_DF,
      'alpha':0.3,
      'fill':True,
      'color':'c'
    },
    intensity_range='absolute',
    vmin=6.8,
    vmax=8
)

def gaussian(xy,x0,y0,sigma2,A):
    x=xy[0]
    y=xy[1]
    result=A*np.exp(-((x-x0)**2+(y-y0)**2)/(2 *sigma2**2))
    return result.ravel()

max_value=np.amax(dc_dfrt)
max_indices=np.unravel_index(np.argmax(dc_dfrt),dc_dfrt.shape)

x0,y0=max_indices[1],max_indices[0]
initial_guess=[np.float32(x0),np.float32(y0),2.,max_value]

x=np.arange(dc_dfrt.shape[1],dtype=np.float32)
y=np.arange(dc_dfrt.shape[0],dtype=np.float32)
X,Y=np.meshgrid(x,y)

params_dfrt,pcov=curve_fit(gaussian,np.array([X,Y]),dc_dfrt.ravel())

# mask the central spot and fit its position using Gaussian function
image_size=200
qx2_DF,qy2_DF = image_size,image_size
r_C = 10
center_x0=round(probe_qx0)
center_y0=round(probe_qy0)
py4DSTEM.visualize.show(
    dc_1[center_x0-image_size:center_x0+image_size,center_y0-image_size:center_y0+image_size],
    scaling='log',
    circle = {
      'center':(qx2_DF,qy2_DF),
      'R': r_C,
      'alpha':0.3,
      'fill':True,
      'color':'r'
    },
    intensity_range='absolute',
    vmin=7,
    vmax=10
)
dc_cspot=dc_1[center_x0-r_C:center_x0+r_C,center_y0-r_C:center_y0+r_C]
py4DSTEM.visualize.show(
    dc_cspot,
    scaling='log',
    intensity_range='absolute',
    vmin=9,
    vmax=11
)

max_value=np.amax(dc_cspot)
max_indices=np.unravel_index(np.argmax(dc_cspot),dc_cspot.shape)

x0,y0=max_indices[1],max_indices[0]
initial_guess=[np.float32(x0),np.float32(y0),2.,max_value]

x=np.arange(dc_cspot.shape[1],dtype=np.float32)
y=np.arange(dc_cspot.shape[0],dtype=np.float32)
X,Y=np.meshgrid(x,y)

params_cspot,pcov=curve_fit(gaussian,np.array([X,Y]),dc_cspot.ravel())

center_cspot=[center_x0-r_C+params_cspot[1],center_y0-r_C+params_cspot[0]]
center_dfrt=[center_x1-r_DF+params_dfrt[1],center_y1-r_DF+params_dfrt[0]]
dx=(center_dfrt[0]-center_cspot[0])
dy=(center_dfrt[1]-center_cspot[1])
dr=(dx**2+dy**2)**0.5
print('dx= %.5f' % dx)
print('dy= %.5f' % dy)
print('dr= %.5f' % dr)

# calculate distence between the masked diffraction spot and the central spot
result_Q2=np.zeros((Rxdim,Rydim))

for i in range(0,Rxdim):
    for j in range(0,Rydim):
        dc=datacube[i][j][:][:]
        dc_dfrt0=dc[center_x1-r_DF:center_x1+r_DF,center_y1-r_DF:center_y1+r_DF]
        dc_cspot0=dc[center_x0-r_C:center_x0+r_C,center_y0-r_C:center_y0+r_C]
        
        max_value=np.amax(dc_cspot0)
        max_indices=np.unravel_index(np.argmax(dc_cspot0),dc_cspot0.shape)
        x0,y0=max_indices[1],max_indices[0]
        initial_guess=[np.float32(x0),np.float32(y0),8.,max_value]

        x=np.arange(dc_cspot0.shape[0],dtype=np.float32)
        y=np.arange(dc_cspot0.shape[1],dtype=np.float32)
        X,Y=np.meshgrid(x,y)

        #params_cspot,pcov=curve_fit(gaussian,np.array([X,Y]),dc_cspot.ravel())
        
        max_value2=np.amax(dc_dfrt0)
        max_indices2=np.unravel_index(np.argmax(dc_dfrt0),dc_dfrt0.shape)
        x0,y0=max_indices2[1],max_indices2[0]
        initial_guess2=[np.float32(x0),np.float32(y0),2.,max_value]

        x2=np.arange(dc_dfrt0.shape[0],dtype=np.float32)
        y2=np.arange(dc_dfrt0.shape[1],dtype=np.float32)
        X2,Y2=np.meshgrid(x2,y2)
        
        try:
            params_dfrt,pcov=curve_fit(gaussian,np.array([X2,Y2]),dc_dfrt0.ravel())
            params_cspot,pcov=curve_fit(gaussian,np.array([X,Y]),dc_cspot0.ravel())
        
        #result_Q[i][j][0:3]=[center_x0-r_C+params_cspot[0],center_y0-r_C+params_cspot[1],center_x1-r_DF+params_dfrt[0],center_y1-r_DF+params_dfrt[1]]
        except RuntimeError:
            #print('No spot found')
            result_Q2[i][j]=0.0049
            pass
        else:
            params_dfrt,pcov=curve_fit(gaussian,np.array([X2,Y2]),dc_dfrt0.ravel())
            params_cspot,pcov=curve_fit(gaussian,np.array([X,Y]),dc_cspot0.ravel())
        
            #result_Q[i][j][0:3]=[center_x0-r_C+params_cspot[0],center_y0-r_C+params_cspot[1],center_x1-r_DF+params_dfrt[0],center_y1-r_DF+params_dfrt[1]]
            center_cspot=[center_x0-r_C+params_cspot[1],center_y0-r_C+params_cspot[0]]
            center_dfrt=[center_x1-r_DF+params_dfrt[1],center_y1-r_DF+params_dfrt[0]]
            dx=(center_dfrt[0]-center_cspot[0])
            dy=(center_dfrt[1]-center_cspot[1])
            dr=(dx**2+dy**2)**0.5
            result_Q[i][j][0]=np.float32(dx)
            result_Q[i][j][1]=np.float32(dy)
            result_Q[i][j][2]=np.float32(dr)
            
            if max_value2>70:
                result_Q2[i][j]=1/np.float32(dr)
            else:
                result_Q2[i][j]=0.0049
    print('%d' %i)

np.amax(result_Q2)
x=np.arange(result_Q2.shape[0],dtype=np.float32)
y=np.arange(result_Q2.shape[1],dtype=np.float32)
X,Y=np.meshgrid(x,y)

# print lattice spacing map
fig=plt.figure()
pos=plt.imshow(
    result_Q2,
    cmap='jet',
    vmin=0.0048,
    vmax=0.0050
)
fig.colorbar(pos)
