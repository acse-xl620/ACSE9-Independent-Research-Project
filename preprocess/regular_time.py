import os, sys
from numpy.lib.npyio import save
import vtk, vtktools
import u2r
import numpy as np
import matplotlib.pyplot as plt
import random
from nirom_dd_tools import *
import copy

def get_clean_vtk_file(filename):
    "Removes fields and arrays from a vtk file, leaving the coordinates/connectivity information."
    vtu_data = vtktools.vtu(filename)
    clean_vtu = vtktools.vtu()
    clean_vtu.ugrid.DeepCopy(vtu_data.ugrid)
    fieldNames = clean_vtu.GetFieldNames()
    # remove all fields and arrays from this vtu
    for field in fieldNames:
        clean_vtu.RemoveField(field)
        fieldNames = clean_vtu.GetFieldNames()
        vtkdata = clean_vtu.ugrid.GetCellData()
        arrayNames = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    for array in arrayNames:
        vtkdata.RemoveArray(array)
    return clean_vtu

def plot_grid(grid_origin, grid_width, nx, ny):
    # include plot of entire coordinates with grid
    # plt.figure(figsize=(9,9))
    plt.plot(coordinates[:, 0], coordinates[:, 1], 'g.', ms=0.3,
             label='angle = {}'.format(random_angle))  # corrdinates


    for d in range(ny + 1):
        if d % 4 == 0:
            plt.plot([grid_origin[0], grid_origin[0] + grid_width[0]],
                     [grid_origin[1] + d * ddx[1], grid_origin[1] + d * ddx[1]], 'k-', lw=1.2)  # horizontal
            if ny == nx:
                plt.plot([grid_origin[0] + d * ddx[1], grid_origin[0] + d * ddx[1]],
                         [grid_origin[1], grid_origin[1] + grid_width[1]], 'k-', lw=1.2)  # vertical

    if ny != nx:
        for d in range(nx + 1):  # vertical
            if d % 4 == 0:
                plt.plot([grid_origin[0] + d * ddx[0], grid_origin[0] + d * ddx[0]],
                         [grid_origin[1], grid_origin[1] + grid_width[1]], 'k-', lw=1.2)  # vertical

def sample_(mesh, grid_origin):
    """
    Returns a snapshot matrix of shape (5,nx*ny) and
    snapshot_ae of shape (5,nx,ny) with given
    mesh and central grid origin for the starshape grid formation
    """


    s_ae = np.zeros((nx, ny))

    value_grid = u2r.simple_interpolate_from_mesh_to_grid(mesh, x_all, x_ndgln, ddx, grid_origin, nx, ny, nz,
                                                          zeros_beyond_mesh, nEl, nloc, nNodes, nScalar, nDim, 1)
    s_ae[:, :] = value_grid.reshape((nx, ny))

    return s_ae

def rotate_mesh(angle):
    theta = np.radians(angle)

    # shift coordinates so that they are centred at (0,0)
    # for i in range(coordinates.shape[0]):
    #     coordinates[i][0] -= 1.5
    #     coordinates[i][1] -= 1.5

    new_mesh = np.zeros(coordinates.shape)

    for i in range(coordinates.shape[0]):
        new_mesh[i][0] = (coordinates[i][0] - 1.5) * np.cos(theta) - (coordinates[i][1] - 1.5) * np.sin(theta)
        new_mesh[i][1] = (coordinates[i][0] - 1.5) * np.sin(theta) + (coordinates[i][1] - 1.5) * np.cos(theta)

    # rotate the velocity field as well

    return new_mesh



# snapshot_data_location = 'F:/new_VTU_POD/time_data_new/'
snapshot_data_location = 'F:/new_VTU_POD/time_data_new/'
snapshot_file_base = 'Flow_past_buildings_'

field_names = ['Velocity', 'VelocityAbsorption']
nFields = len(field_names)

xlength = 3.0
ylength = 3.0

grid_width = [0.48, 0.48]
# spacing inside small grid
nx = int(grid_width[0] * 100)
ny = nx
nz = 1
ddx = np.array((0.01, 0.01))

# get a vtu file (any will do as the mesh is not adapted)
#filename = snapshot_data_location + snapshot_file_base + '400.vtu'


nTime = 100
start = 50
Ssnapshot_ae = np.zeros((nTime,6 * 6, nx, ny, 2))
Ssnapshot_ae_info = np.zeros((nTime,6 * 6, nx, ny, 1))

for iTime in range(nTime):
    filename = snapshot_data_location + snapshot_file_base + str(start+iTime)+'.vtu'
    print(filename)
    representative_vtu = vtktools.vtu(filename)
    coordinates_org = representative_vtu.GetLocations()  # coordinates of the nodes
    coordinates = coordinates_org

    nNodes = coordinates.shape[0]  # vtu_data.ugrid.GetNumberOfPoints()
    nEl = representative_vtu.ugrid.GetNumberOfCells()
    # print('nEl', nEl, type(nEl), 'nNodes', nNodes)
    nloc = 3  # number of local nodes, ie three nodes per element (in 2D)
    # nScalar = 2 # dimension of fields , 2 = u and v
    nScalar = 1
    nDim = 2


    # get global node numbers
    x_ndgln = np.zeros((nEl * nloc), dtype=int)
    for iEl in range(nEl):
        n = representative_vtu.GetCellPoints(iEl) + 1
        x_ndgln[iEl * nloc:(iEl + 1) * nloc] = n


    # min_x = min(coordinates[:, 0])
    # max_x = max(coordinates[:, 0])
    # min_y = min(coordinates[:, 1])
    # max_y = max(coordinates[:, 1])
    # print(min_x,max_x,min_y,max_y)


    print("Generating starshape grids")
    zeros_beyond_mesh = 0
    start_point = [0, 2.5]

    nIter = 6

    for iIter in range(nIter):

        velocity_field = representative_vtu.GetField(field_names[0])[:, :nDim]  # field name 0 is velocity field
        start_point[0] = 0

        for iGrid in range(6):

            # rotate mesh and velocity
            coordinates = coordinates_org
            #coordinates = rotate_mesh(rangle)
            x_all = np.transpose(coordinates[:, 0:nDim])
            #velocity_field_g = rotate_vel(rangle, velocity_field)
            va_field = representative_vtu.GetField(field_names[1])[:,0]

            # starshape for u_vel
            u_mesh = np.zeros((1, nNodes, 1))
            u_mesh[:, :, 0] = np.transpose(velocity_field[:, 0])
            u_sae = sample_(u_mesh, start_point)
            # Ssnapshots_matrix[:nx*ny,5*iGrid:5*iGrid+5] = u_smatrix
            Ssnapshot_ae[iTime, iIter * 6 + iGrid, :, :, 0] = u_sae

            # starshape for v_vel
            v_mesh = np.zeros((1, nNodes, 1))
            v_mesh[:, :, 0] = np.transpose(velocity_field[:, 1])
            v_sae = sample_(v_mesh, start_point)
            # Ssnapshots_matrix[nx*ny:nx*ny*2,5*iGrid:5*iGrid+5] = v_smatrix
            Ssnapshot_ae[iTime, iIter * 6 + iGrid, :, :, 1] = v_sae

            # starshape for va
            va_mesh = np.zeros((1, nNodes, 1))
            va_mesh[:, :, 0] = np.transpose(va_field)
            va_sae = sample_(va_mesh, start_point)
            # Ssnapshots_matrix[nx*ny*2:nx*ny*3,5*iGrid:5*iGrid+5] = va_smatrix
            Ssnapshot_ae_info[iTime, iIter * 6 + iGrid, :, :, 0] = va_sae

            start_point[0] = start_point[0] + 0.48

        start_point[1] = start_point[1] - 0.48

# Scale values
# min_vel = np.amin(Ssnapshots_matrix[:nx*ny*2,:]) #minimum among u and v velocity
# max_vel = np.amax(Ssnapshots_matrix[:nx*ny*2,:])
min_vel = np.amin(Ssnapshot_ae)
max_vel = np.amax(Ssnapshot_ae)
vel_scaling = 1 / (max_vel - min_vel)
min_va = np.amin(Ssnapshot_ae_info)
max_va = np.amax(Ssnapshot_ae_info)

va_scaling = 1 / (max_va - min_va)
# scale snapshots for ae
Ssnapshot_ae = vel_scaling * (Ssnapshot_ae - min_vel)
Ssnapshot_ae_info = va_scaling* (Ssnapshot_ae_info - min_va)
print("min_vel = ", min_vel, " , max_vel = ", max_vel, " , vel_scaling = ", vel_scaling)

print("min_va = ", min_va, " , max_va = ", max_va, " , va_scaling = ", va_scaling)

# np.save("regular_time_velocity.npy", Ssnapshot_ae)
# np.save("regular_time_info.npy",Ssnapshot_ae_info)