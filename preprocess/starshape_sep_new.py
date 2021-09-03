import os,sys
import vtk,vtktools

sys.path.append("E:\new_VTU_POD")
import u2r
import numpy as np
import matplotlib.pyplot as plt
from nirom_dd_tools import *
import random

def write_sing_values(s_values):
    f= open('singular_values.dat',"w+")
    f.write('# index, s_values, normalised s_values, cumulative energy \n' )
    for k in range(len(s_values)):
        #f.write('# field: %s\n' % field[k])
        total = 0.0
        s_values = s_values[k]
        for i in range(len(s_values)):
            total = total + s_values[i]*s_values[i]

        running_total = 0.0
        for i in range(len(s_values)):
            running_total = running_total + s_values[i]*s_values[i]
            f.write ('%d %g %g %18.10g \n' % (i, s_values[i], s_values[i]/s_values[0], running_total/total) )
    f.close()
    return


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
        vtkdata=clean_vtu.ugrid.GetCellData()
        arrayNames = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    for array in arrayNames:
        vtkdata.RemoveArray(array)
    return clean_vtu

#(nNodes, reconstruction_on_mesh[iTime*nScalar:(iTime+1)*nScalar,:], template_vtu, original_data[0][iTime*nDim:(iTime+1)*nDim], iTime)
def create_vtu_file(path, nNodes, value_mesh_twice_interp, filename, orig_vel, iTime):
    velocity_field = np.zeros((nNodes,3))
    velocity_field[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) # streamwise component only

    difference = np.zeros((nNodes,3))
    difference[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) - orig_vel # streamwise component only
    difference = difference / np.max(velocity_field)

    clean_vtk = get_clean_vtk_file(filename)
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = path + 'recon_' + str(iTime) + '.vtu'
    new_vtu.AddField('Velocity',velocity_field)
    new_vtu.AddField('Original',orig_vel)
    new_vtu.AddField('Velocity_diff',difference)
    new_vtu.Write()
    return

def create_vtu_file_timelevel(nNodes, value_mesh_twice_interp, template_vtu, iTime):
    velocity_field = np.zeros((nNodes,3))
    velocity_field[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) # streamwise component only

#    difference = np.zeros((nNodes,3))
#    difference[:,0:nDim] = np.transpose(value_mesh_twice_interp[0:nDim,:]) - orig_vel # streamwise component only

    clean_vtk = get_clean_vtk_file(template_vtu)
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = 'reconstructed_' + str(iTime) + '.vtu'
    new_vtu.AddField('Velocity',velocity_field)
    #new_vtu.AddField('Velocity_diff',difference)
    new_vtu.Write()
    return


def get_grid_end_points(grid_origin,grid_width,iGrid ):
    return np.array(( grid_origin[0]+iGrid*grid_width[0], grid_origin[1] +iGrid*grid_width[1]))

def get_grid_end_points_new(grid_origin,grid_width ):
    return np.array((grid_origin[0], grid_origin[1]))


def plot_grid(grid_origin, grid_width,nx,ny,coordinates):
    plt.plot(coordinates[:,0],coordinates[:,1],'b.',ms = 0.3 )
    for d in range(ny):
        if d%5==0:
            plt.plot([grid_origin[0], grid_origin[0] + grid_width[0]],[grid_origin[1]+d*ddx[1], grid_origin[1]+d*ddx[1]],'k-', lw = 1.2)
            if ny==nx:
                plt.plot([grid_origin[0]+d*ddx[1], grid_origin[0] +d*ddx[1]],[grid_origin[1], grid_origin[1]+ grid_width[1]],'k-', lw = 1.2)

    if ny!=nx:
        for d in range(nx):
            plt.plot([grid_origin[0]+d*ddx[0], grid_origin[0] +d*ddx[0]],[grid_origin[1], grid_origin[1]+ grid_width[1]],'k-', lw = 1.2)
    #plt.grid(':')
    plt.tight_layout()

def random_select(coordinates):
    min_x = min(coordinates[:, 0])
    max_x = max(coordinates[:, 0])
    min_y = min(coordinates[:, 1])
    max_y = max(coordinates[:, 1])

    plt.plot(min_x + 0.5, min_y + 0.5)
    plt.plot(min_x + 0.5, max_y - 1.0)
    plt.plot(max_x - 1.0, min_y + 0.5)
    plt.plot(max_x - 1.0, max_y - 1.0)
    grid_orgin = [3,3]

    while np.sqrt(grid_orgin[0]**2+grid_orgin[1]**2>=1.5):
        grid_orgin = [random.uniform(min_x+0.5,max_x-1.0),random.uniform(min_y+0.5,max_y-1.0)]

    return grid_orgin

def rotate(angle,coordinate):
    theta = np.radians(angle)

    new_mesh = np.zeros(coordinate.shape)
    for i in range(coordinate.shape[0]):
        new_mesh[i][0] = (coordinate[i][0] - 1.5)*np.cos(theta) - (coordinate[i][1]-1.5)*np.sin(theta)
        new_mesh[i][1] = (coordinate[i][0] - 1.5)*np.sin(theta) + (coordinate[i][1]-1.5)*np.cos(theta)
    return new_mesh

def rotate_velo(angle,velocity_field,coordinate):
    theta = np.radians(angle)
    new_mesh = np.zeros(velocity_field.shape)

    for i in range(coordinate.shape[0]):
        new_mesh[i][0] = (velocity_field[i][0] * np.cos(theta) - velocity_field[i][1] * np.sin(theta))
        new_mesh[i][1] = (velocity_field[i][0] * np.sin(theta) + velocity_field[i][1] * np.cos(theta))

    return new_mesh
## settings

snapshot_data_location = 'F:/new_VTU_POD/FPI/'
snapshot_file_base = 'Flow_past_buildings_Re1_training_'

# snapshot_data_location = 'F:/new_VTU_POD/FPC_Re3900_2D_CG_old/'
# snapshot_file_base = 'fpc_2D_Re3900_CG_'

nTime = 1
field_names = ['Velocity', 'VelocityAbsorption']
nFields = len(field_names)

# xlength = 2.2  # should work this out from coordinates, below
# ylength = 0.41 #
xlength = 3.0
ylength = 3.0




# grid_origin = [0.0,0.0]
grid_width = [0.47,0.47]
nx = int(grid_width[0]*100+1)
ny = nx
nz = 1
ddx = np.array((0.01,0.01))
#ddx = np.array((xlength/(nGrids*(nx-1)),ylength/(ny-1)))
print ('ddx', ddx)


#get a vtu file (any will do as the mesh is not adapted)
# filename = snapshot_data_location + snapshot_file_base + '0.vtu'
filename = snapshot_data_location + snapshot_file_base + '400.vtu'
representative_vtu = vtktools.vtu(filename)
coordinates_org = representative_vtu.GetLocations()
#coordinates = coordinates_org
print("shape",coordinates_org.shape)
nNodes = coordinates_org.shape[0] # vtu_data.ugrid.GetNumberOfPoints()
nEl = representative_vtu.ugrid.GetNumberOfCells()
print('nEl', nEl, type(nEl), 'nNodes', nNodes)
print(nNodes)

#nNodes = 3571
#nEl = 6850
nloc = 3 # number of local nodes, ie three nodes per element (in 2D)
nScalar = 2 # dimension of fields
nDim = 2 # dimension of problem (no need to interpolate in the third dimension)
nScalar_test = 1



# get global node numbers
x_ndgln = np.zeros((nEl*nloc), dtype=int)
for iEl in range(nEl):
    n = representative_vtu.GetCellPoints(iEl) + 1
    x_ndgln[iEl*nloc:(iEl+1)*nloc] = n






# rangle_list= np.load(("F:/new_VTU_POD/rangle_list.npy"))
# origin_list= np.load(("F:/new_VTU_POD/rangle_list.npy"))


nStar = 5000
snapshot_matrix_recon = np.zeros((nx*ny*3, nStar*5))
snapshot_matrix_recon_velocity = np.zeros((nx*ny*2, nStar*5))
snapshot_matrix_recon_info = np.zeros((nx*ny, nStar*5))
for iStar in range (nStar):
    print("istar :",iStar)
    random_angle = random.randint(0, 360)
    coordinates_new = rotate(random_angle, coordinates_org)
    velocity_field = representative_vtu.GetField(field_names[0])[:, :nDim]
    velocity_field_new = rotate_velo(random_angle, velocity_field, coordinates_new)
    value_field = representative_vtu.GetField(field_names[1])[:, 0]
    grid_origin = random_select(coordinates_new)

    origin_array = []
    origin_array.append(grid_origin)
    origin_array.append([grid_origin[0] - grid_width[0], grid_origin[1]])
    origin_array.append([grid_origin[0] + grid_width[0], grid_origin[1]])
    origin_array.append([grid_origin[0], grid_origin[1] + grid_width[1]])
    origin_array.append([grid_origin[0], grid_origin[1] - grid_width[1]])
    #print("origin_array :",origin_array)
    for i in range (5):
        block_x_start = get_grid_end_points_new(origin_array[i], grid_width)
        zeros_on_mesh = 0
        x_all = np.transpose(coordinates_new[:, 0:nDim])

        # interpolate velocity
        u_mesh = np.zeros((1, nNodes, 1))
        u_mesh[:, :, 0] = np.transpose(velocity_field[:, 0])

        u_value_grid = u2r.simple_interpolate_from_mesh_to_grid(u_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_on_mesh, nEl, nloc, nNodes, nScalar_test, nDim, 1)
        snapshot_matrix_recon[:nx * ny * 1, 5*iStar + i] = u_value_grid.reshape(-1)
        snapshot_matrix_recon_velocity[:nx * ny * 1, 5 * iStar + i] = u_value_grid.reshape(-1)

        v_mesh = np.zeros((1, nNodes, 1))
        v_mesh[:, :, 0] = np.transpose(velocity_field[:, 1])

        v_value_grid = u2r.simple_interpolate_from_mesh_to_grid(v_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_on_mesh, nEl, nloc, nNodes, nScalar_test, nDim, 1)
        snapshot_matrix_recon[nx * ny:nx * ny * 2, 5*iStar + i] = v_value_grid.reshape(-1)
        snapshot_matrix_recon_velocity[nx * ny:nx * ny * 2, 5 * iStar + i] = v_value_grid.reshape(-1)

        info_mesh = np.zeros((1, nNodes, 1))
        info_mesh[:, :, 0] = np.transpose(value_field)

        info_value_grid = u2r.simple_interpolate_from_mesh_to_grid(info_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_on_mesh, nEl, nloc, nNodes, nScalar_test, nDim, 1)
        snapshot_matrix_recon[nx * ny * 2:nx * ny * 3, 5*iStar + i] = info_value_grid.reshape(-1)
        snapshot_matrix_recon_info[:nx * ny * 1, 5 * iStar + i] = info_value_grid.reshape(-1)


minv = np.amin(snapshot_matrix_recon[:nx*ny*2,:])
maxv = np.amax(snapshot_matrix_recon[:nx*ny*2,:])
vel_scaling = 1/(maxv-minv)
va_scaling = 1e-5

snapshot_matrix_recon_velocity[:nx*ny*2,:] = vel_scaling*(snapshot_matrix_recon_velocity[:nx*ny*2,:]-minv)
snapshot_matrix_recon_info[:nx*ny,:] = va_scaling*snapshot_matrix_recon_info[:nx*ny,:]


star_velocity = np.zeros((5*nStar, nx, ny, 2))
for iIter in range(5*nStar):

    star_velocity[iIter, :, :, 0] = snapshot_matrix_recon_velocity[:nx * ny * 1, iIter].reshape(48, 48)
    star_velocity[iIter, :, :, 1] = snapshot_matrix_recon_velocity[nx * ny:nx * ny * 2, iIter].reshape(48, 48)

star_info = np.zeros((5*nStar, nx, ny, 1))
for iIter in range(5*nStar):

    star_info[iIter, :, :, 0] = snapshot_matrix_recon_info[:nx * ny * 1, iIter].reshape(48, 48)

np.save("star_velocity.npy", np.array(star_velocity))
np.save("star_info.npy", np.array(star_info))