import os,sys
import vtk,vtktools

sys.path.append("E:\new_VTU_POD")
import u2r
import numpy as np
import matplotlib.pyplot as plt
from nirom_dd_tools import * 
import random
# from keras.losses import mean_squared_error

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
snapshot_file_base = 'Flow_past_implicit_block_'

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
filename = snapshot_data_location + snapshot_file_base + '84.vtu'
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


nGrids = 500


snapshot_matrix = np.zeros((nx*ny*3,nGrids))
snapshot_output_data = np.zeros((nGrids, nx, ny, 3))

# segmentation
snapshot_velocity = np.zeros((nGrids, nx, ny, 2))
snapshot_info = np.zeros((nGrids, nx, ny, 1))

for iGrid in range (nGrids):
    #snapshot_matrix = np.zeros((nx * ny * 3, nGrids))


    random_angle = random.randint(0,360)
    #print("random angle is : ", random_angle)
    #print("org",coordinates_org)
    coordinates_new = rotate(random_angle,coordinates_org)
    #print("new",coordinates_new)
    velocity_field = representative_vtu.GetField(field_names[0])[:,:nDim]
    #print("fitst velocity",velocity_field[0])
    velocity_field_new = rotate_velo(random_angle,velocity_field,coordinates_new)
    #print("rotate velocity", velocity_field_new[0])
    value_field = representative_vtu.GetField(field_names[1])[:,0]
    #print("value shape", value_field.shape)

    # plt.figure(figsize=(9,9))
    # plt.rcParams.update({'font.size':18})
    grid_origin = random_select(coordinates_new)
    # print("grid origin : ", grid_origin)
    # plot_grid(grid_origin, grid_width, nx, ny,coordinates_new)
    # plt.plot(coordinates_new[:, 0], coordinates_new[:, 1], 'g.', ms=0.3)
    # plt.tight_layout()
    # plt.legend(loc='best',fontsize = 18)

    #variables for interpolation
    block_x_start = get_grid_end_points_new(grid_origin, grid_width)
    zeros_on_mesh = 0
    x_all = np.transpose(coordinates_new[:, 0:nDim])

    #interpolate velocity
    u_mesh = np.zeros((1,nNodes,1))
    u_mesh[:,:,0] = np.transpose(velocity_field[:,0])

    u_value_grid = u2r.simple_interpolate_from_mesh_to_grid(u_mesh,x_all,x_ndgln,ddx,block_x_start,nx,ny,nz,zeros_on_mesh, nEl,nloc,nNodes,nScalar_test,nDim,1)
    snapshot_matrix[:nx*ny*1,iGrid] = u_value_grid.reshape(-1)
    #snapshot_output_data[iGrid, :, :, 0] = u_value_grid.reshape((nx,ny))
    snapshot_velocity[iGrid, :, :, 0] = u_value_grid.reshape((nx,ny))

    v_mesh = np.zeros((1,nNodes,1))
    v_mesh[:,:,0] = np.transpose(velocity_field[:,1])

    v_value_grid = u2r.simple_interpolate_from_mesh_to_grid(v_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_on_mesh, nEl, nloc, nNodes, nScalar_test, nDim, 1)
    snapshot_matrix[nx * ny:nx * ny * 2, iGrid] = v_value_grid.reshape(-1)
    #snapshot_output_data[iGrid, :, :, 1] = v_value_grid.reshape((nx,ny))
    snapshot_velocity[iGrid, :, :, 1] = v_value_grid.reshape((nx,ny))

    info_mesh = np.zeros((1,nNodes,1))
    info_mesh[:,:,0] = np.transpose(value_field)

    info_value_grid = u2r.simple_interpolate_from_mesh_to_grid(info_mesh, x_all, x_ndgln, ddx, block_x_start, nx, ny, nz, zeros_on_mesh, nEl, nloc, nNodes, nScalar_test, nDim, 1)
    snapshot_matrix[nx * ny * 2:nx * ny * 3, iGrid] = info_value_grid.reshape(-1)
    #snapshot_output_data[iGrid, :, :, 2] = info_value_grid.reshape((nx,ny))
    snapshot_info[iGrid, :, :, 0] = info_value_grid.reshape((nx,ny))
# plt.show()





minv = np.amin(snapshot_matrix[:nx*ny*2,:])
maxv = np.amax(snapshot_matrix[:nx*ny*2,:])
vel_scaling = 1/(maxv-minv)
va_scaling = 1e-5

snapshot_matrix[:nx*ny*2,:] = vel_scaling*(snapshot_matrix[:nx*ny*2,:]-minv)
snapshot_matrix[nx*ny*2:nx*ny*3,:] = va_scaling*snapshot_matrix[nx*ny*2:nx*ny*3,:]

snapshot_output_data[:,:,:,:2] = vel_scaling*(snapshot_output_data[:,:,:,:2]-minv)
snapshot_output_data[:,:,:,2] = va_scaling*snapshot_output_data[:,:,:,2]


np.save("snaphsots_field.npy", np.array(snapshot_output_data))
np.save("snaphsots_velocity.npy", np.array(snapshot_velocity))
np.save("snaphsots_info.npy", np.array(snapshot_info))
# ---------------------------------------------------------------------------------------
# apply POD to the snapshots
# some POD truncation settings


cumulative_tol = 0.99
#nPOD = [nTime] # len(nPOD) = nFields
nPOD = [100]
#nPOD = [10] # 100 50 10
bases = []
singular_values = []
print("pod -------------",nPOD)

nrows, ncols = snapshot_matrix.shape

if nrows > ncols:
    SSmatrix = np.dot(snapshot_matrix.T, snapshot_matrix)
else:
    SSmatrix = np.dot(snapshot_matrix, snapshot_matrix.T)
    print('WARNING - CHECK HOW THE BASIS FUNCTIONS ARE CALCULATED WITH THIS METHOD')

print('SSmatrix', SSmatrix.shape)
eigvalues, v = np.linalg.eigh(SSmatrix)
print("v:",v)
eigvalues =  eigvalues[::-1]
# get rid of small negative eigenvalues (there shouldn't be any as the eigenvalues of a real, symmetric
# matrix are non-negative, but sometimes very small negative values do appear)
eigvalues[eigvalues<0] = 0
s_values = np.sqrt(eigvalues)
#print('s values', s_values[0:20])

singular_values.append(s_values)

cumulative_info = np.zeros(len(eigvalues))
for j in range(len(eigvalues)):
    if j==0:
        cumulative_info[j] = eigvalues[j]
    else:
        cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

cumulative_info = cumulative_info / cumulative_info[-1]
nAll = len(eigvalues)


if nPOD[0] == -1:
    # SVD truncation - percentage of information captured or number
    #cumulative_tol = nirom_options.compression.cumulative_tol[iField]
    nPOD_iField = sum(cumulative_info <= cumulative_tol) #tolerance
    nPOD[0] = nPOD_iField
elif nPOD[0] == -2:
    nPOD_iField = nAll
    nPOD[0] = nPOD_iField
else:
    nPOD_iField = nPOD[0]

print("retaining", nPOD_iField, "basis functions of a possible", len(eigvalues))


basis_functions = np.zeros((nx*ny*nz*3,nPOD_iField)) # nDim should be nScalar?
for j in reversed(range(nAll-nPOD_iField,nAll)):
    Av = np.dot(snapshot_matrix,v[:,j])
    basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)

bases.append(basis_functions)

write_sing_values(singular_values)
print("basis functions shape:", basis_functions.shape)
print("singular values:", singular_values)
print("basis functions:", bases)
np.save("basis_function.npy",np.array(basis_functions))



nStar = 5000
snapshot_matrix_recon = np.zeros((nx*ny*3, nStar*5))
snapshot_matrix_recon_velocity = np.zeros((nx*ny*2, nStar*5))
for iStar in range (nStar):
    random_angle = random.randint(0, 360)
    coordinates_new = rotate(random_angle, coordinates_org)
    velocity_field = representative_vtu.GetField(field_names[0])[:, :nDim]
    velocity_field_new = rotate_velo(random_angle, velocity_field, coordinates_new)
    value_field = representative_vtu.GetField(field_names[1])[:, 0]
    grid_origin = random_select(coordinates_new)

    # plt.figure(figsize=(9,9))
    # plt.rcParams.update({'font.size':18})
    # print("grid origin : ", grid_origin)
    # plot_grid(grid_origin, grid_width, nx, ny,coordinates_new)
    # plt.plot(coordinates_new[:, 0], coordinates_new[:, 1], 'g.', ms=0.3)
    # plt.tight_layout()
    # plt.legend(loc='best',fontsize = 18)

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
#plt.show()

reconstruction_data = []
basis = bases[0]

minv = np.amin(snapshot_matrix_recon[:nx*ny*2,:])
maxv = np.amax(snapshot_matrix_recon[:nx*ny*2,:])
vel_scaling = 1/(maxv-minv)
va_scaling = 1e-5

snapshot_matrix_recon[:nx*ny*2,:] = vel_scaling*(snapshot_matrix_recon[:nx*ny*2,:]-minv)
snapshot_matrix_recon[nx*ny*2:nx*ny*3,:] = va_scaling*snapshot_matrix_recon[nx*ny*2:nx*ny*3,:]


# alphas = np.zeros((snapshot_matrix_recon.shape[1], nPOD_iField))
# for i in range(snapshot_matrix_recon.shape[1]):
#     alphas[i] = np.dot(basis_functions.T, snapshot_matrix_recon[:,i])
#
# np.save("POD_coe.npy", np.array(alphas))
# print("pod shape",alphas.shape)
#
# reconstructions = np.zeros(snapshot_matrix_recon.shape)
# mse_errors = []
# for i in range(snapshot_matrix_recon.shape[1]):
#     reconstructions[:,i] = np.dot(basis_functions , alphas[i])
#     #mse_errors.append(mean_squared_error(snapshot_matrix_recon[:,i],reconstructions[:,i]))
#
# reconstruction_per_grid = np.zeros((5*nStar, nx, ny, 3))
# for iIter in range(snapshot_matrix_recon.shape[1]):
#     reconstruction_per_grid[iIter, :, :, 0] = reconstructions[:nx * ny * 1, iIter].reshape(41, 41)
#     reconstruction_per_grid[iIter, :, :, 1] = reconstructions[nx * ny:nx * ny * 2, iIter].reshape(41, 41)
#     reconstruction_per_grid[iIter, :, :, 2] = reconstructions[nx * ny * 2:nx * ny * 3, iIter].reshape(41, 41)




reconstruction_per_org = np.zeros((5*nStar, nx, ny, 3))
for iIter in range(5*nStar):

    reconstruction_per_org[iIter, :, :, 0] = snapshot_matrix_recon[:nx * ny * 1, iIter].reshape(48, 48)
    reconstruction_per_org[iIter, :, :, 1] = snapshot_matrix_recon[nx * ny:nx * ny * 2, iIter].reshape(48, 48)
    reconstruction_per_org[iIter, :, :, 2] = snapshot_matrix_recon[nx * ny * 2:nx * ny * 3, iIter].reshape(48, 48)


#print(reconstruction_per_grid.shape)
np.save("snaphsots_res_org_48.npy", np.array(reconstruction_per_org))
#np.save("snaphsots_res.npy", np.array(reconstruction_per_grid))
