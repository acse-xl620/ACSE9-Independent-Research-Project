import numpy as np
import vtk, vtktools
import os
import u2r
import matplotlib.pyplot as plt

def get_clean_vtu_file(filename):
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

# find a snapshot that has the same mesh / building configuration as the results you want to put in the vtu (ie. you have to have merged your results from the central grids of the star-shaped domains to one big grid and then mapped from this grid to the unstructured mesh using Chris's fortran subroutines. If your large grid is smaller than the origianl domain then you will have to set an option so that it gives you zeros for nodes on the unstructured mesh that are outside your large grid I think. )

snapshots_filename = './time_data/Flow_past_buildings_20.vtu' # any snapshot with the same geometry / building config as the results you want to write to vtu
results_folder = 'prediction_vtus' # folder for results

# transfer information about the mesh/geometry from snapshots_filename to a new vtu data structure called clean_vtu
clean_vtu = get_clean_vtu_file(snapshots_filename)

nNodes = clean_vtu.ugrid.GetNumberOfPoints()
#nEl = clean_vtu.ugrid.GetNumberOfCells()

velocity = np.zeros((nNodes,3)) #paraview always writes 3D fields, even for 2D
buildings = np.zeros((nNodes,1)) # for scalar fields paraview accepts 1D arrays

# I'm not sure what shape array your predictions are stored in? But ultimately if you have mapedd from your large grid to the unstructured mesh it will probably be (nscalar,nonods,ntime), ie for buildings, (1,nNodes, 1) or (1,nNodes, nTime) if you do all the time levels at once, and (2,nNodes, 1) or (2,nNodes, nTime) for velocity. This will be determined by Chris's routine which maps from the structured grids to the unstructured meshes

# let's assume you have velocity_grid (2, nNodes, nTime) and buildings_grid (1, nNodes, nTime)
nTime = 399 # for example
# set some values to test
velocity_grid = np.zeros((2,nNodes,nTime))
buildings_grid = np.zeros((1,nNodes,nTime))

result_velocity = np.load("./prediction/time_seen_original_velocity.npy")
result_info = np.load("./prediction/time_seen_info.npy")
print(result_velocity.shape)

result_info = result_info.transpose(0,2,1,3)
result_info = result_info[:,:,::-1,:]

result_velocity = result_velocity.transpose(0,2,1,3)
result_velocity = result_velocity[:,:,::-1,:]


plt.imshow(result_velocity[1,:,:,0])
plt.show()

for iTime in range(nTime):
    zeros_beyond_grid = 0  # 0 extrapolate solution; 1 gives zeros for nodes outside grid
    nScalar = 2
    nDim = 2
    ddx = np.array((0.01, 0.01))
    nx = 288
    ny = nx
    nz = 1
    coordinates = clean_vtu.GetLocations()
    x_all = np.transpose(coordinates[:, 0:nDim])
    block_x_start = [0, 0]

    reconstruction_on_velocity = u2r.interpolate_from_grid_to_mesh(result_velocity[:, :, :, iTime],
                                                                             block_x_start, ddx, x_all,
                                                                             zeros_beyond_grid, nScalar, nx, ny, nz,
                                                                             nNodes, nDim, 1)

    velocity_grid[:,:,iTime] = reconstruction_on_velocity[:,:,0]

    nScalar = 1
    nx = 288
    ny = nx
    nz = 1
    reconstruction_on_building_info = u2r.interpolate_from_grid_to_mesh(result_info[:, :, :, iTime],
                                                                             block_x_start, ddx, x_all,
                                                                             zeros_beyond_grid, nScalar, nx, ny, nz,
                                                                             nNodes, nDim, 1)
    buildings_grid[:,:,iTime] = reconstruction_on_building_info[:,:,0]

print(buildings_grid)
cwd = os.getcwd()
if not os.path.isdir(results_folder):
    os.mkdir(results_folder)
os.chdir(results_folder) # will overwrite files in results

for i in range(nTime):

    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtu.ugrid)
    # new_vtu.filename = 'prediction_' + str(i) + '.vtu'
    new_vtu.filename = 'original_' + str(i) + '.vtu'
    # new_vtu.filename = 'CFD_unseen_' + str(i) + '.vtu'
    field = "Velocity_prediction"
    # nDim = 2 # for velocity
    nDim = 1
    # velocity[:,0:nDim] = velocity_grid[0,:,i].reshape((nNodes,nDim),order='F')
    # velocity[:, 1:2] = velocity_grid[1, :, i].reshape((nNodes, nDim), order='F')
    velocity[:, 0] = velocity_grid[0, :, i]
    velocity[:, 1] = velocity_grid[1, :, i]
    new_vtu.AddField(field,velocity)

    field = "Buildings_prediction"
    nDim = 1 # for buildings
    buildings[:,0:nDim] = buildings_grid[:,:,i].reshape((nNodes,nDim),order='F')
    new_vtu.AddField(field,buildings)

    new_vtu.Write()

os.chdir(cwd)

  


