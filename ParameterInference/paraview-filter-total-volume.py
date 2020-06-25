import numpy as np
from paraview.vtk.numpy_interface import dataset_adapter as dsa

# Script to get total volume of a thresholded tumor of all timesteps
# Pipeline: Input -> Threshold -> GroupTimeSteps -> ProgrammableFilter with this script
# Output: Outputs the #cells and the total volume of these cells per time-step

input = inputs[0]
print("Found input with keys ", dsa.WrapDataObject(input).PointData.keys())
if input.IsA("vtkCompositeDataSet"):
    numTimesteps = input.GetNumberOfBlocks()
    
    # process block
    for i in range(numTimesteps):
        print("Process timestep ", i)
        # get block
        block = input.GetBlock(i)                        # type vtkCommonDataModelPython.vtkUnstructuredGrid
        block_np = dsa.WrapDataObject(input.GetBlock(i)) # type paraview.vtk.numpy_interface.dataset_adapter.UnstructuredGrid
        # get volumes
        numTets = block.GetNumberOfCells()
        volumeArray = np.empty( numTets,dtype = np.float64 )
        for j in range( numTets ):
            cell = block.GetCell(j)
            p1 = block.GetPoint( cell.GetPointId(0) )
            p2 = block.GetPoint( cell.GetPointId(1) )
            p3 = block.GetPoint( cell.GetPointId(2) )
            p4 = block.GetPoint( cell.GetPointId(3) )
            volumeArray[j] = abs(vtk.vtkTetra.ComputeVolume(p1 ,p2 ,p3 , p4))
        # total volume
        totalVolume = 0
        for v in volumeArray:
            totalVolume = totalVolume + v
        # print results
        print("\tNumber cells = ", len(volumeArray))
        print("\tTotal volume = ", totalVolume)

    print("finished")
else:
    print("not vtkCompositDataSet")
