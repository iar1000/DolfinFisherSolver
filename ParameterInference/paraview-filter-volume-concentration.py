import numpy as np
from paraview.vtk.numpy_interface import dataset_adapter as dsa
 
# Script to get total volume and total concentration over all timesteps
# Pipeline: Input -> Threshold -> GroupTimeSteps -> PointToCellData -> ProgrammableFilter with this script
# Output: Outputs the #cells, the total volume, and total concentration of these cells (average of point values) per time-step to csv file

output_path_from_paraview_bin = "../../../Projects/DolfinFisherSolver Version 2/ParameterInference/"
output_meshname = "lh-plial-3mio"
output_D = 0.13
output_rho = 0.025
output_filename = "paraview-volume-concentration-" + output_meshname + "-" + str(output_D) + "-" + str(output_rho)

input = inputs[0]
print("Found input with point keys ", dsa.WrapDataObject(input).PointData.keys())
print("Found input with cell keys ", dsa.WrapDataObject(input).CellData.keys())

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

        # sum up all cell concentrations
        totalConcentration = sum(block_np.CellData['u'])
        
        # print results
        print("\tNumber cells = ", len(volumeArray))
        print("\tTotal volume = ", totalVolume)
        print("\tTotal concentration = ", totalConcentration) 
        
        print("Writing to csv")
        with open(output_path_from_paraview_bin + output_filename + ".csv", 'a') as csvfile:
           csvfile.write("{},{},{},{},{},{}\n".format(i,output_meshname,output_D, output_rho,totalVolume,totalConcentration))

    print("finnished")
else:
    print("not vtkCompositDataSet")
