# Compute the sum of all tetrahedra volumes of the input
pdi = self.GetInput()
pdo = self.GetOutput()

newData = vtk.vtkDoubleArray()
newData.SetName("Total Volume")
newData.SetNumberOfComponents(1)

numTets = pdi.GetNumberOfCells()
cellData = pdi.GetCellData()

totalVolume = 0
for i in range(numTets):
       volume = cellData.GetArray("Volume").GetValue(i)
       totalVolume = totalVolume + volume

print("Total volume of ", numTets, " input tets is ", totalVolume)
newData.InsertNextValue(totalVolume)
pdo.GetCellData().AddArray(newData)