#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>


int main()
{
    std::random_device rd;
    std::normal_distribution<double> dis(0.0, 1.0);
    std::array<double, 3> center{0.0, 0.0, 0.0};
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    for (size_t i = 0; i < 5000; ++i)
    {
        double x = center[0] + dis(rd);
        double y = center[1] + dis(rd);
        double z = center[2] + dis(rd);
        dsb.AddPoint({x, y, z});
    }
    
    vtkm::cont::DataSet dataSet = dsb.Create();
    vtkm::io::VTKDataSetWriter writer("scattered.vtk");
    writer.WriteDataSet(dataSet);
}
