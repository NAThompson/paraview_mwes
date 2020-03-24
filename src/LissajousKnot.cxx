#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>


vtkm::cont::DataSet lissajous_knot(int nx, int ny, int nz, double phix, double phiy)
{
    using std::cos;
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids;
    int samples = 2048;
    for (size_t i = 0; i < samples; ++i)
    {
        double t = 2*M_PI*double(i)/double(samples);
        vtkm::Id pid = dsb.AddPoint({cos(nx*t + phix), cos(ny*t + phiy), cos(nz*t)});
        ids.push_back(pid);
    }
    dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
    return dsb.Create();
}

int main(int argc, char** argv)
{
    vtkm::io::writer::VTKDataSetWriter writer("lissajous_knot.vtk");
    writer.WriteDataSet(lissajous_knot(3,5,7, 0.7, 1.0));
}
