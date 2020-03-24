#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>


vtkm::cont::DataSet lissajous(int n)
{
    using std::cos;
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids;
    int samples = 512;
    for (size_t i = 0; i < samples; ++i)
    {
        double t = double(i)/double(samples);
        vtkm::Id pid = dsb.AddPoint({-cos(M_PI*(n+1)*t), -cos(M_PI*n*t), -cos(M_PI*(n+2)*t)});
        ids.push_back(pid);
    }
    dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
    return dsb.Create();
}

void random_rectilinear(int samples)
{
    // Padua points, first family:
    // https://en.wikipedia.org/wiki/Padua_points
    // 
    std::random_device rd;
    std::uniform_real_distribution<double> dis(0.0, 1.0/samples);
    std::vector<double> x(samples);
    std::vector<double> y(samples);
    std::vector<double> z(samples);
    for (size_t i = 0; i < samples; ++i)
    {
        x[i] = dis(rd);
        y[i] = dis(rd);
        z[i] = 0;
    }
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    vtkm::cont::DataSetBuilderRectilinear dsb;
    vtkm::cont::DataSet ds = dsb.Create(x, y, z);
    vtkm::io::writer::VTKDataSetWriter writer("random_rectilinear.vtk");
    writer.WriteDataSet(ds);
}

int main(int argc, char** argv)
{
    int samples = 8;
    random_rectilinear(samples);

    // Now the Lissajous curve:
    vtkm::io::writer::VTKDataSetWriter writer("lissajous.vtk");
    writer.WriteDataSet(lissajous(7));

}
