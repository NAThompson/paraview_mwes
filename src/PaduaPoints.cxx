#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/PartitionedDataSet.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>


vtkm::cont::DataSet lissajous(int n)
{
    using std::cos;
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids;
    int samples = 2048;
    for (size_t i = 0; i < samples; ++i)
    {
        double t = M_PI*double(i)/double(samples);
        vtkm::Id pid = dsb.AddPoint({-cos((n+1)*t), -cos(n*t), 0.0});
        ids.push_back(pid);
    }
    dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
    return dsb.Create();
}
vtkm::cont::PartitionedDataSet padua_points(int n)
{
    using std::cos;
    vtkm::cont::DataSetBuilderRectilinear dsb;

    std::vector<double> x1(n+1, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> y1(n+1, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> z1(n+1, std::numeric_limits<double>::quiet_NaN());

    std::vector<double> x2(n+1, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> y2(n+1, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> z2(n+1, std::numeric_limits<double>::quiet_NaN());

    // The Padua points are the union of two Chebyshev grids, see:
    // https://www.math.unipd.it/~marcov/pdf/poster5ecm.pdf
    // C_{n+1} := { cos(jπ/n), j = 0,...n }
    // TODO: Fix this x1 is size n/2.
    for (int i = 0; i < n+1; ++i)
    {
        if (i&1)
        {
            x1[i] = cos(i*M_PI/n);
            y2[i] = cos(i*M_PI/(n+1));
        }
        else
        {
            x2[i] = cos(i*M_PI/n);
            y1[i] = cos(i*M_PI/(n+1));
        }
        z1[i] = 0;
        z2[i] = 0;
    }
    vtkm::cont::DataSet ds1 = dsb.Create(x1, y1, z1);
    vtkm::cont::DataSet ds2 = dsb.Create(x2, y2, z2);
    vtkm::cont::PartitionedDataSet pds;
    pds.AppendPartitions({ds1, ds2});
    return pds;
}



int main(int argc, char** argv)
{
    int n = 12;
    auto pds = padua_points(n);
    for (int i = 0; i < pds.GetNumberOfPartitions(); ++i){
        auto & p = pds.GetPartition(i);
        vtkm::io::writer::VTKDataSetWriter writer("padua_" + std::to_string(i+1) + ".vtk");
        writer.WriteDataSet(p);
    }

    vtkm::io::writer::VTKDataSetWriter writer("lissajous.vtk");
    writer.WriteDataSet(lissajous(n));
}
