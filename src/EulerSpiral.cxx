#include <iostream>
#include <string>
#include <cmath>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/constants/constants.hpp>


void write_1d_unstructured(std::vector<std::array<double, 2>> const & spiral_data)
{
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids;
    for (size_t i = 0; i < spiral_data.size(); ++i)
    {
        auto point = spiral_data[i];
        vtkm::Id pid = dsb.AddPoint({point[0], point[1], 0.0});
        ids.push_back(pid);
    }
    dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
    vtkm::cont::DataSet dataSet = dsb.Create();
    vtkm::io::writer::VTKDataSetWriter writer("euler_spiral.vtk");
    writer.WriteDataSet(dataSet);
}


void build_euler_spiral(int samples, double t_max)
{
    std::vector<std::array<double, 2>> spiral_data(samples);
    double dt = 2*t_max/samples;
    
    auto integrator = boost::math::quadrature::tanh_sinh<double>();
    auto x_coord = [](double s) { return std::cos(s*s); };
    auto y_coord = [](double s) { return std::sin(s*s); };
    for (size_t i = 0; i < samples; ++i)
    {
        double t = -t_max + i*dt;
        if (t == 0)
        {
            spiral_data[i] = {0, 0};   
        }
        else
        {
            double x = integrator.integrate(x_coord, 0.0, std::abs(t));
            double y = integrator.integrate(y_coord, 0.0, std::abs(t));
            if (t < 0)
            {
                spiral_data[i] = {-x, -y};   
            }
            else
            {
                spiral_data[i] = {x, y};
            }
        }
    }
    write_1d_unstructured(spiral_data);
}

int main(int argc, char** argv)
{
    int samples = 1024;
    double t_max = 9;
    build_euler_spiral(samples, t_max);
}
