#include <iostream>
#include <string>
#include <cmath>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/constants/constants.hpp>

#include <Eigen/Dense>


double rhs(double x)
{
    using boost::math::constants::pi;
    return std::exp(-x*x)*std::sin(5*pi<double>()*x);   
}


void write_1d_unstructured(Eigen::VectorXd const & y, std::vector<double> const & nodes)
{
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids;
    
    std::vector<vtkm::Vec3f_64> coords(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        vtkm::Id pid = dsb.AddPoint({nodes[i], 0.0, 0.0});
        ids.push_back(pid);
    }
    
    dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
    vtkm::cont::DataSet dataSet;
    
    dataSet = dsb.Create();
    vtkm::cont::DataSetFieldAdd dsf;
    dsf.AddPointField(dataSet, "y", y.data(), y.size());
    vtkm::io::writer::VTKDataSetWriter writer("fredholm.vtk");
    writer.WriteDataSet(dataSet);
}

void fredholm()
{
    // http://eqworld.ipmnet.ru/en/solutions/ie/ie0401.pdf
    // This code can be verified by taking the limit lambda -> 0, at which point y -> f.
    double lambda = 2.5;
    // take a = -1, b = 1.
    auto gauss_data = boost::math::quadrature::gauss<double, 256>();
    std::vector<double> nodes;
    std::vector<double> weights;
    if (gauss_data.abscissa()[0] == 0)
    {   
        nodes.resize(2*gauss_data.abscissa().size() - 1);
        weights.resize(nodes.size());
        int n = gauss_data.abscissa().size() - 1;
        for(int i = 1; i < gauss_data.abscissa().size(); ++i)
        {
            nodes[n + i] = gauss_data.abscissa()[i];
            weights[n+i] = gauss_data.weights()[i];
            nodes[n - i] = -gauss_data.abscissa()[i];
            weights[n-i] = gauss_data.weights()[i];
        }
    }
    else
    {
        nodes.resize(2*gauss_data.abscissa().size());
        weights.resize(nodes.size());
        int n = gauss_data.abscissa().size();
        for(int i = 0; i < gauss_data.abscissa().size(); ++i)
        {
            nodes[n + i] = gauss_data.abscissa()[i];
            weights[n+i] = gauss_data.weights()[i];
            nodes[n -1 - i] = -gauss_data.abscissa()[i];
            weights[n - 1 - i] = gauss_data.weights()[i];
        }
    }
    
    Eigen::VectorXd f(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        f[i] = rhs(nodes[i]);
    }
    
    Eigen::MatrixXd M(nodes.size(), nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        for (size_t j = 0; j < nodes.size(); ++j)
        {
            if (i == j)
            {
                M(i,j) = 1;   
            }
            else
            {
                M(i,j) = -lambda*(nodes[i] - nodes[j])*weights[j];   
            }
        }
    }
    Eigen::VectorXd y = M.fullPivLu().solve(f);
    
    double relative_error = (M*y - f).norm() / f.norm(); // norm() is L2 norm
    std::cout << "The relative error is " << relative_error << std::endl;
    // The values of the y vector gives the solution at Gaussian quadrature nodes.
    
    // To validate graphically, take lambda -> 0 and replace y by f:
    write_1d_unstructured(y, nodes);
}

int main(int argc, char** argv)
{
    fredholm();
}
