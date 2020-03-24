#include <chrono>
#include <complex>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

std::complex<double> jacobi_theta(std::complex<double> const & z)
{
    const double pi = 4*std::atan(1.0);
    std::complex<double> i{0.0, 1.0};
    std::complex<double> q = 0.1*std::exp(0.1*i*pi);
    std::complex<double> u = i*pi*z;
    std::complex<double> theta = std::sin(u);
    int64_t n = 1;
    std::complex<double> term = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()} ;
    while (std::abs(term) > std::numeric_limits<double>::epsilon())
    {
        term = std::pow(q, n*(n+1))*std::sin((2*n+1.0)*u);
        if (n&1)
        {
            theta -= term;
        }
        else
        {
            theta += term;          
        }
        ++n;
        if (std::isnan(theta.real()))
        {
            std::cerr << "Theta function evaluation diverged.\n";
            std::cerr << "term = " << term << ", theta = " << theta << ", n = " << n << ", n*n = " << n*n << ", z = " << z << "\n";
            throw 1;
        }
    }
    return 2.0*std::pow(q, 0.25)*theta;
}


vtkm::cont::DataSet jacobi_theta_image_dataset(double s_min, double s_max, double t_min, double t_max, int s_samples)
{
    vtkm::cont::DataSetBuilderUniform dsb;
    vtkm::Id2 dims(s_samples, s_samples);
    vtkm::Vec2f_64 origin(t_min, s_min);
    double ds = (s_max - s_min)/double(dims[0] - 1);
    double dt = (t_max - t_min)/double(dims[1] - 1);
    
    vtkm::Vec2f_64 spacing(dt, ds);
    
    vtkm::cont::DataSet dataSet = dsb.Create(dims, origin, spacing);
    vtkm::cont::DataSetFieldAdd dsf;
    vtkm::Id nVerts = s_samples*s_samples;
    std::vector<double> r(nVerts);
    std::vector<double> theta(nVerts);

    vtkm::Id idx = 0;
    for (vtkm::Id y = 0; y < dims[0]; ++y)
    {
        for (vtkm::Id x = 0; x < dims[1]; ++x)
        {
          double wy = s_min + y*ds;
          double wx = t_min + x*dt;
          
          std::complex<double> z = jacobi_theta(std::complex<double>(wx, wy));
          r[idx] = std::sqrt(std::norm(z));
          theta[idx] = std::arg(z);
          idx++;
        }
    } 
    
    dsf.AddPointField(dataSet, "r", r.data(), r.size());
    dsf.AddPointField(dataSet, "θ", theta.data(), theta.size());
    return dataSet;
}


int main(int argc, char* argv[])
{
    vtkm::cont::DataSet input = jacobi_theta_image_dataset(-3, 3, -3, 3, 8);
    vtkm::io::writer::VTKDataSetWriter writer("jacobi_theta_image.vtk");
    writer.WriteDataSet(input);
    return 0;
}
