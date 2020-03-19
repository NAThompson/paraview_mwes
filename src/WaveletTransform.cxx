#include <chrono>
#include <complex>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <boost/math/quadrature/wavelet_transforms.hpp>


vtkm::cont::DataSet scalogram_log_s_dataset(double log_s_min, double log_s_max, double t_min, double t_max, int s_samples)
{
    vtkm::cont::DataSetBuilderUniform dsb;
    vtkm::Id2 dims(s_samples, s_samples);
    vtkm::Vec2f_64 origin(t_min, log_s_min);
    double dlogs = (log_s_max - log_s_min)/double(dims[0] - 1);
    double dt = (t_max - t_min)/double(dims[1] - 1);
    
    vtkm::Vec2f_64 spacing(dt, dlogs);

    vtkm::cont::DataSet dataSet = dsb.Create(dims, origin, spacing);
    vtkm::cont::DataSetFieldAdd dsf;
    vtkm::Id nVerts = s_samples*s_samples;
    std::vector<double> pointvar(nVerts);
    
    auto f = [](double x)
    {
        return std::cos(10*x)*std::exp(-abs(x));
    };
    auto Wf = boost::math::quadrature::daubechies_wavelet_transform<decltype(f), double, 8>(f);

    vtkm::Id idx = 0;
    auto start =  std::chrono::steady_clock::now();
    for (vtkm::Id y = 0; y < dims[0]; ++y)
    {
      for (vtkm::Id x = 0; x < dims[1]; ++x)
      {
        // x is the index, t is the translation.
        // y is the index, s is the scale.
        double log_s = log_s_min + y*dlogs;
        double t = t_min + x*dt;
        auto z = Wf(std::pow(10, log_s), t);
        pointvar[idx] = z*z;
        idx++;
      }
    } 
    auto end = std::chrono::steady_clock::now();
    double pixels_per_second = double(pointvar.size())/std::chrono::duration_cast<std::chrono::seconds>(end - start).count() ;
    std::cout << "computed " << pointvar.size() << " pixels in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds\n";
    std::cout << " or " << pixels_per_second << " pixels per second.\n";
    
    dsf.AddPointField(dataSet, "|W[f](s,t)|^2", pointvar.data(), pointvar.size());
    return dataSet;
}

int main(int argc, char* argv[])
{
    vtkm::cont::Initialize(argc, argv, vtkm::cont::InitializeOptions::Strict);
    vtkm::cont::DataSet input = scalogram_log_s_dataset(-4, 0, -2.0, 2.0, 512);
    vtkm::io::writer::VTKDataSetWriter writer("scalogram_log_s.vtk");
    writer.WriteDataSet(input);
    return 0;
}