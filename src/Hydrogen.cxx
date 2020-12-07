#include <complex>
#include <cmath>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/constants/constants.hpp>

#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

template<typename Real>
std::complex<Real> hydrogen_wavefunction(int64_t n, int64_t l, int64_t m, Real r, Real theta, Real phi)
{
    using std::pow;
    using std::sqrt;
    using boost::math::factorial;
    using boost::math::spherical_harmonic;
    using boost::math::laguerre;
    Real scale_sq = pow(Real(2)/n, 3)*factorial<Real>(n - l -1)/(factorial<Real>(n-l)*2*n);
    std::complex<Real> Y = spherical_harmonic(l, m, theta, phi);
    Real rho = 2*r/n;
    Real z = exp(-rho/2)*pow(rho, l)*laguerre(n-l-1, 2*l+1, rho);
    return sqrt(scale_sq)*Y*z;
}

void write_dataset(int64_t n, int64_t l, int64_t m, int64_t cbrt_samples)
{
    int64_t samples = std::pow(cbrt_samples, 3);
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids(samples);
    std::vector<double> rsq(samples, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> phase(samples, std::numeric_limits<double>::quiet_NaN());
    double r_max = 3;
    double dr = r_max/cbrt_samples;
    double dtheta = boost::math::constants::pi<double>()/cbrt_samples;
    double dphi = boost::math::constants::two_pi<double>()/cbrt_samples;
    int64_t idx = 0;
    // TODO: Use vktm::cont::CoordinateSystem to do this in spherical coordinates directly.
    for (int64_t i = 0; i < cbrt_samples; ++i)
    {
        double r = (i+1)*dr;
        for (int64_t j = 0; j < cbrt_samples; ++j)
        {
            double theta = j*dtheta;
            for (int64_t k = 0; k < cbrt_samples; ++k)
            {
                double phi = k*dphi;
                std::complex<double> psi = hydrogen_wavefunction(n, l, m, r, theta, phi);
                double x = r*sin(theta)*sin(phi);
                double y = r*sin(theta)*cos(phi);
                double z = r*cos(theta);
                dsb.AddPoint({x, y, z});
                rsq.at(idx) = std::norm(psi);
                phase.at(idx) = std::arg(psi);
                ++idx;
            }
        }
    }
    vtkm::cont::DataSet dataSet = dsb.Create();
    dataSet.AddPointField("r^2", rsq.data(), rsq.size());
    dataSet.AddPointField("phase", phase.data(), phase.size());
    vtkm::io::VTKDataSetWriter writer("hydrogen.vtk");
    writer.WriteDataSet(dataSet);
}

int main()
{
    write_dataset(2, 1, 0, 96);
}
