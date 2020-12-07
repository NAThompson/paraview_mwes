#include <complex>
#include <cmath>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/VTKDataSetWriter.h>

template<typename Real>
vtkm::Vec<Real, 3> breather_surface(Real u, Real v)
{
    using std::sqrt;
    Real a = Real(2)/Real(5);
    Real coshau = std::cosh(a*u);
    Real sinhau = std::sinh(a*u);
    Real t = 1-a*a;
    Real sintv = std::sin(sqrt(t)*v);
    Real sinv = std::sin(v);
    Real cosv = std::cos(v);
    Real costv = std::cos(sqrt(t)*v);

    vtkm::Vec<Real, 3> p;
    Real denom = a*(t*coshau*coshau + a*a*sintv*sintv);
    p[0] = -u + 2*t*coshau*sinhau/denom;
    Real ynum = 2*sqrt(t)*coshau*(-sqrt(t)*cosv*costv - sinv*sintv );
    p[1] = ynum/denom;
    Real znum = 2*sqrt(t)*coshau*(-sqrt(t)*sinv*costv + cosv*sintv);
    p[2] = znum/denom;
    return p;
}

void write_dataset(int64_t u_samples, int64_t v_samples)
{
    double u_max = 14;
    double u_min = -u_max;
    double v_max = 37.4;
    double v_min = -v_max;
    vtkm::cont::DataSetBuilderExplicitIterative dsb;
    std::vector<vtkm::Id> ids(u_samples*v_samples);
    int64_t idx = 0;
    double du = (u_max-u_min)/(u_samples-1);
    double dv = (v_max-v_min)/(u_samples-1);
    for (int64_t i = 0; i < u_samples; ++i)
    {
        double u = u_min + i*du;
        for (int64_t j = 0; j < v_samples; ++j)
        {
            double v = v_min + j*dv;
            vtkm::Vec<double, 3> p = breather_surface(u, v);
            dsb.AddPoint(p);
            ++idx;
        }
    }
    vtkm::cont::DataSet dataSet = dsb.Create();

    vtkm::io::VTKDataSetWriter writer("breather_surface.vtk");
    writer.WriteDataSet(dataSet);
}

int main()
{
    write_dataset(64, 64);
}
