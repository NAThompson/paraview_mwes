#include <iostream>
#include <string>
#include <cmath>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>


// This momentum is a conserved quantity of the numerical method.
// Use it to sanity check the solution:
template<typename Real>
Real momentum(std::vector<Real> const & u)
{
    Real p = 0;
    for (int64_t i = 0; i < u.size(); ++i)
    {
        p += u[i];
    }
    return p;
}

template<typename Real>
void write_step(std::vector<Real> const & u, int64_t i, Real dx)
{
    vtkm::cont::DataSetBuilderUniform dsb;
    Real origin{0};
    
    vtkm::cont::DataSet dataSet = dsb.Create(u.size(), origin, dx);
    vtkm::cont::DataSetFieldAdd dsf;
    dsf.AddPointField(dataSet, "u(x,t)", u.data(), u.size());
    
    vtkm::io::writer::VTKDataSetWriter writer("kdv_" + std::to_string(i) + ".vtk");
    writer.WriteDataSet(dataSet);
}

template<typename Real>
void KdV(int64_t N, Real dt, Real t_max)
{
    using std::ceil;
    using std::cos;
    using std::abs;
    using std::isnan;
    // Sad! Looking forward to std::numeric::pi;
    const Real pi = 4*std::atan(Real(1));

    if (N <= 0)
    {
        throw std::domain_error("N > 0 is required");
    }
    if (dt > 1)
    {
        throw std::domain_error("time step is too big");
    }
    if (dt <= 0)
    {
        throw std::domain_error("dt > 0 is required");
    }

    Real dx = Real(2)/(N);

    int64_t M = ceil(t_max/dt);
    std::cout << "Solving the KdV equation for dx = " << dx << ", dt = " << dt << " and t_max = " << t_max << "\n";
    std::vector<Real> u0(N);
    for (int64_t i = 0; i < N; ++i)
    {
        u0[i] = cos(pi*i*dx);
    }

    write_step<Real>(u0, 0, dx);
    Real original_momentum = momentum(u0);
    std::cout << "j = 0 momentum = " << original_momentum << "\n";

    std::vector<Real> u1(N);
    for (int64_t i = 0; i < N; ++i)
    {
        Real cdt = cos(pi*i*dx)*dt;
        u1[i] = cos(pi*(i*dx - cdt));
    }
    std::cout << "j = 1 momentum = " << momentum(u1) << "\n";

    // Now the update:
    Real k1 = dt/(3*dx);
    Real delta = 0.022;
    Real k2 = delta*delta*dt/(dx*dx*dx);
    Real t1, t2;
    std::vector<Real> u2(N);
    for (int64_t j = 1; j < M - 1; ++j)
    {
        t1 = (u1[1] + u1[0] + u1[N-1])*(u1[1] - u1[N-1]);
        t2 = u1[2]- 2*u1[1] + 2*u1[N-1] - u1[N-2];
        u2[0] = u0[0] - k1*t1 - k2*t2;
        t1 = (u1[2] + u1[1] + u1[0])*(u1[2] - u1[0]);
        t2 = u1[3] - 2*u1[2] + 2*u1[0] - u1[N-1];
        u2[1] = u0[1] - k1*t1 - k2*t2;
        for (int64_t i = 2; i < N-2; ++i)
        {
        t1 = (u1[i+1] + u1[i] + u1[i-1])*(u1[i+1] - u1[i-1]);
        t2 = u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2];
        u2[i] = u0[i] - k1*t1 - k2*t2;
        }
        u2[N-2] = u0[N-2] - k1*(u1[N-1] + u1[N-2] + u1[N-3])*(u1[N-1] - u1[N-3]) - k2*(u1[0] - 2*u1[N-1] + 2*u1[N-3] - u1[N-4]);
        u2[N-1] = u0[N-1] - k1*(u1[0  ] + u1[N-1] + u1[N-2])*(u1[0  ] - u1[N-2]) - k2*(u1[1] - 2*u1[0  ] + 2*u1[N-2] - u1[N-3]);

        Real p = momentum(u2);
        if (abs(p) > sqrt(std::numeric_limits<Real>::epsilon()) || isnan(p))
        {
            std::cout << "Solution diverged at t = " << (j+1)*dt << "\n";
            std::cout << "Momentum = " << p << "\n";
            return;
        }
        if ( (j+1) % 10000*4000 == 0)
        {
            write_step<Real>(u2, j+1, dx);
        }
        u0 = u1;
        u1 = u2;
    }
    Real p = momentum(u2);
    std::cout << "Final momentum = " << p << "\n";
}


int main(int argc, char** argv)
{
    using Real = double;
    int64_t N = 256;
    Real t_max = 5;
    if (argc > 1)
    {
        std::string dx_str = argv[1];
        if (dx_str == "-h")
        {
            std::cout << "Usage: ./KdV.x N t_max, where N is number of spacial gridpoints (dt chosen from dx via stability conditions) and t_max is max simulation time; e.g., ./KdV.x 512 10\n";
            return 0;
        }
        N = std::stoi(dx_str);
    }
    if (argc > 2)
    {
        t_max = Real(std::stod(argv[2]));
    }

    Real dx = Real(1)/N;
    Real dt = 27*dx*dx*dx/4;
    KdV<Real>(N, dt, t_max);
}
