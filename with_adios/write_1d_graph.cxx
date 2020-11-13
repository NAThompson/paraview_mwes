#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <adios2.h>

double chirp(double t)
{
    double phi0 = 0.7;
    double k = 1.2;
    double f = 3.4;
    return std::sin(phi0 + k*t*t + f*t)*std::exp(-t*t/2);
}

int main(int argc, char** argv)
{
    size_t n = 256;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    double t0 = -3;
    double tmax = 3;
    double dt = (tmax-t0)/(n-1);
    std::vector<double> chirp_vec(n, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n; ++i)
    {
        double t = t0 + i*dt;
        chirp_vec[i] = chirp(t);
    }

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("myio");
    auto chirp_variable = io.DefineVariable<double>("chirp", {n}, {0}, {n}, adios2::ConstantDims);

    io.DefineAttribute<double>("t0", t0);
    io.DefineAttribute<double>("dt", dt);
    double origin[2] = {0.3, 3.7};
    io.DefineAttribute<double>("origin", origin, 2);
    io.DefineAttribute<std::string>("interpretation", "Equispaced");
    adios2::Engine bp_file_writer = io.Open("chirp_graph.bp", adios2::Mode::Write);
    bp_file_writer.Put(chirp_variable, chirp_vec.data());
    bp_file_writer.Close();
}
