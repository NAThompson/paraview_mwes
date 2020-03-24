/*
 * The Gray-Scott model of reaction diffusion.
 * See:
 *    Pearson, John E. "Complex patterns in a simple system." Science 261.5118 (1993): 189-192.
 *
 * This code attempts to reproduce the notation in that paper down to the casing.
 */
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

constexpr const struct parameters {
    double Du = 2e-5;
    double Dv = 1e-5;
    double x_max = 2.5;
    double y_max = 2.5;
    double noise_level = 0.01;
    double k = 0.05;
    double F = 0.01;
    double t_max = 10.0;
    int64_t n = 256;
} gs_params;


void initialize_UV(Eigen::MatrixXd & U, Eigen::MatrixXd & V)
{
    // "Initially, the entire system was placed in the trivial state (U = 1,V = 0)."
    assert(U.rows() == V.rows());
    assert(U.cols() == V.cols());
    assert(U.rows() == U.cols());
    int64_t n = U.rows();
    for (int64_t i = 0; i < n; ++i)
    {
        for (int64_t j = 0; j < n; ++j)
        {
            U(i,j) = 1;
            V(i,j) = 0;
        }
    }
    //  "The 20 by 20 mesh point area located symmetrically about the center of the grid was then perturbed to (U = 1/2,V = 1/4)."
    //  We'll use a size 20/256 to spot check:
    int64_t l = 10*n/256.0;
    for (int64_t i = n/2 - l; i < n/2 + l; ++i)
    {
        for (int64_t j = n/2 - l; j < n/2 + l; ++j)
        {
            U(i,j) = 0.5;
            V(i,j) = 0.25;
        }
    }
    // "These conditions were then perturbed with + 1% random noise in order to break the square symmetry"
    std::random_device rd;
    std::uniform_real_distribution dis(-1.0, 1.0);
    for (int64_t i = 0; i < n; ++i)
    {
        for (int64_t j = 0; j < n; ++j)
        {
            U(i,j) += gs_params.noise_level*dis(rd);
            V(i,j) += gs_params.noise_level*dis(rd);
        }
    }
}

void write_data(Eigen::MatrixXd const & U, Eigen::MatrixXd const & V, int64_t k)
{
    std::cout << "Writing the data at step " << k << "\n";
    vtkm::cont::DataSetBuilderUniform dsb;
    vtkm::Id2 dims(gs_params.n, gs_params.n);
    vtkm::Vec2f_64 origin(0, 0);
    double dx = gs_params.x_max/gs_params.n;
    vtkm::Vec2f_64 spacing(dx, dx);    
    vtkm::cont::DataSet dataSet = dsb.Create(dims, origin, spacing);
    vtkm::cont::DataSetFieldAdd dsf;
    dsf.AddPointField(dataSet, "U", U.data(), U.size());
    dsf.AddPointField(dataSet, "V", V.data(), V.size());
    vtkm::io::writer::VTKDataSetWriter writer("gray_scott_" + std::to_string(k) + ".vtk");
    writer.WriteDataSet(dataSet);
}

void solve()
{
    int64_t n = gs_params.n;
    double dx = gs_params.x_max/n;
    // Stability condition: dt < dx*dx/2:
    double dt = dx*dx/4;
    // But that takes forever to run.
    // Instead, we have to live dangerously: (I suspect this is why D_u and D_v are taken so small in the reference)
    dt = dx/5;
    
    Eigen::MatrixXd U0(n,n);
    Eigen::MatrixXd U1(n,n);
    Eigen::MatrixXd V0(n,n);
    Eigen::MatrixXd V1(n,n);
    for (int64_t i = 0; i < n; ++i) {
        for (int64_t j = 0; j < n; ++j) {
            U1(i,j) = std::numeric_limits<double>::quiet_NaN();
            V1(i,j) = std::numeric_limits<double>::quiet_NaN();
        }
    }
    initialize_UV(U0, V0);
    int64_t k = 0;
    int64_t step_max = gs_params.t_max/dt;
    step_max = 2000000;
    double rhsU, rhsV;
    while (k < step_max)
    {
        // "The boundary conditions are periodic."
        // 5 point Laplacian stencil under periodic BCs.
        // This gives 4 corners (i,j) = (0,0), (0, n-1), (n-1, 0), (n-1, n-1) as special cases:
        // i = 0, j = 0: 
        rhsU = gs_params.Du*(U0(1,0) + U0(0,1) - 4*U0(0,0) + U0(n-1, 0) + U0(0, n-1))/(dx*dx) - U0(0,0)*V0(0,0)*V0(0,0) + gs_params.F*(1-U0(0,0));
        rhsV = gs_params.Dv*(V0(1,0) + V0(0,1) - 4*V0(0,0) + V0(n-1, 0) + V0(0, n-1))/(dx*dx) + U0(0,0)*V0(0,0)*V0(0,0) - (gs_params.F+gs_params.k)*V0(0,0);
        U1(0, 0) = U0(0, 0) + dt*rhsU;
        V1(0, 0) = V0(0, 0) + dt*rhsV;

        // i = 0, j = n-1:
        rhsU = gs_params.Du*(U0(1,n-1) + U0(0,0) - 4*U0(0,n-1) + U0(n-1, n-1) + U0(0, n-2))/(dx*dx) - U0(0,n-1)*V0(0,n-1)*V0(0,n-1) + gs_params.F*(1-U0(0,n-1));
        rhsV = gs_params.Dv*(V0(1,n-1) + V0(0,0) - 4*V0(0,n-1) + V0(n-1, n-1) + V0(0, n-2))/(dx*dx) + U0(0,n-1)*V0(0,n-1)*V0(0,n-1) - (gs_params.F+gs_params.k)*V0(0,n-1);
        U1(0,n-1) = U0(0,n-1) + dt*rhsU;
        V1(0,n-1) = V0(0,n-1) + dt*rhsV;

        // i = n-1, j = 0:
        rhsU = gs_params.Du*(U0(0,0) + U0(n-1,1) - 4*U0(n-1,0) + U0(n-2, 0) + U0(n-1, n-1))/(dx*dx) - U0(n-1,0)*V0(n-1,0)*V0(n-1,0) + gs_params.F*(1-U0(n-1,0));
        rhsV = gs_params.Dv*(V0(0,0) + V0(n-1,1) - 4*V0(n-1,0) + V0(n-2, 0) + V0(n-1, n-1))/(dx*dx) + U0(n-1,0)*V0(n-1,0)*V0(n-1,0) - (gs_params.F+gs_params.k)*V0(n-1,0);
        U1(n-1, 0) = U0(n-1, 0) + dt*rhsU;
        V1(n-1, 0) = V0(n-1, 0) + dt*rhsV;

        // i = n-1, j = n-1:
        rhsU = gs_params.Du*(U0(0,n-1) + U0(n-1,0) - 4*U0(n-1,n-1) + U0(n-2, n-1) + U0(n-1, n-2))/(dx*dx) - U0(n-1,n-1)*V0(n-1,n-1)*V0(n-1,n-1) + gs_params.F*(1-U0(n-1,n-1));
        rhsV = gs_params.Dv*(V0(0,n-1) + V0(n-1,0) - 4*V0(n-1,n-1) + V0(n-2, n-1) + V0(n-1, n-2))/(dx*dx) + U0(n-1,n-1)*V0(n-1,n-1)*V0(n-1,n-1) - (gs_params.F+gs_params.k)*V0(n-1,n-1);
        U1(n-1, n-1) = U0(n-1, n-1) + dt*rhsU;
        V1(n-1, n-1) = V0(n-1, n-1) + dt*rhsV;

        // And 4 sides as special cases:
        // i = 0:
        for (int64_t j = 1; j < n - 1; ++j)
        {
            rhsU = gs_params.Du*(U0(1,j) + U0(0,j+1) - 4*U0(0,j) + U0(n-1, j) + U0(0, j-1))/(dx*dx) - U0(0,j)*V0(0,j)*V0(0,j) + gs_params.F*(1-U0(0,j));
            rhsV = gs_params.Dv*(V0(1,j) + V0(0,j+1) - 4*V0(0,j) + V0(n-1, j) + V0(0, j-1))/(dx*dx) + U0(0,j)*V0(0,j)*V0(0,j) - (gs_params.F+gs_params.k)*V0(0,j);
            U1(0,j) = U0(0,j) + dt*rhsU;
            V1(0,j) = V0(0,j) + dt*rhsV;
        }
        // i = n-1:
        for (int64_t j = 1; j < n - 1; ++j)
        {
            rhsU = gs_params.Du*(U0(0,j) + U0(n-1,j+1) - 4*U0(n-1,j) + U0(n-2, j) + U0(n-1, j-1))/(dx*dx) - U0(n-1,j)*V0(n-1,j)*V0(n-1,j) + gs_params.F*(1-U0(n-1,j));
            rhsV = gs_params.Dv*(V0(0,j) + V0(n-1,j+1) - 4*V0(n-1,j) + V0(n-2, j) + V0(n-1, j-1))/(dx*dx) + U0(n-1,j)*V0(n-1,j)*V0(n-1,j) - (gs_params.F+gs_params.k)*V0(n-1,j);
            U1(n-1,j) = U0(n-1,j) + dt*rhsU;
            V1(n-1,j) = V0(n-1,j) + dt*rhsV;
        }
        // j = 0:
        for (int64_t i = 1; i < n - 1; ++i)
        {
            rhsU = gs_params.Du*(U0(i+1,0) + U0(i,1) - 4*U0(i,0) + U0(i-1, 0) + U0(i, n-1))/(dx*dx) - U0(i,0)*V0(i,0)*V0(i,0) + gs_params.F*(1-U0(i,0));
            rhsV = gs_params.Dv*(V0(i+1,0) + V0(i,1) - 4*V0(i,0) + V0(i-1, 0) + V0(i, n-1))/(dx*dx) + U0(i,0)*V0(i,0)*V0(i,0) - (gs_params.F+gs_params.k)*V0(i,0);
            U1(i,0) = U0(i,0) + dt*rhsU;
            V1(i,0) = V0(i,0) + dt*rhsV;
        }
        // j = n-1:
        for (int64_t i = 1; i < n - 1; ++i)
        {
            rhsU = gs_params.Du*(U0(i+1,n-1) + U0(i,0) - 4*U0(i,n-1) + U0(i-1, n-1) + U0(i, n-2))/(dx*dx) - U0(i,n-1)*V0(i,n-1)*V0(i,n-1) + gs_params.F*(1-U0(i,n-1));
            rhsV = gs_params.Dv*(V0(i+1,n-1) + V0(i,0) - 4*V0(i,n-1) + V0(i-1, n-1) + V0(i, n-2))/(dx*dx) + U0(i,n-1)*V0(i,n-1)*V0(i,n-1) - (gs_params.F+gs_params.k)*V0(i,n-1);
            U1(i,n-1) = U0(i,n-1) + dt*rhsU;
            V1(i,n-1) = V0(i,n-1) + dt*rhsV;                
        }
        
        // And now the interior:
        // This is a great target for an openmp parallel for, but I don't want to make this have any more dependencies than it already has!
        //#pragma omp parallel for
        for (int64_t i = 1; i < n - 1; ++i)
        {
            for (int64_t j = 1; j < n - 1; ++j)
            {
                rhsU = gs_params.Du*(U0(i+1,j) + U0(i,j+1) - 4*U0(i,j) + U0(i-1, j) + U0(i, j-1))/(dx*dx) - U0(i,j)*V0(i,j)*V0(i,j) + gs_params.F*(1-U0(i,j));
                rhsV = gs_params.Dv*(V0(i+1,j) + V0(i,j+1) - 4*V0(i,j) + V0(i-1, j) + V0(i, j-1))/(dx*dx) + U0(i,j)*V0(i,j)*V0(i,j) - (gs_params.F+gs_params.k)*V0(i,j);
                U1(i,j) = U0(i,j) + dt*rhsU;
                V1(i,j) = V0(i,j) + dt*rhsV;                
            }
        }

        int64_t step_skip = 5000;
        if (k % step_skip == 0)
        {
            write_data(U0, V0, k);
        }
        U1.swap(U0);
        V1.swap(V0);
        ++k;
    }
    write_data(U0, V0, k);
}

int main(int argc, char** argv)
{
    solve();
}
