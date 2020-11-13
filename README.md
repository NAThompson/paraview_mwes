## Intro to the Paraview Ecosystem

---

## What is Paraview?

A 3D graphics GUI!

(But fewer dimensions ok too!)

---

## What is the ecosystem?

The set of programming tools built around Paraview to achieve scientific goals: VTK-m, Fides, ADIOS.

---

## My two cents:

Paraview is the crown jewel of the DOE software ecosystem.

It shows we can meet our unique challenges with performant and usable software.

---

## How do I get it?

[paraview.org/download](https://www.paraview.org/download/)

---

## Source build?

Easy (on Ubuntu)

```bash
$ sudo apt-get install qt5-default libqt5x11extras5-dev \
  libqt5help5 qttools5-dev qtxmlpatterns5-dev-tools libqt5svg5-dev
$ git clone --recursive https://gitlab.kitware.com/paraview/paraview.git
$ cd paraview; mkdir build
paraview/build$ cmake ../ -G Ninja
paraview/build$ ninja
paraview/build$ ./bin/paraview
```


---

## MWE: Compute Gaussian curvature

Idiom: Take a *source*, and apply a *filter*, create a render.

---


![](figures/CurvatureWorkflow.mov)

---

![](figures/CurvatureOfCone.png)

---

## That was cheating.

Yes, we used canned data packaged with the software.

Let's upload some real data.

---

## Uploading data

Easiest thing: A `.csv`:

```
➜  paraview_tutorial$ cat daubechies_5_scaling_convergence.csv | more
r, linear, quadratic_b_spline, cubic_b_spline, quintic_b_spline, cubic_hermite, pchip, makima, quintic_hermite, septic_hermite
2, 0.045726192045372316, 0.024828473821310204, 0.023277651982070491, 0.020662363191588429, 0.003992590683903730, 0.022484271673629985, 0.029204046982060805, 0.010163942298858419, 0.018448258714584442
3, 0.016022890862746997, 0.005755604628244537, 0.005230814883100732, 0.004517849453574918, 0.001041719154716125, 0.012099678303480466, 0.007022457175191787, 0.002578238856418169, 0.004660243256272567
4, 0.004326423048900463, 0.001479662333354892, 0.001403577051129257, 0.001258471400971933, 0.000264701848608107, 0.003008722563990818, 0.002047077677412190, 0.000660995684770599, 0.001198057778804051
5, 0.001309431207276557, 0.000373346431293220, 0.000339254208968076, 0.000291630399409548, 0.000067783984772363, 0.000465714198778056, 0.000491306619344600, 0.000168562863492538, 0.000305161839582041
```

---

![](figures/CSVWorkflow.mov)


---

![35%](figures/InterpolatorConvergence.png)

---

## Plot a 1D graph

```
➜  paraview_tutorial$ cat data/1dgraph.vtk
# vtk DataFile Version 3.0
vtk output
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS 96 1 1 POINTS 96 float
0 0 0
0.0208333 0 0
0.0416667 0 0
...
1.97917 0 0
POINT_DATA 96
SCALARS u double 1
LOOKUP_TABLE default
0.84768
0.787256
...
0.575502
```

---

![](figures/1dgraph.mov)

---

![45%](figures/KdVState.png)

---

## How about a graph that evolves in time?

Make a series of `.vtk` files, and then label them with an increasing index `foo_1.vtk`, `foo_3.vtk`, `foo_5.vtk`.

Then Paraview knows that this is a dataset that evolves in time.

---

## Korteweg-de Vries equation

$$
\partial_t u + \delta^2\partial_x^3 u  + u \partial_x u = 0, \quad u(x, 0) = \cos(\pi x)
$$

Let's use this as an example of how to watch the evolution of a simulation in Paraview.

---

## Discretize KdV via [Zabusky & Krustal](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.15.240)

$$
u_{i}^{j+1} = u_{i}^{j-1}  - \frac{\Delta t}{3 \Delta x}\left(u_{i+1}^{j} + u_{i}^{j} + u_{i-1}^{j} \right)(u_{i+1}^{j} - u_{i-1}^{j})  - \frac{\delta^2 \Delta t}{\Delta x^3} (u_{i+2}^{j} - 2 u_{i+1}^{j}  + 2u_{i-1}^{j} - u_{i-2}^{j})
$$

where $$\delta$$ is an adjustable parameter, here chosen as 0.022.

---

![](figures/kdv_simulation.mov)

---

![](figures/kdv.avi)

---


## Let's level up

Given a wavelet transform

$$
\mathcal{W}[f](s,t) := \frac{1}{\sqrt{|s|}} \int_{-\infty}^{\infty} f(x) \psi\left( \frac{x-t}{s}  \right)^{*} \, \mathrm{d}x
$$

let's compute a *scalogram*, a graph of $$z = |\mathcal{W}[f](s, t)|^2$$.

---

## Need a `.vtk` image file

```
# vtk DataFile Version 3.0
vtk output
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS 2 2 1 POINTS 4 float 
-2 1e-05 0
2 1e-05 0
-2 1 0
2 1 0
POINT_DATA 4
SCALARS scalogram double 1
LOOKUP_TABLE default
6.86438e-06
6.86117e-06
0.000390553
0.000382475
```

---


![](figures/Scalogram.mov)

---

![40%](figures/log_of_scalogram.png)

---

## Warp by Scalar

2D image with a field on it can be represented by a colortable, or by a surface.

Let's see how to create a surface.

---

![](figures/WarpByScalar.mov)

---

![100%](figures/warp_by_scalar.png)

---


## Compound datasets

Let's look at the Jacobi $$\vartheta$$ function

$$
\vartheta_{3}(z, q) := 2q^{1/4}\sum_{n=0}^{\infty} (-1)^nq^{n(n+1)}\sin((2n+1)i\pi z)
$$

We might want to look at magnitude, phase, real part, or imaginary part of this function, for various values of $$q$$.

---

## Quick taste of VTK-m

```cpp
vtkm::cont::DataSet dataSet = ...
vtkm::cont::DataSetFieldAdd dsf;
// ...
dsf.AddPointField(dataSet, "r", r.data(), r.size());
dsf.AddPointField(dataSet, "theta", theta.data(), theta.size());
dsf.AddPointField(dataSet, "Re(z)", real_part.data(), real_part.size());
dsf.AddPointField(dataSet, "Im(z)", imag_part.data(), imag_part.size());
```

---

![](figures/jacobi_theta_phase.mov)

---

![80%](figures/jacobi_theta_phase.png)

---

## Unstructured

Unstructured gets hard in a hurry. Let's solve an easy 1D problem to wrap our heads around it.

---

## Fredholm integral equation of the second kind

$$
y(x) - \lambda \int_{-1}^{1} (x-t)y(t) \, \mathrm{d}t = f(x)
$$

Easy to validate: As $$\lambda \to 0$$, $$y \to f$$.

---

## Numerics

Let $$\{x_k\}$$ be Gaussian quadrature nodes and $$\{w_k\}$$ be quadrature weights. Let $$y_k := y(x_k)$$. Then the discretization is

$$
y_k - \lambda \sum_{j} w_{j}(x_k - x_j)y_j  = f_k
$$

This is a dense linear system $$My = f$$, which gives the values of $$y$$ on an irregular grid (the roots of the Legendre polynomials).

---

## Numerics -> C++

```cpp
double lambda = 2.5;
auto gauss_data = boost::math::quadrature::gauss<double, 256>();
std::vector<double> nodes = ... ;
std::vector<double> weights = ... ;
Eigen::VectorXd f(nodes.size());
for (size_t i = 0; i < nodes.size(); ++i) {
  f[i] = rhs(nodes[i]);
}

Eigen::MatrixXd M(nodes.size(), nodes.size());
for (size_t i = 0; i < nodes.size(); ++i) {
  for (size_t j = 0; j < nodes.size(); ++j) {
    if (i == j) {
      M(i,j) = 1;   
    }
    else {
      M(i,j) = -lambda*(nodes[i] - nodes[j])*weights[j];   
    }
  }
}
// Solution y at Gauss nodes:
Eigen::VectorXd y = M.fullPivLu().solve(f);

```

---


## VTK-m

```cpp
vtkm::cont::DataSetBuilderExplicitIterative dsb;
std::vector<vtkm::Id> ids;
for (size_t i = 0; i < nodes.size(); ++i) {
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
```

---

![](figures/Fredholm.mov)

---

![35%](figures/fredholm.png)

---

## Can't stop the fun

Let's build an Euler spiral:

$$
x(t) := \int_{0}^{t} \cos(s^2) \, \mathrm{d}s \quad y(t) := \int_{0}^{t} \sin(s^2) \, \mathrm{d}s
$$

---

## Numerics

```cpp
#include <boost/math/quadrature/tanh_sinh.hpp>
auto integrator = boost::math::quadrature::tanh_sinh<double>();
auto x_coord = [](double s) { return std::cos(s*s); };
auto y_coord = [](double s) { return std::sin(s*s); };
// ...
double x = integrator.integrate(x_coord, 0.0, t);
double y = integrator.integrate(y_coord, 0.0, t);
```

---

## VTK-m: No field data, just a path:

```cpp
vtkm::cont::DataSetBuilderExplicitIterative dsb;
std::vector<vtkm::Id> ids;
for (size_t i = 0; i < spiral_data.size(); ++i) {
  auto point = spiral_data[i];
  vtkm::Id pid = dsb.AddPoint({point[0], point[1], 0.0});
  ids.push_back(pid);
}
dsb.AddCell(vtkm::CELL_SHAPE_POLY_LINE, ids);
vtkm::cont::DataSet dataSet = dsb.Create();
vtkm::io::writer::VTKDataSetWriter writer("euler_spiral.vtk");
writer.WriteDataSet(dataSet);
```

---

![](figures/EulerSpiral.mov)

---

![45%](figures/euler_spiral.png)

---

## 2D video: The [Gray-Scott](https://arxiv.org/pdf/patt-sol/9304003.pdf) Initial Value Problem

$$
\partial_t U = D_{u}\nabla^2 U - UV^2 + F(1-U)
$$
$$
\partial_t V = D_{v}\nabla^2V + UV^2 - (F+k)V
$$

Subject to

$$U(x,y,0) = U_0(x,y)$$, $$V(x,y,0) = V_0(x,y)$$.

---

## Numerics

$$
\frac{U_{i,j}^{k+1} - U_{i,j}^{k}}{\Delta t} = D_{u} \frac{U_{i+1,j}^{k} + U_{i,j+1}^{k} - 4U_{i,j}^{j} + U_{i-1,j}^{k} + U_{i,j-1}^{k} }{\Delta x^2} - U_{i,j}^{k} (V_{i,j}^{k})^{2} + F(1-U_{i,j}^{k})
$$

$$
\frac{V_{i,j}^{k+1} - V_{i,j}^{k}}{\Delta t} = D_{v} \frac{V_{i+1,j}^{k} + V_{i,j+1}^{k} - 4V_{i,j}^{j} + V_{i-1,j}^{k} + V_{i,j-1}^{k} }{\Delta x^2} + U_{i,j}^{k} (V_{i,j}^{k})^{2} - (F+k)V_{i,j}^{k}
$$



Stability condition $$\Delta t < \Delta x^2/2$$, but that leads to an incredibly expensive simulation.

Perfect use case for *watching* you sim.

---

## Numerics -> C++

```cpp
Eigen::MatrixXd U0(n,n);
Eigen::MatrixXd U1(n,n);
Eigen::MatrixXd V0(n,n);
Eigen::MatrixXd V1(n,n);
// ...
for (int64_t i = 1; i < n - 1; ++i)
{
  for (int64_t j = 1; j < n - 1; ++j)
  {
    rhsU = Du*(U0(i+1,j) + U0(i,j+1) - 4*U0(i,j) + U0(i-1, j) + U0(i, j-1))/(dx*dx) - U0(i,j)*V0(i,j)*V0(i,j) + F*(1-U0(i,j));
    rhsV = Dv*(V0(i+1,j) + V0(i,j+1) - 4*V0(i,j) + V0(i-1, j) + V0(i, j-1))/(dx*dx) + U0(i,j)*V0(i,j)*V0(i,j) - (F+k)*V0(i,j);
    U1(i,j) = U0(i,j) + dt*rhsU;
    V1(i,j) = V0(i,j) + dt*rhsV;
  }
}
```

---

## VTK-m

```cpp
vtkm::cont::DataSetBuilderUniform dsb;
vtkm::Id2 dims(n, n);
vtkm::Vec2f_64 origin(0, 0);
double dx = gs_params.x_max/n;
vtkm::Vec2f_64 spacing(dx, dx);
vtkm::cont::DataSet dataSet = dsb.Create(dims, origin, spacing);
vtkm::cont::DataSetFieldAdd dsf;
dsf.AddPointField(dataSet, "U", U.data(), U.size());
dsf.AddPointField(dataSet, "V", V.data(), V.size());
vtkm::io::writer::VTKDataSetWriter writer("gray_scott_" + std::to_string(step_index) + ".vtk");
writer.WriteDataSet(dataSet);
```

---


![](figures/GrayScottRealTime.mov)

---

![](figures/GrayScott.avi)

---

![](figures/GrayScottWarpedByScalar.avi)

---

## Onto 3D: A Scatterplot

Very simple VTK-m as no connectivity info is needed:

```cpp
std::random_device rd;
std::normal_distribution<double> dis(0.0, 1.0);
std::array<double, 3> center{0.0, 0.0, 0.0};
vtkm::cont::DataSetBuilderExplicitIterative dsb;
for (size_t i = 0; i < 500; ++i)
{
    double x = center[0] + dis(rd);
    double y = center[1] + dis(rd);
    double z = center[2] + dis(rd);
    dsb.AddPoint({x, y, z});
}

vtkm::cont::DataSet dataSet = dsb.Create();
vtkm::io::writer::VTKDataSetWriter writer("scattered.vtk");
writer.WriteDataSet(dataSet);
```

---

![](figures/ScatterPlot.mov)

---

![](figures/ScatterPlot.png)

---

## 3D Volume Rendering

Start easy with data provided by Paraview.


---

![34%](figures/VolumeRenderingWavelet.mov)

---

Further into the VTK-m data model: Rectilinear grids

![75%](figures/random_rectilinear.png)

---

## VTK-m

```cpp
std::random_device rd;
std::uniform_real_distribution<double> dis(0.0, 1.0/samples);
std::vector<double> x(samples);
std::vector<double> y(samples);
std::vector<double> z(samples);
for (size_t i = 0; i < samples; ++i)
{
  x[i] = dis(rd);
  y[i] = dis(rd);
  z[i] = 0;
}
std::sort(x.begin(), x.end());
std::sort(y.begin(), y.end());

vtkm::cont::DataSetBuilderRectilinear dsb;
vtkm::cont::DataSet ds = dsb.Create(x, y, z);
vtkm::io::writer::VTKDataSetWriter writer("random_rectilinear.vtk");
writer.WriteDataSet(ds);
```

---

## Padua Points

The Padua points are the union of two Chebyshev grids.

Chebyshev nodes:

$$
C_{n+1} := \{\cos(j\pi/n), j=0,\ldots n-1\}
$$

---

## Padua points

$$
\mathrm{Pad}_{n} := (C_{n+1}^{\mathrm{odd}} \times C_{n+2}^{\mathrm{even}}) \cup (C_{n+1}^{\mathrm{even}} \times C_{n+2}^{\mathrm{odd}})
$$

All points lie on the Lissajous curve 

$$g(t) := (-\cos((n+1)t), -\cos(nt)), t \in [0, \pi].$$

---

## VTK-m data model

Two `Rectilinear` grids; perfect use case for `PartitionedDataSet`.

---

## VTK-m

```cpp
vtkm::cont::DataSetBuilderRectilinear dsb;

std::vector<double> x1((n+1)/2, std::numeric_limits<double>::quiet_NaN());
std::vector<double> y1((n+2)/2, std::numeric_limits<double>::quiet_NaN());
std::vector<double> z1(1, 0);
      
std::vector<double> x2((n+2)/2, std::numeric_limits<double>::quiet_NaN());
std::vector<double> y2((n+1)/2, std::numeric_limits<double>::quiet_NaN());
std::vector<double> z2(1, 0);
// initialize with Padua points . . .
vtkm::cont::DataSet ds1 = dsb.Create(x1, y1, z1);
vtkm::cont::DataSet ds2 = dsb.Create(x2, y2, z2);
vtkm::cont::PartitionedDataSet pds;
pds.AppendPartitions({ds1, ds2});
for (int i = 0; i < pds.GetNumberOfPartitions(); ++i)
{
  auto & p = pds.GetPartition(i);
  vtkm::io::writer::VTKDataSetWriter writer("padua_" + std::to_string(i+1) + ".vtk");
  writer.WriteDataSet(p);
}
              
```

---

![100%](figures/padua.png)

---

## Field associations: Staggered grids

In fluid dynamics, it is common to store pressure as a cell variable, and velocity at a point variable.

Hence we have *field associations* in VTK-m.

```
vtkm::cont::Field::Association::POINTS
vtkm::cont::Field::Association::CELL_SET
```


---

## Visualizing big data

The `.vtk` files are nice, intuitive ASCII.

But they quickly overwhelm the Paraview parser, so we need a binary format to move forward.

---

## VTK-m rendering

```cpp
vtkm::rendering::Actor actor(dataSet.GetCellSet(),
                             dataSet.GetCoordinateSystem(),
                             dataSet.GetField("U"));
vtkm::rendering::Scene scene;
scene.AddActor(actor);
vtkm::rendering::MapperRayTracer mapper;
vtkm::rendering::CanvasRayTracer canvas(1920, 1080);
vtkm::rendering::View2D view(scene, mapper, canvas);
view.Initialize();
view.Paint();
view.SaveAs("gray_scott_u_" + std::to_string(k) + ".pnm");
```

---

VTK-m rendering is for sanity checks; it isn't a full-featured renderer like Paraview.

In addition, customizing the graphics in C++ is incredibly painful.

It's an option, but let's move on . . .

---

## ADIOS

ADIOS is a high-performance binary file format with Paraview support.

We'll begin by learning how to do a simple KdV simulation.

---

```cpp
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <adios2.h>

double chirp(double t) {
    double phi0 = 0.7;
    double k = 1.2;
    double f = 3.4;
    return std::sin(phi0 + k*t*t + f*t)*std::exp(-t*t/2);
}

int main(int argc, char** argv) {
    size_t n = 256;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    double t0 = -3;
    double tmax = 3;
    double dt = (tmax-t0)/(n-1);
    std::vector<double> chirp_vec(n, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n; ++i) {
        chirp_vec[i] = chirp(t0 + i*dt);
    }

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("myio");
    auto chirp_variable = io.DefineVariable<double>("chirp", {n}, {0}, {n}, adios2::ConstantDims);

    io.DefineAttribute<double>("t0", t0);
    io.DefineAttribute<double>("dt", dt);
    io.DefineAttribute<std::string>("interpretation", "Equispaced");
    adios2::Engine bp_file_writer = io.Open("chirp_graph.bp", adios2::Mode::Write);
    bp_file_writer.Put(chirp_variable, chirp_vec.data());
    bp_file_writer.Close();
}
```

---

## Did we get what we expected?

```
$ bpls chirp_graph.bp -lav
File info:
  of variables:  1
  of attributes: 3
  statistics:    Min / Max 

  double   chirp           {256} = -0.696846 / 0.974067
  double   dt              attr   = 0.0235294
  string   interpretation  attr   = "Equispaced"
  double   t0              attr   = -3
```

---

## NOTES:

- To get a png to render in paraview, turn of "Map Scalars"

---

