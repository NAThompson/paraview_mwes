## Intro to the Paraview Ecosystem

---

## What is Paraview?

A 3D graphics GUI!

(But fewer dimensions ok too!)

---

## What is the ecosystem?

The set of programming tools built around Paraview to achieve scientific goals: VTK-m, ADIS, ADIOS.

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

Idiom: Take a *source*, and apply a *filter*, export to paper.

---


![](CurvatureWorkflow.mov)

---

![](CurvatureOfCone.png)

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

![](CSVWorkflow.mov)


---

![35%](InterpolatorConvergence.png)

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

![](1dgraph.mov)

---

![45%](KdVState.png)

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

![](kdv_simulation.mov)

---

![](kdv.avi)

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


![](Scalogram.mov)

---

![40%](log_of_scalogram.png)

---

## Warp by Scalar

2D image with a field on it can be represented by a colortable, or by a surface.

Let's see how to create a surface.

---

![](WarpByScalar.mov)

---

![100%](warp_by_scalar.png)

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

![](jacobi_theta_phase.mov)

---

![80%](jacobi_theta_phase.png)

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

![](Fredholm.mov)

---

![35%](fredholm.png)

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

![](EulerSpiral.mov)

---

![45%](euler_spiral.png)

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


![](GrayScottRealTime.mov)

---

![](GrayScott.avi)

---

![](GrayScottWarpedByScalar.avi)

---











