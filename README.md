# TricubicInterpolations.jl
Julia implementation of a tricubic interpolator in three dimensions. The scheme is based on [Lekien and Marsden (2005), "Tricubic interpolation in three dimensions," Int. J. Numer. Meth. Eng. 63, 455](https://doi.org/10.1002/nme.1296).

## Usage
`TricubicInterpolations.jl` is written using pure `Julia` and has no additional dependencies.

Here is a simple example to get you started. We start with
```
using TricubicInterpolations
```
We will consider the following function:
$$f(x, y, z) = - x^3 + x + y^2 - z.$$
The interpolator object `Tricubic` accepts four inputs `(X, Y, Z, F)`, which are the samples of the three independent variables $(x, y, z)$ and the one dependent variable $f$. These can be generated for our function as
```
f(x, y, z) = - x^3 + x + y^2 - z

X = Y = Z = LinRange(-1, 1, 21)
F = [f(x, y, z) for x=X, y=Y, z=Z]
```
Then the interpolator is initialised as
```
t = Tricubic(X, Y, Z, F)
```
The interpolator can be called at a point, say $(0.5, -0.1, 0.3)$, for an estimate of the function
```
t(0.5, -0.1, 0.3)
```
and its derivatives
```
partial_derivative_x(t, 0.5, -0.1, 0.3)
partial_derivative_y(t, 0.5, -0.1, 0.3)
partial_derivative_z(t, 0.5, -0.1, 0.3)
```
