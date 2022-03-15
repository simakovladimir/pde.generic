# Generic PDE Solver

## Summary

For a given two- or unidimensional Partial Differential Equation (PDE) specified
in generic form whose derivative does not exceed degree of two and followed by
initial and (generic) boundary conditions for a (generic) geometric
configuration, `pde.solver()` function declared within the source file
`pde.kernel.R` applies Galerkin method with generic (user-defined) weighting
functions and returns numeric solution for reduced system of Ordinary
Differential Equations (ODE) using Runge-Kutta method

## Important Notion

This tool has been primarily designed to be a **Proof-of-Concept** illustrating
that R code can be organized and arranged in a *purely functional manner*. This
styling paradigm is applied in further projects. However, the tool itself is
universal enough (although, maybe not well-optimized) for real-use applications

## Syntax

```
pde.solver(k, a, f, s, W, v, xmin, xmax, ymin, ymax, xpar, ypar,
           pmin, pmax, sector, dirich, reltol, abstol, t, init)
```

where

- `k`, `a`, `f`, `s`, `W` are factor functions for generic PDE, with the number
  of arguments corresponding to the dimension of the prolem;

- `v` is a list of base functions used for Galerkin decomposition with the
  number of arguments corresponding to the dimension of the prolem;

- `xmin`, `xmax`, `ymin`, `ymax` are numeric values specifying boundaries
  (NULL values are acceptable);

- `xpar`, `ypar` are time-dependent functions specifying initial conditions on
  boundaries (NULL values are acceptable);

- `pmin`, `pmax` are time variable constraints;

- `sector` is a logical value specifying whether the problem is being defined in
  polar reference frame;

- `dirich` is a logical value specifying whether the boundary value problem is a
  Dirichlet problem;

- `reltol`, `abstol` are numeric values specifying relative and absolute
  tolerances respectively;

- `t` is a vector specifying a grid for the time variable located between `pmin`
  and `pmax`;

- `init` is a list of functions whose dimensions are agreed with the problem's
  dimension specifying initial conditions if the problem is time-dependent
  (`init = list()` should be defined otherwise)

## Output

Vector of proper dimension whose cells represent the value of supposed solution
for appropriate time and space coordinate. Result is being returned in invisible
form

## Demo

For the sake of sample demonstrations, use `pde.launch.run(n = 1L, src)`
function located within `pde.launch.R`. This function wraps the process of
generating solution with progress messaging and time measurement, and further
graphical constructor (if specified) is being invoked to build and display
animated movie visualizing obtained numeric solution

Arguments:

- `n` is supposed type of the equation (`1L` means elliptic, `2L` means
  parabolic, `3L` means hyperbolic);

- `src` is a source R-file containing implementation of internal
  `.pde.sample.deploy()` function intended to deploy appropriate environment.
  See self-explanatory `pde.sample.dirichlet-2d.R` and `pde.sample.robin-2d.R`
  as working examples for various boundary value problems

## License

MIT License
