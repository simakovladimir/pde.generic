# Copyright 2015 Vladimir Simakov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

library(pracma)
library(Deriv)
library(deSolve)

.pde.deploy.1d.dirichlet <- function(env) {
  
  delayedAssign(
    "scalar.L",
    as.function(
      alist(
        g = , h = ,
        integral(
          Vectorize(
            as.function(
              alist(
                x = ,
                g(x) * h(x)))),
          xmin, xmax, , , , reltol, abstol))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "scalar.H",
    as.function(
      alist(
        g = , h = ,
        as.function(
          alist(
            grad.g = Deriv(g), grad.h = Deriv(h),
            integral(
              Vectorize(
                as.function(
                  alist(
                    x = ,
                    k(x) * grad.g(x) * grad.h(x) + a(x) * g(x) * h(x)))),
              xmin, xmax, , , , reltol, abstol)))())),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "bv.apply",
    as.function(
      alist(
        g = ,
        -scalar.H(g, W))),
    eval.env = env, assign.env = env)
}

.pde.deploy.2d.dirichlet <- function(env) {
  
  delayedAssign(
    "scalar.L",
    as.function(
      alist(
        g = , h = ,
        `$`(
          integral2(
            Vectorize(
              as.function(
                alist(
                  x = , y = ,
                  g(x, y) * h(x, y)))),
            xmin, xmax, ymin, ymax, sector, reltol, abstol,
            singular = TRUE),
          Q))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "scalar.H",
    as.function(
      alist(
        g = , h = ,
        as.function(
          alist(
            grad.g = Deriv(g), grad.h = Deriv(h),
            `$`(
              integral2(
                Vectorize(
                  as.function(
                    alist(
                      x = , y = ,
                      `+`(
                        `*`(
                          sum(
                            `*`(
                              unname(grad.g(x, y)),
                              unname(grad.h(x, y)))),
                          k(x, y)),
                        a(x, y) * g(x, y) * h(x, y))))),
                xmin, xmax, ymin, ymax, sector, reltol, abstol,
                singular = TRUE),
              Q)))())),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "bv.apply",
    as.function(
      alist(
        g = ,
        -scalar.H(g, W))),
    eval.env = env, assign.env = env)
}

.pde.deploy.1d.neumann <- function(env) {
  
  delayedAssign(
    "scalar.L",
    as.function(
      alist(
        g = , h = ,
        integral(
          Vectorize(
            as.function(
              alist(
                x = ,
                g(x) * h(x)))),
          xmin, xmax, , , , reltol, abstol))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "scalar.H",
    as.function(
      alist(
        g = , h = ,
        as.function(
          alist(
            grad.g = Deriv(g), grad.h = Deriv(h),
            `+`(
              integral(
                Vectorize(
                  as.function(
                    alist(
                      x = ,
                      k(x) * grad.g(x) * grad.h(x) + a(x) * g(x) * h(x)))),
                xmin, xmax, , , , reltol, abstol),
              `-`(
                k(xmax) * s(xmax) * g(xmax) * h(xmax),
                k(xmin) * s(xmin) * g(xmin) * h(xmin)))))())),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "bv.apply",
    as.function(
      alist(
        g = ,
        k(xmax) * W(xmax) * g(xmax) - k(xmin) * W(xmin) * g(xmin))),
    eval.env = env, assign.env = env)
}

.pde.deploy.2d.neumann <- function(env) {
  
  delayedAssign(
    "scalar.L",
    as.function(
      alist(
        g = , h = ,
        `$`(
          integral2(
            Vectorize(
              as.function(
                alist(
                  x = , y = ,
                  g(x, y) * h(x, y)))),
            xmin, xmax, ymin, ymax, sector, reltol, abstol,
            singular = TRUE),
          Q))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "scalar.H",
    as.function(
      alist(
        g = , h = ,
        as.function(
          alist(
            grad.g = Deriv(g), grad.h = Deriv(h),
            `+`(
              `$`(
                integral2(
                  Vectorize(
                    as.function(
                      alist(
                        x = , y = ,
                        `+`(
                          `*`(
                            sum(
                              `*`(
                                unname(grad.g(x, y)),
                                unname(grad.h(x, y)))),
                            k(x, y)),
                          a(x, y) * g(x, y) * h(x, y))))),
                  xmin, xmax, ymin, ymax, sector, reltol, abstol,
                  singular = TRUE),
                Q),
              integral(
                Vectorize(
                  as.function(
                    alist(
                      t = ,
                      `*`(
                        k(xpar(t), ypar(t)) * s(xpar(t), ypar(t)),
                        g(xpar(t), ypar(t)) * h(xpar(t), ypar(t)))))),
                pmin, pmax, , , , reltol, abstol))))())),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "bv.apply",
    as.function(
      alist(
        g = ,
        integral(
          Vectorize(
            as.function(
              alist(
                t = ,
                `*`(
                  k(xpar(t), ypar(t)) * W(xpar(t), ypar(t)),
                  g(xpar(t), ypar(t)))))),
          pmin, pmax, , , , reltol, abstol))),
    eval.env = env, assign.env = env)
}

.pde.deploy.generic <- function(env) {
  
  switch(
    1L + (!is.null(env$ymin) && !is.null(env$ymax)) + 2L * !env$dirich,
    .pde.deploy.1d.dirichlet(env),
    .pde.deploy.2d.dirichlet(env),
    .pde.deploy.1d.neumann(env),
    .pde.deploy.2d.neumann(env))
  
  delayedAssign(
    "n",
    length(v),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H",
    outer(v, v, Vectorize(scalar.H)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "H.inv",
    chol2inv(chol(H)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "L",
    outer(v, v, Vectorize(scalar.L)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "L.inv",
    chol2inv(chol(L)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "A",
    L.inv %*% H,
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "A.hyp",
    rbind(
      cbind(matrix(0, n, n), -diag(n)),
      cbind(A, matrix(0, n, n))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "R",
    sapply(v, scalar.L, f) + sapply(v, bv.apply),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "B",
    drop(L.inv %*% R),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "B.hyp",
    c(rep(0, n), B),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "w.0",
    as.numeric(
      do.call(
        c,
        lapply(
          init,
          as.function(
            alist(
              fun = ,
              drop(L.inv %*% sapply(v, scalar.L, fun))))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "fun.par",
    as.function(
      alist(
        t = , y = , parms = ,
        list(B - drop(A %*% y)))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "fun.hyp",
    as.function(
      alist(
        t = , y = , parms = ,
        list(B.hyp - drop(A.hyp %*% y)))),
    eval.env = env, assign.env = env)
}

pde.solver <- function(k, a, f, s, W, v, xmin, xmax, ymin, ymax, xpar, ypar,
                       pmin, pmax, sector, dirich, reltol, abstol, t, init) {
  
  .pde.deploy.generic(environment())
  
  w <- switch(
    1L + length(w.0) %/% n,
    base::t(forceAndCall(2L, replicate, length(t), drop(H.inv %*% R))),
    lsoda(w.0, t, fun.par)[, 1L + seq_len(n)],
    lsoda(w.0, t, fun.hyp)[, 1L + seq_len(n)])
  
  return(w)
}
