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

.pde.sample.deploy <- function(env) {
  
  delayedAssign(
    "keyword",
    "dirichlet-2d",
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "args.generic",
    list(
      k = as.function(alist(x = , y = , 1)),
      a = as.function(alist(x = , y = , 0)),
      f = as.function(alist(x = , y = , 1 / sqrt(x ^ 2L + y ^ 2L))),
      s = NULL,
      W = as.function(alist(x = , y = , 0)),
      v = polybase(5L),
      xmin = 0,
      xmax = 2 * pi,
      ymin = 0,
      ymax = 1,
      xpar = cos,
      ypar = sin,
      pmin = 0,
      pmax = 2 * pi,
      sector = TRUE,
      dirich = TRUE,
      reltol = 1e-6,
      abstol = 1e-6,
      t = seq(0, 5, 0.05)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "args.specific",
    list(
      list(init = list()),
      list(init = list(gauss(1.0, rep(0.0, 2L), rep(0.25, 2L)))),
      list(init = list(gauss(1.0, rep(0.0, 2L), rep(0.25, 2L)),
                       gauss(0.5, rep(0.0, 2L), rep(0.25, 2L))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "grid",
    list(
      x = seq(-1, 1, 0.01),
      y = seq(-1, 1, 0.01)),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "prepare.pattern",
    as.function(
      alist(
        fun = ,
        outer(grid$x, grid$y, fun))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "exterior",
    outer(grid$x ^ 2L, grid$y ^ 2L, `+`) >= 1,
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "visualize",
    as.function(
      alist(
        i = , u = ,
        filled.contour(grid$x, grid$y, u[[i]], col = brewer.pal(9L, "Blues"),
                       levels = pretty(range(u, na.rm = TRUE, finite = TRUE),
                                       9L)))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "record.movie",
    as.function(
      alist(
        u = ,
        saveGIF(
          sapply(seq_len(length(u)), visualize, u),
          movie.name = paste("pde.simulation", keyword, n, types[n], "gif",
                             sep = "."),
          interval = 0.12,
          outdir = getwd(),
          verbose = TRUE,
          autobrowse = TRUE))),
    eval.env = env, assign.env = env)
}
