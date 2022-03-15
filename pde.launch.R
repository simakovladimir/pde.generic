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

library(ggplot2)
library(RColorBrewer)
library(animation)

source("./pde.kernel.R")

.pde.launch.deploy <- function(env) {
  
  assign(
    "origin",
    proc.time()["elapsed"],
    envir = env)
  
  delayedAssign(
    "types",
    c(
      "elliptic",
      "parabolic",
      "hyperbolic"),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "msg",
    c(
      "complete\n",
      "-->\n",
      "Galerkin PDE Solver\n",
      "Solving linear system ",
      "Preparing numeric pattern ",
      "Generating solution ",
      "Producing animation frames ",
      "Elapsed time: %d:%02d:%02d\n",
      "Grabbed size: %d bytes\n"),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "msg.box",
    as.function(
      alist(
        num = , complete = , ... = ,
        `{`(
          cat(
            switch(
              1L + is.null(complete),
              paste(
                paste(
                  rep(".", 45L - nchar(sprintf(msg[num], ...))),
                  collapse = ""),
                msg[complete]),
              sprintf(msg[num], ...)),
            sep = ""),
          flush.console()))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "gauss",
    as.function(
      alist(
        a = , mu = , sigma = ,
        forceAndCall(
          1L,
          as.function,
          c(
            eval(
              parse(
                "", NULL,
                paste0(
                  "alist(",
                  paste(
                    paste("x", seq_along(mu), sep = "."), "= ",
                    collapse = ", "),
                  ")"))),
            bquote(
              `*`(
                exp(
                  `*`(
                    sum(
                      `^`(
                        `/`(
                          `-`(
                            unlist(
                              lapply(
                                as.list(match.call())[-1],
                                eval,
                                parent.frame()),
                              FALSE, FALSE),
                            mu),
                          sigma),
                        2L)),
                    -0.5)),
                a)))))),
    eval.env = env, assign.env = env)
  
  delayedAssign(
    "polybase",
    as.function(
      alist(
        deg = ,
        as.function(
          alist(
            template = paste("as.function(alist(x = , y = ,",
                             "(1 - x ^ 2L - y ^ 2L) * x ^ %dL * y ^ %dL))"),
            do.call(
              c,
              lapply(
                seq(0L, deg),
                as.function(
                  alist(
                    i = ,
                    lapply(
                      seq(0L, i),
                      as.function(
                        alist(
                          j = ,
                          forceAndCall(
                            1L,
                            eval,
                            parse(
                              "", NULL,
                              sprintf(template, j, i - j))))))))))))())),
    eval.env = env, assign.env = env)
}

pde.launch.run <- function(n = 1L, src = "./pde.sample.dirichlet-2d.R") {
  
  source(src)
  .pde.launch.deploy(environment())
  .pde.sample.deploy(environment())
  
  msg.box(3L, NULL)
  
  msg.box(4L, NULL)
  w <- do.call(pde.solver, c(args.generic, args.specific[[n]]))
  msg.box(4L, 1L)
  
  msg.box(5L, NULL)
  ptrn <- lapply(
    lapply(args.generic$v, prepare.pattern),
    replace, exterior, 0)
  msg.box(5L, 1L)
  
  msg.box(6L, NULL)
  u <- lapply(
    lapply(
      seq_len(nrow(w)),
      as.function(alist(i = , Reduce(`+`, Map(`*`, w[i, ], ptrn))))),
    `+`, prepare.pattern(Vectorize(args.generic$W)))
  msg.box(6L, 1L)
  
  msg.box(7L, NULL)
  msg.box(2L, NULL)
  record.movie(u)
  msg.box(7L, NULL)
  msg.box(7L, 1L)
  
  forceAndCall(
    1L,
    as.function(
      alist(
        d = ,
        msg.box(8L, NULL, d %/% 3600L, (d %/% 60L) %% 60L, d %% 60L))),
    as.integer(round(proc.time()["elapsed"] - origin)))
  
  msg.box(9L, NULL, object.size(u))
  
  return(invisible(u))
}
