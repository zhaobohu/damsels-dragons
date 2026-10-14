# Stochastic character mapping without rejection sampling: endpoint-conditioned branch paths
# by uniformization (Hobolth & Stone 2009). Returns a phytools-style simmap ($maps, $mapped.edge).

suppressPackageStartupMessages({ library(ape); library(expm) })

simmap_uniform <- function(tree, x, Q, pi, nsim = 1) {
  stopifnot(is.rooted(tree))
  # floor zero-length branches at a tiny fraction of tree depth (P(t) = I would block sampling)
  floor_len <- 1e-6 * max(node.depth.edgelength(tree))
  tree$edge.length[tree$edge.length < floor_len] <- floor_len

  # compute in postorder, but return maps in the caller's edge order
  orig <- tree
  tree <- reorder(tree, "postorder")
  ekey <- function(e) paste(e[, 1], e[, 2])
  back <- match(ekey(orig$edge), ekey(tree$edge))   # orig edge i == postorder edge back[i]
  stopifnot(!anyNA(back), !anyDuplicated(back))
  ntip <- Ntip(tree); nn <- ntip + tree$Nnode
  S  <- ncol(Q); st <- colnames(Q)
  stopifnot(identical(colnames(x), st))
  if (!is.null(names(pi))) { stopifnot(setequal(names(pi), st)); pi <- pi[st] }
  stopifnot(length(pi) == S, abs(sum(pi) - 1) < 1e-8)
  # align tip rows to the tree's tip order by name
  stopifnot(setequal(rownames(x), tree$tip.label))
  x  <- x[tree$tip.label, , drop = FALSE]

  # pruning likelihoods; P_e = expm(Q * t_e) cached by unique branch length
  el  <- tree$edge.length
  ul  <- unique(el)
  Pl  <- lapply(ul, function(t) expm(Q * t))
  Pe  <- Pl[match(el, ul)]
  L   <- matrix(0, nn, S)
  L[seq_len(ntip), ] <- x
  # postorder (children before parents), rescaled at every node to avoid underflow
  init <- rep(FALSE, nn)
  for (k in seq_len(nrow(tree$edge))) {
    p <- tree$edge[k, 1]; c <- tree$edge[k, 2]
    if (!init[p]) { L[p, ] <- 1; init[p] <- TRUE }
    L[p, ] <- L[p, ] * as.vector(Pe[[k]] %*% L[c, ])
    m <- max(L[p, ]); if (m > 0) L[p, ] <- L[p, ] / m
  }
  root <- ntip + 1L

  # uniformization constants
  mu <- max(-diag(Q))
  R  <- diag(S) + Q / mu

  Rpow <- function(n) {          # R^n with a small cache
    key <- as.character(n)
    if (!is.null(.rc[[key]])) return(.rc[[key]])
    v <- if (n == 0) diag(S) else R %*% Rpow(n - 1)
    .rc[[key]] <- v; v
  }
  .rc <- new.env(parent = emptyenv())

  # endpoint-conditioned path on one branch: returns (states, times)
  path <- function(a, b, t, Pab) {
    lam <- mu * t
    # number of virtual jumps N | a, b, t
    if (a == b && Pab <= 0) return(list(s = a, len = t))
    # cap on virtual jumps, so the Poisson tail beyond it is negligible for this branch
    ncap <- max(50L, qpois(1 - 1e-12, lam) + 20L)
    n <- 0L; cum <- 0; u <- runif(1)
    repeat {
      pn <- dpois(n, lam) * Rpow(n)[a, b] / Pab
      cum <- cum + pn
      if (cum >= u || n >= ncap) break
      n <- n + 1L
    }
    N <- n
    if (N == 0L) return(list(s = a, len = t))
    # states at the virtual jumps
    s <- integer(N + 1L); s[1] <- a; s[N + 1L] <- b
    for (k in seq_len(N - 1L)) {
      w <- R[s[k], ] * Rpow(N - k)[, b]
      s[k + 1L] <- sample.int(S, 1, prob = w)
    }
    tm <- sort(runif(N, 0, t))
    # collapse virtual self-jumps into runs of one state
    bounds <- c(0, tm, t)
    K <- which(c(TRUE, diff(s) != 0))
    starts <- bounds[K]
    ends   <- c(bounds[K[-1]], t)
    list(s = s[K], len = ends - starts)
  }

  out <- vector("list", nsim)
  for (sim in seq_len(nsim)) {
    # sample node states top-down (root from pi * L)
    ns <- integer(nn)
    pr <- pi * L[root, ]; ns[root] <- sample.int(S, 1, prob = pr / sum(pr))
    for (k in rev(seq_len(nrow(tree$edge)))) {     # preorder = reverse postorder
      p <- tree$edge[k, 1]; c <- tree$edge[k, 2]
      w <- Pe[[k]][ns[p], ] * L[c, ]
      ns[c] <- sample.int(S, 1, prob = w / sum(w))
    }

    # sample the path on every edge
    maps <- vector("list", nrow(tree$edge))
    for (k in seq_len(nrow(tree$edge))) {
      a <- ns[tree$edge[k, 1]]; b <- ns[tree$edge[k, 2]]
      pt <- path(a, b, el[k], Pe[[k]][a, b])
      v <- pt$len; names(v) <- st[pt$s]
      maps[[k]] <- v
    }

    # back to the caller's edge order, plus mapped.edge and node states
    maps <- maps[back]
    m <- orig; m$maps <- maps
    me <- matrix(0, nrow(orig$edge), S, dimnames = list(NULL, st))
    for (k in seq_along(maps)) for (nm in names(maps[[k]]))
      me[k, nm] <- me[k, nm] + maps[[k]][[nm]]
    m$mapped.edge <- me
    m$node.states <- st[ns]                # node numbering is unaffected by reorder
    class(m) <- c("simmap", "phylo")
    out[[sim]] <- m
  }
  if (nsim == 1) out[[1]] else { class(out) <- c("multiSimmap", "multiPhylo"); out }
}
