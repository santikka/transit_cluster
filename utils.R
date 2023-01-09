children <- function(x, g, v = igraph::V(g)) {
    ch_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "out"))
    v[ch_ind]$name
}

parents <- function(x, g, v = igraph::V(g)) {
    pa_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "in"))
    v[pa_ind]$name
}

descendants <- function(x, g, v = igraph::V(g)) {
    de_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "out"))
    v[de_ind]$name
}

ancestors <- function(x, g, v = igraph::V(g)) {
    an_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "in"))
    v[an_ind]$name
}

neighbors_ <- function(x, g, v = igraph::V(g)) {
    ne_ind <- unlist(igraph::neighborhood(g, order = 1, nodes = x, mode = "all"))
    v[ne_ind]$name
}

connected <- function(x, g, v = igraph::V(g)) {
    co_ind <- unlist(igraph::neighborhood(g, order = length(v), nodes = x, mode = "all"))
    v[co_ind]$name
}

uu <- function(x) {
    if (length(x)) unique(unlist(x))
    else character(0)
}

edge_subgraph <- function(g, incoming, outgoing) {
    # Setting from and to to NULL to satisfy CRAN if we end up making a package
    # R thinks these are global bindings, but they are igraph-operators for edges
    .to <- .from <- NULL
    e <- igraph::E(g)
    e_inc <- e[.to(incoming)]
    e_out <- e[.from(outgoing)]
    igraph::subgraph.edges(g, e[setdiff(e, union(e_inc, e_out))], delete.vertices = FALSE)
}

# Convert an igraph graph using causaleffect syntax into a dag
# with explicit latent variables
to_dag <- function(g) {
    out <- g
    unobs_edges <- which(igraph::edge.attributes(g)$description == "U")
    if (length(unobs_edges)) {
        e <- igraph::get.edges(g, unobs_edges)
        e <- e[e[ ,1] > e[ ,2], , drop = FALSE]
        e_len <- nrow(e)
        new_nodes <- paste0("U[", 1:e_len, "]")
        g <- igraph::set.vertex.attribute(g, name = "description", value = "")
        g <- g + igraph::vertices(new_nodes, description = rep("U", e_len))
        v <- igraph::get.vertex.attribute(g, "name")
        g <- g + igraph::edges(c(rbind(new_nodes, v[e[ ,1]]), rbind(new_nodes, v[e[ ,2]])))
        obs_edges <- setdiff(igraph::E(g), igraph::E(g)[unobs_edges])
        out <- igraph::subgraph.edges(g, igraph::E(g)[obs_edges], delete.vertices = FALSE)
    }
    out
}

# Convert a dag with explicit latent variables into an 
# acyclic directed mixed graph with causaleffect igraph syntax 
to_admg <- function(g) {
    out <- g
    unobs_vars <- which(igraph::vertex.attributes(g)$description == "U")
    obs_vars <- setdiff(1:length(igraph::V(g)), unobs_vars)
    u <- length(unobs_vars)
    if (u) {
        e <- igraph::E(g)
        g_obs <- igraph::induced_subgraph(g, obs_vars)
        for (i in 1:u) {
            unobs_edges <- igraph::get.edges(g, e[.from(unobs_vars[i])])
            if (nrow(unobs_edges) == 2) {
                g_obs <- g_obs + igraph::edges(c(unobs_edges[1:2,2], unobs_edges[2:1,2]), description = "U")
            }
        }
        out <- g_obs
    }
    out
}

random_dag <- function(n, pedge = 0.5, force = TRUE) {
    if (force) {
        while (TRUE) {
            dag <- random_dag(n, pedge, force = FALSE)
            G <- dag$E
            if (!any(rowSums(dag$E) == 0 &
                     colSums(dag$E) == 0)) {
                break
            }
        }
        out_E <- igraph::graph_from_adjacency_matrix(t(dag$E))
        out_E <- igraph::set.vertex.attribute(out_E, "name", value = dag$names)
        return(out_E)
    }
    dag <- list(E = array(0, c(n, n)))
    dag$order <- sample(n, n)
    dag$E <- array(sample(c(0, 1), n^2, prob = c(1 - pedge, pedge), replace = TRUE), c(n, n))
    dag$E[upper.tri(dag$E)] <- 0
    diag(dag$E) <- 0
    if (n < 18) {
        dag$names <- c("X", "Y", "Z", "W", "V", "A", 
                       "B", "C", "D", "E", "F", "G",
                       "H", "I", "J", "K", "L", "M")[1:n]
    } else {
        dag$names <- paste("X_", 1:n, sep = "")
    }
    dag
}

random_causal_diagram <- function(n, pedge = 0.5, pconf = 0.25, force = TRUE) {
    if (force) {
        link <- FALSE
        while (!link) {
            dag <- random_causal_diagram(n, pedge, pconf, force = FALSE)
            G <- dag$E
            # Ensure that graph is connected
            if (any(rowSums(dag$E) == 0 &
                    colSums(dag$E) == 0 &
                    colSums(dag$B) == 0)) {
                next
            }
            diag(G) <- 1
            for (i in 1:n) {
                G <- G %*% G
            }
            # Ensure path X -> ... -> Y
            link <- G[2, 1] != 0
        }
        out_E <- igraph::graph_from_adjacency_matrix(t(dag$E))
        out_E <- igraph::set.vertex.attribute(out_E, "name", value = dag$names)
        out_E <- igraph::set.edge.attribute(out_E, "description", value = "")
        out_B <- igraph::graph_from_adjacency_matrix(t(dag$B))
        out_B <- igraph::set.vertex.attribute(out_B, "name", value = dag$names)
        #return(list(out_E, out_B))
        return(out_E + igraph::edges(t(igraph::get.edgelist(out_B)), description = "U"))
    }
    dag <- list(E = array(0, c(n, n)), B = array(0, c(n, n)))
    dag$order <- sample(n, n)
    dag$E <- array(sample(c(0, 1), n^2, prob = c(1 - pedge, pedge), replace = TRUE), c(n, n))
    dag$E[upper.tri(dag$E)] <- 0
    diag(dag$E) <- 0
    dag$E[dag$order, dag$order] <- dag$E
    dag$B <- array(sample(c(0, 1), n^2, prob = c(1 - pconf, pconf), replace = TRUE), c(n, n))
    dag$B[upper.tri(dag$B)] <- 0
    diag(dag$B) <- 0
    dag$B <- dag$B + t(dag$B)
    if (n < 18) {
        dag$names <- c("X", "Y", "Z", "W", "V", "A", 
                       "B", "C", "D", "E", "F", "G",
                       "H", "I", "J", "K", "L", "M")[1:n]
    } else {
        dag$names <- paste("X_", 1:n, sep = "")
    }
    dag
}

