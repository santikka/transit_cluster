source("clustering.R")
library(igraph)

g1 <- graph.formula(X -+ Y, Z_1 -+ X, Z_1 -+ Y, Z_2 -+ X, Z_2 -+ Y)

g2 <- graph.formula(U_1 -+ X_1, U_1 -+ R_1, U_1 -+ R_2, X_2 -+ R_1,
                    X_2 -+ R_2, X_2 -+ Y_2, X_1 -+ R_1, X_1 -+ R_2,
                    R_1 -+ A_3, A_3 -+ R_2, A_2 -+ R_2, A_2 -+ E_2,
                    R_2 -+ E_2, R_1 -+ A_1, A_1 -+ E_1, E_1 -+ A_4,
                    A_4 -+ E_2, E_1 -+ Y_1, E_1 -+ Y_2, E_2 -+ Y_1,
                    E_2 -+ Y_2, E_1 -+ E_2)

g3 <- graph.formula(U_1 -+ X_1, U_1 -+ R_1, U_1 -+ R_2, X_2 -+ R_1,
                    X_2 -+ R_2, X_2 -+ Y_2, X_1 -+ R_1, X_1 -+ R_2,
                    R_1 -+ A_3, A_3 -+ R_2, A_2 -+ R_2, A_2 -+ E_2,
                    R_2 -+ E_2, R_1 -+ A_1, A_1 -+ E_1, E_1 -+ A_4,
                    A_4 -+ E_2, E_1 -+ Y_1, E_1 -+ Y_2, E_2 -+ Y_1,
                    E_2 -+ Y_2, E_1 -+ E_2, R_1 -+ A_5, A_5 -+ A_1)

g4 <- graph.formula(U_1 -+ X_1, U_1 -+ R_1, U_1 -+ R_2, X_2 -+ R_1,
                    X_2 -+ R_2, X_2 -+ Y_2, X_1 -+ R_1, X_1 -+ R_2,
                    R_1 -+ A_3, A_3 -+ R_2, A_2 -+ R_2, A_2 -+ E_2,
                    R_2 -+ E_2, R_1 -+ A_1, A_1 -+ E_1, E_1 -+ A_4,
                    A_4 -+ E_2, E_1 -+ Y_1, E_1 -+ Y_2, E_2 -+ Y_1,
                    E_2 -+ Y_2, E_1 -+ E_2, A_5 -+ A_1, A_1 -+ A_6)

g5 <- graph.formula(X_1 -+ R_1 -+ E_1 -+ Y_1,
                    X_2 -+ R_2 -+ E_2 -+ Y_2,
                    X_3 -+ R_3 -+ E_3 -+ Y_3,
                    X_1 -+ R_2, X_2 -+ R_3, X_3 -+ R_1,
                    E_1 -+ Y_2, E_2 -+ Y_3, E_3 -+ Y_1)

g6 <- graph.formula(X_1 -+ R_1 -+ E_1 -+ Y_1,
                    X_1 -+ R_2 -+ E_2 -+ Y_1,
                    X_2 -+ R_2, X_2 -+ Y_1)

g7 <- graph.formula(X -+ V_1 -+ V_2 -+ Y, V_1 +- U -+ V_2)

g8 <- graph.formula(X -+ R_3 -+ R_2 -+ R_1 -+ E_1 -+ E_2 -+ Y,
                    X -+ R_2, X -+ R_1, U -+ R_3, U -+ E_2, E_1 -+ Y)

g9 <- graph.formula(X -+ R_1 -+ I_1 -+ E_1 -+ Y,
                    X -+ R_2 -+ I_2 -+ E_2 -+ Y,
                    U -+ I_2, U -+ Y)

g10 <- graph.formula(X -+ R_3 -+ R_2 -+ R_1 -+ E_1 -+ Y, U -+ R_3, U -+ Y,
                     X -+ R_2, X -+ R_1, R_2 -+ R_4, R_4 -+ E_2 -+ Y)

g11 <- graph.formula(A -+ B -+ C -+ D)

g12 <- graph.formula(A -+ B -+ C -+ D, A -+ C, A -+ D, B -+ D)

g13 <- graph.formula(A -+ B_1 -+ C_1 -+ D, A -+ B_2 -+ C_2 -+ D)

# Two components, but can be merged post hoc
find_transit_components(g1, prohibit = c("X", "Y"))

# Variations of the identifiability proof graph sketch
find_transit_components(g2, prohibit = c("X_1", "X_2", "Y_1", "Y_2", "U_1"))
find_transit_components(g3, prohibit = c("X_1", "X_2", "Y_1", "Y_2", "U_1"))
find_transit_components(g4, prohibit = c("X_1", "X_2", "Y_1", "Y_2", "U_1"))

# More groups found compared to algorithm v2!
find_transit_components(g5)

# Old counterexample, now found!
find_transit_components(g6, prohibit = c("X_1", "X_2", "Y_1"))

# Example by Jouni, now found!
find_transit_components(g7, prohibit = c("X", "Y"))
find_transit_components(g7)

# Example, where there is a maximal grouping, and its subset
# that is also a grouping
find_transit_components(g8, prohibit = c("X", "Y"))

# Check neighborhood for components only, otherwise R_1 -> I_1 -> E_1 not found!
# NOW FIXED!
find_transit_components(g9, prohibit = c("X", "Y", "U"))

# Two components
find_transit_components(g10, prohibit = c("X", "Y"))

# Maximum amounts
find_transit_components(g11, singletons = TRUE)
find_transit_components(g12, singletons = TRUE)

# Two components found by one pair of representatives
find_transit_components(g13, singletons = TRUE)
