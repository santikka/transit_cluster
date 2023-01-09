library(causaleffect)
library(igraph)
# grouping
source("clustering.R")

# Sangiovese graph, renamed for shorter equations
# original_name  new_name
#      Treatment        X
#        SproutN      Z_1
#         BunchN      Z_4
#         GrapeW        Y
#          WoodW      Z_7
#         SPAD06      Z_2
#         NDVI06      Z_3
#         SPAD08      Z_5
#         NDVI08      Z_6
#           Acid     Z_12
#         Potass      Z_9
#           Brix     Z_10
#             pH     Z_13
#         Anthoc      Z_8
#         Polyph     Z_11

g <- graph.formula(Z_11 -+ Z_13, Z_11 -+ Z_12, Z_8 -+ Z_11, Z_8 -+ Z_13, 
  Z_8 -+ Z_10, Z_8 -+ Z_9, Z_8 -+ Z_12, Z_8 -+ Y, Z_13 -+ Y, Z_10 -+ Z_11, 
  Z_10 -+ Z_13, Z_10 -+ Z_12, Z_10 -+ Y, Z_9 -+ Z_13, Z_12 -+ Z_13, Z_12 -+ Y, 
  Z_6 -+ Z_11, Z_6 -+ Z_8, Z_6 -+ Z_12, Z_6 -+ Z_7, Z_6 -+ Y, Z_5 -+ Z_6, 
  Z_5 -+ Z_7, Z_3 -+ Z_11, Z_3 -+ Z_12, Z_3 -+ Z_6, Z_3 -+ Z_5, Z_3 -+ Y, 
  Z_2 -+ Z_13, Z_2 -+ Z_9, Z_2 -+ Z_12, Z_2 -+ Z_5, Z_2 -+ Z_3, Z_2 -+ Z_7, 
  Z_7 -+ Z_8, Z_7 -+ Z_13, Z_7 -+ Y, Z_4 -+ Z_11, Z_4 -+ Z_8, Z_4 -+ Z_9, 
  Z_4 -+ Z_12, Z_4 -+ Z_7, Z_4 -+ Y, Z_1 -+ Z_13, Z_1 -+ Z_12, Z_1 -+ Z_6, 
  Z_1 -+ Z_3, Z_1 -+ Z_2, Z_1 -+ Z_7, Z_1 -+ Y, Z_1 -+ Z_4, X -+ Z_10, 
  X -+ Z_2, X -+ Z_4, X -+ Z_1, X -+ Y, Y -+ X, simplify = FALSE)
g <- set.edge.attribute(g, "description", 56:57, "U")

causal.effect(G = g, y = "Y", x = "X")
"\\sum_{Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8,Z_10,Z_9,Z_11,Z_12,Z_13}
\\left(\\sum_{X}P(Y|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8,Z_10,Z_9,Z_11,Z_12,Z_13)
P(X)\\right)P(Z_13|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8,Z_10,Z_9,Z_11,Z_12)
P(Z_12|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8,Z_10,Z_11)
P(Z_11|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8,Z_10)
P(Z_9|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8)
P(Z_10|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7,Z_8)P(Z_8|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6,Z_7)
P(Z_7|X,Z_1,Z_2,Z_4,Z_3,Z_5,Z_6)P(Z_6|X,Z_1,Z_2,Z_3,Z_5)P(Z_5|X,Z_1,Z_2,Z_3)
P(Z_3|X,Z_1,Z_2)P(Z_4|X,Z_1)P(Z_2|X,Z_1)P(Z_1|X)"

# This simplifies to
"\\sum_{Z_1,Z_4,Z_3,Z_6,Z_7,Z_8,Z_10,Z_12,Z_13}
P(Z_1,Z_4,Z_3,Z_6,Z_7,Z_8,Z_10,Z_12,Z_13 | X)
\\left(\\sum_{X}P(Y|X,Z_1,Z_4,Z_3,Z_6,Z_7,Z_8,Z_10,Z_12,Z_13)P(X)\\right)"
# ie standard frontdoor with Z = (Z_1,Z_4,Z_3,Z_6,Z_7,Z_8,Z_10,Z_12,Z_13)...

# explicit U :
g2 <- graph.formula(Z_11 -+ Z_13, Z_11 -+ Z_12, Z_8 -+ Z_11, Z_8 -+ Z_13, 
  Z_8 -+ Z_10, Z_8 -+ Z_9, Z_8 -+ Z_12, Z_8 -+ Y, Z_13 -+ Y, Z_10 -+ Z_11, 
  Z_10 -+ Z_13, Z_10 -+ Z_12, Z_10 -+ Y, Z_9 -+ Z_13, Z_12 -+ Z_13, Z_12 -+ Y, 
  Z_6 -+ Z_11, Z_6 -+ Z_8, Z_6 -+ Z_12, Z_6 -+ Z_7, Z_6 -+ Y, Z_5 -+ Z_6, 
  Z_5 -+ Z_7, Z_3 -+ Z_11, Z_3 -+ Z_12, Z_3 -+ Z_6, Z_3 -+ Z_5, Z_3 -+ Y, 
  Z_2 -+ Z_13, Z_2 -+ Z_9, Z_2 -+ Z_12, Z_2 -+ Z_5, Z_2 -+ Z_3, Z_2 -+ Z_7, 
  Z_7 -+ Z_8, Z_7 -+ Z_13, Z_7 -+ Y, Z_4 -+ Z_11, Z_4 -+ Z_8, Z_4 -+ Z_9, 
  Z_4 -+ Z_12, Z_4 -+ Z_7, Z_4 -+ Y, Z_1 -+ Z_13, Z_1 -+ Z_12, Z_1 -+ Z_6, 
  Z_1 -+ Z_3, Z_1 -+ Z_2, Z_1 -+ Z_7, Z_1 -+ Y, Z_1 -+ Z_4, X -+ Z_10, 
  X -+ Z_2, X -+ Z_4, X -+ Z_1, U -+ Y, U -+ X, simplify = FALSE)

find_transit_components(g2, c("X", "Y"))
# [[1]]
# [[1]]$vertices
# [1] "Z_11" "Z_13" "Z_12" "Z_8"  "Z_10" "Z_9"  "Z_6"  "Z_7"  "Z_5"  "Z_3"  "Z_2"  "Z_4"  "Z_1" 
# 
# [[1]]$receivers
# [1] "Z_10" "Z_2"  "Z_4"  "Z_1" 
# 
# [[1]]$emitters
# [1] "Z_13" "Z_12" "Z_8"  "Z_10" "Z_6"  "Z_7"  "Z_3"  "Z_4"  "Z_1" 

# Graph obtained from clustering:
g3 <- graph.formula(X -+ Z, Z-+ Y, X -+ Y, Y -+ X, simplify = FALSE)
g3 <- set.edge.attribute(g3, "description", 3:4, "U")

## Timing:
microbenchmark::microbenchmark(
  original = causal.effect(G = g, y = "Y", x = "X"),
  simplify = causal.effect(G = g, y = "Y", x = "X", simp = TRUE), 
  both = causal.effect(G = g, y = "Y", x = "X", simp = TRUE, prune = TRUE), 
  clustering = {
    find_transit_components(g2, c("X", "Y")) 
    causal.effect(G = g3, y = "Y", x = "X")
    }, 
  times = 10, unit = "seconds")
# Unit: seconds
# expr         min          lq        mean      median          uq         max neval
#   original   0.1040187   0.1053718   0.1325634   0.1075394   0.1120852   0.3421086    10
#   simplify 119.8897374 126.9114463 134.4376348 132.3178575 143.3740094 149.5324991    10
#       both 118.8405364 124.2968720 133.5803521 134.3901111 142.5337139 146.7986841    10
# clustering   0.2488937   0.2646434   0.4298052   0.4814837   0.5202207   0.6245915    10

causal.effect(G = g3, y = "Y", x = "X")
