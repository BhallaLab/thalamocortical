graph [
  directed 1
  doc "Celltype-based connectivity data. count of node *n* is the number of cells of type *n* that are present in the model. weight of edge (a, b) is the number of cells of type *a* that connect to each cell of type *b*."
  node [
    id 0
    label "SupBasket"
    count 90
    index 12
  ]
  node [
    id 1
    label "DeepBasket"
    count 100
    index 1
  ]
  node [
    id 2
    label "TuftedIB"
    count 800
    index 8
  ]
  node [
    id 3
    label "TCR"
    count 100
    index 3
  ]
  node [
    id 4
    label "SupPyrRS"
    count 1000
    index 4
  ]
  node [
    id 5
    label "TuftedRS"
    count 200
    index 5
  ]
  node [
    id 6
    label "SupAxoaxonic"
    count 90
    index 6
  ]
  node [
    id 7
    label "SupPyrFRB"
    count 50
    index 7
  ]
  node [
    id 8
    label "DeepAxoaxonic"
    count 100
    index 2
  ]
  node [
    id 9
    label "DeepLTS"
    count 100
    index 9
  ]
  node [
    id 10
    label "SpinyStellate"
    count 240
    index 10
  ]
  node [
    id 11
    label "SupLTS"
    count 90
    index 11
  ]
  node [
    id 12
    label "NontuftedRS"
    count 500
    index 0
  ]
  node [
    id 13
    label "nRT"
    count 100
    index 13
  ]
  edge [
    source 0
    target 12
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 1
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 8
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 4
    ps_comps [1, 2, 3, 4, 5, 6, 7, 8, 9, 38, 39]
    weight 20
  ]
  edge [
    source 0
    target 5
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 0
    target 7
    ps_comps [1, 2, 3, 4, 5, 6, 7, 8, 9, 38, 39]
    weight 20
  ]
  edge [
    source 0
    target 2
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 9
    ps_comps []
    weight 0
  ]
  edge [
    source 0
    target 10
    ps_comps [1, 2, 15, 28, 41]
    weight 20
  ]
  edge [
    source 0
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 0
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 0
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 12
    ps_comps [1, 2, 3, 4, 5, 6, 35, 36]
    weight 20
  ]
  edge [
    source 1
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 1
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 1
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 4
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 5
    ps_comps [1, 2, 3, 4, 5, 6, 35, 36]
    weight 20
  ]
  edge [
    source 1
    target 6
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 7
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 2
    ps_comps [1, 2, 3, 4, 5, 6, 35, 36]
    weight 20
  ]
  edge [
    source 1
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 1
    target 10
    ps_comps [1, 2, 15, 28, 41]
    weight 20
  ]
  edge [
    source 1
    target 11
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 0
    ps_comps []
    weight 0
  ]
  edge [
    source 1
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 2
    target 12
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    weight 20
  ]
  edge [
    source 2
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 2
    target 4
    ps_comps [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
    weight 2
  ]
  edge [
    source 2
    target 5
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 20
  ]
  edge [
    source 2
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 7
    ps_comps [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
    weight 2
  ]
  edge [
    source 2
    target 2
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 50
  ]
  edge [
    source 2
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 2
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 3
    target 12
    ps_comps [40, 41, 42, 43, 44]
    weight 10
  ]
  edge [
    source 3
    target 1
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 20
  ]
  edge [
    source 3
    target 8
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 10
  ]
  edge [
    source 3
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 3
    target 4
    ps_comps [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 10
  ]
  edge [
    source 3
    target 5
    ps_comps [47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 10
  ]
  edge [
    source 3
    target 6
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 10
  ]
  edge [
    source 3
    target 7
    ps_comps [45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 10
  ]
  edge [
    source 3
    target 2
    ps_comps [47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 10
  ]
  edge [
    source 3
    target 9
    ps_comps []
    weight 0
  ]
  edge [
    source 3
    target 10
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 0
  ]
  edge [
    source 3
    target 11
    ps_comps []
    weight 0
  ]
  edge [
    source 3
    target 0
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 10
  ]
  edge [
    source 3
    target 13
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 25
  ]
  edge [
    source 4
    target 12
    ps_comps [38, 39, 40, 41, 42, 43, 44]
    weight 3
  ]
  edge [
    source 4
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 30
  ]
  edge [
    source 4
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 30
  ]
  edge [
    source 4
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 4
    target 4
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 10, 11, 12, 13, 22, 23, 24, 25, 34, 35, 36, 37]
    weight 50
  ]
  edge [
    source 4
    target 5
    ps_comps [39, 40, 41, 42, 43, 44, 45, 46]
    weight 60
  ]
  edge [
    source 4
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 90
  ]
  edge [
    source 4
    target 7
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 10, 11, 12, 13, 22, 23, 24, 25, 34, 35, 36, 37]
    weight 50
  ]
  edge [
    source 4
    target 2
    ps_comps [39, 40, 41, 42, 43, 44, 45, 46]
    weight 60
  ]
  edge [
    source 4
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 30
  ]
  edge [
    source 4
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 3
  ]
  edge [
    source 4
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 90
  ]
  edge [
    source 4
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 90
  ]
  edge [
    source 4
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 5
    target 12
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    weight 20
  ]
  edge [
    source 5
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 5
    target 4
    ps_comps [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
    weight 2
  ]
  edge [
    source 5
    target 5
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 10
  ]
  edge [
    source 5
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 7
    ps_comps [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
    weight 2
  ]
  edge [
    source 5
    target 2
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 20
  ]
  edge [
    source 5
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 5
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 12
    ps_comps []
    weight 5
  ]
  edge [
    source 6
    target 1
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 8
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 4
    ps_comps []
    weight 20
  ]
  edge [
    source 6
    target 5
    ps_comps []
    weight 5
  ]
  edge [
    source 6
    target 6
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 7
    ps_comps []
    weight 20
  ]
  edge [
    source 6
    target 2
    ps_comps []
    weight 5
  ]
  edge [
    source 6
    target 9
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 10
    ps_comps []
    weight 5
  ]
  edge [
    source 6
    target 11
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 0
    ps_comps []
    weight 0
  ]
  edge [
    source 6
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 7
    target 12
    ps_comps [38, 39, 40, 41, 42, 43, 44]
    weight 1
  ]
  edge [
    source 7
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 3
  ]
  edge [
    source 7
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 3
  ]
  edge [
    source 7
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 7
    target 4
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 10, 11, 12, 13, 22, 23, 24, 25, 34, 35, 36, 37]
    weight 5
  ]
  edge [
    source 7
    target 5
    ps_comps [39, 40, 41, 42, 43, 44, 45, 46]
    weight 3
  ]
  edge [
    source 7
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 5
  ]
  edge [
    source 7
    target 7
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33, 10, 11, 12, 13, 22, 23, 24, 25, 34, 35, 36, 37]
    weight 5
  ]
  edge [
    source 7
    target 2
    ps_comps [39, 40, 41, 42, 43, 44, 45, 46]
    weight 3
  ]
  edge [
    source 7
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 3
  ]
  edge [
    source 7
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 1
  ]
  edge [
    source 7
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 5
  ]
  edge [
    source 7
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 5
  ]
  edge [
    source 7
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 12
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 1
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 8
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 4
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 5
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 6
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 7
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 2
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 9
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 10
    ps_comps []
    weight 5
  ]
  edge [
    source 8
    target 11
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 0
    ps_comps []
    weight 0
  ]
  edge [
    source 8
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 9
    target 12
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44]
    weight 20
  ]
  edge [
    source 9
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 9
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 9
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 9
    target 4
    ps_comps [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 10
  ]
  edge [
    source 9
    target 5
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 20
  ]
  edge [
    source 9
    target 6
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 10
  ]
  edge [
    source 9
    target 7
    ps_comps [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 10
  ]
  edge [
    source 9
    target 2
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 20
  ]
  edge [
    source 9
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 9
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 9
    target 11
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 10
  ]
  edge [
    source 9
    target 0
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 10
  ]
  edge [
    source 9
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 10
    target 12
    ps_comps [37, 38, 39, 40, 41]
    weight 20
  ]
  edge [
    source 10
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 10
    target 4
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33]
    weight 20
  ]
  edge [
    source 10
    target 5
    ps_comps [7, 8, 9, 10, 11, 12, 36, 37, 38, 39, 40, 41]
    weight 20
  ]
  edge [
    source 10
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 7
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 31, 32, 33]
    weight 20
  ]
  edge [
    source 10
    target 2
    ps_comps [7, 8, 9, 10, 11, 12, 36, 37, 38, 39, 40, 41]
    weight 20
  ]
  edge [
    source 10
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 30
  ]
  edge [
    source 10
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 20
  ]
  edge [
    source 10
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 11
    target 12
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44]
    weight 20
  ]
  edge [
    source 11
    target 1
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 20
  ]
  edge [
    source 11
    target 8
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 20
  ]
  edge [
    source 11
    target 3
    ps_comps []
    weight 0
  ]
  edge [
    source 11
    target 4
    ps_comps [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 20
  ]
  edge [
    source 11
    target 5
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 20
  ]
  edge [
    source 11
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 11
    target 7
    ps_comps [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68]
    weight 20
  ]
  edge [
    source 11
    target 2
    ps_comps [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    weight 20
  ]
  edge [
    source 11
    target 9
    ps_comps [8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 34, 35, 36, 37, 38, 47, 48, 49, 50, 51]
    weight 20
  ]
  edge [
    source 11
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 11
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 11
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 20
  ]
  edge [
    source 11
    target 13
    ps_comps []
    weight 0
  ]
  edge [
    source 12
    target 12
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
    weight 20
  ]
  edge [
    source 12
    target 1
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 8
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 3
    ps_comps [6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 32, 33, 34, 35, 36, 37, 38, 39, 40, 45, 46, 47, 48, 49, 50, 51, 52, 53, 58, 59, 60, 61, 62, 63, 64, 65, 66, 71, 72, 73, 74, 75, 76, 77, 78, 79, 84, 85, 86, 87, 88, 89, 90, 91, 92, 97, 98, 99, 100, 101, 102, 103, 104, 105, 110, 111, 112, 113, 114, 115, 116, 117, 118, 123, 124, 125, 126, 127, 128, 129, 130, 131]
    weight 20
  ]
  edge [
    source 12
    target 4
    ps_comps [41, 42, 43, 44]
    weight 10
  ]
  edge [
    source 12
    target 5
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 10
  ]
  edge [
    source 12
    target 6
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 7
    ps_comps [41, 42, 43, 44]
    weight 10
  ]
  edge [
    source 12
    target 2
    ps_comps [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]
    weight 10
  ]
  edge [
    source 12
    target 9
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 10
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 11
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 0
    ps_comps [5, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23, 31, 32, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49]
    weight 10
  ]
  edge [
    source 12
    target 13
    ps_comps [2, 3, 4, 15, 16, 17, 28, 29, 30, 41, 42, 43]
    weight 20
  ]
  edge [
    source 13
    target 12
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 1
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 8
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 3
    ps_comps [1, 2, 15, 28, 41, 54, 67, 80, 93, 106, 119]
    weight 15
  ]
  edge [
    source 13
    target 4
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 5
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 6
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 7
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 2
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 9
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 10
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 11
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 0
    ps_comps []
    weight 0
  ]
  edge [
    source 13
    target 13
    ps_comps [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]
    weight 10
  ]
]