\ENCODING=ISO-8859-1
\Problem name: IloCplex

Maximize
 obj: - 6 x1 - 8 x2 - 10 x3 + 4 x4 + 7 x5 + 6 x6 + 7 x7
Subject To
 c1:  - x1 + x4 <= 1
 c2:  - x1 + x5 <= 2
 c3:  - x1 + x6 <= 3
 c4:  - x1 + x7 <= 4
 c5:  - x2 + x4 <= 4
 c6:  - x2 + x5 <= 3
 c7:  - x2 + x6 <= 2
 c8:  - x2 + x7 <= 4
 c9:  - x3 + x4 <= 0
 c10: - x3 + x5 <= 2
 c11: - x3 + x6 <= 2
 c12: - x3 + x7 <= 1
Bounds
      x1 >= 0
      x2 >= 0
      x3 >= 0
      x4 >= 0
      x5 >= 0
      x6 >= 0
      x7 >= 0
Generals
 x1  x2  x3  x4  x5  x6  x7 
End
