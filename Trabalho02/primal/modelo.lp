\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj: 6 x1 + 5 x2 + 8 x3 + 13 x4 + 12 x5 + x6 + 7 x7 + 9 x8 + 5 x9 + 10 x10
      + 6 x11 + 4 x12
Subject To
 c1: x1 + x2 + x3 <= 10
 c2: x4 + x5 + x6 <= 20
 c3: x7 + x8 + x9 <= 12
 c4: x10 + x11 + x12 <= 13
 c5: x1 + x4 + x7 + x10 >= 8
 c6: x2 + x5 + x8 + x11 >= 32
 c7: x3 + x6 + x9 + x12 >= 15
End
