from matrix import matrix

A = matrix([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]])

print(2*A/3 + A@A - A)