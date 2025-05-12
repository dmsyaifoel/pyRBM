from matrix import Matrix
from functions import norm
from auto import syms

a, b, c = syms('a, b, c')

M = Matrix([a, b, c, 1])
print(norm(M))