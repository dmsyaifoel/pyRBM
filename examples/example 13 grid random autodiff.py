import functions
functions.use_custom_instead()
from constants import cm, mm
from grid import Grid
from random import random
from sym import symlist, parsedic
from functions import zeros
from matrix import matrix

N = 7
m, n = N, N
dx, dy = 10*cm/N, 10*cm/N

f = Grid(m, n, dx, dy)
f.set_attribute('fix', 'x<1mm')
f.set_attribute('force', '.06<x,.06<y', (1, 3))

nflex = f.nflex()
ndof = f.ndof()

# f.show()

h = 2*mm
t = [mm*random() for i in range(nflex)]
A = [h*t_ for t_ in t]
E = 1650e6
I = [h*t_**3/12 for t_ in t]

Elist = ndof*[E]

f.show(True, t, 0, mm)
q = symlist('q', ndof)
f.move(q)
en = f.energy(A, Elist, I)

def fg(q):
  return en.valgrad(parsedic({'q':list(q)}), ndof, {})

def f(q):
  return en.val(parsedic({'q':list(q)}))

counter = 0
def cb(x):
  global counter
  print(counter)
  counter += 1

# print(g(zeros(ndof)))

from opt import line_search
sol = line_search(f, fg, zeros(ndof))
print(sol)
#
# f.show(True, t, 0, mm)