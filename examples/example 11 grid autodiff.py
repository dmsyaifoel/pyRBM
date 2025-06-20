import functions
functions.use_custom_instead()
from constants import cm, mm
from grid import Grid
from sym import symlist, parsedic
from functions import zeros

N = 7
m, n = N, N
dx, dy = 10*cm/N, 10*cm/N
h = 2*mm
t = .1*mm
A = h*t
E = 1650e6
I = h*t**3/12

f = Grid(m, n, dx, dy)
f.set_attribute('fix', 'x<1mm')
f.set_attribute('force', '.06<x,.06<y', (1, 3))
nflex = f.nflex()
ndof = f.ndof()

Alist = nflex*[A]
Elist = nflex*[E]
Ilist = nflex*[I]

q = symlist('q', ndof)
f.move(q)
en = f.energy(Alist, Elist, Ilist)

def fg(q):
  return en.valgrad(parsedic({'q':list(q)}), ndof, {})

def f(q):
  return en.val(parsedic({'q':list(q)}))

from opt import line_search
sol = line_search(f, fg, zeros(ndof))
print(sol)
