from constants import cm, mm
from grid import Grid

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

# f.show(True)

Alist = f.ndof()*[A]
Elist = f.ndof()*[E]
Ilist = f.ndof()*[I]

def func(x):
  return f.x_to_energy(x, Alist, Elist, Ilist)

counter = 0
def cb(x):
  global counter
  print(counter)
  counter += 1

from scipy.optimize import minimize
sol = minimize(func, f.ndof()*[0], callback=cb)
print(sol)
func(sol.x)
f.show(True)