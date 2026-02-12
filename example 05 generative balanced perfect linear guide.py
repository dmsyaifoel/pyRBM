from prbm import PRBM

mm = 1e-3 # m
cm = 1e-2 # m

import scipy.optimize as so
import numpy as np
import matplotlib.pyplot as mp

t = 1*mm
h = 3*mm
A = t*h
E = 1650e6
I = h*t**3/12

nsamples = 6

def f(x, retl=False):
  z = PRBM(2)
  z.add_body('A', (0, 0))
  z.add_body('B', (0, 2*cm))
  z.add_body('F', (0, 0))

  z.add_body('C', (-x[0], x[1]))
  z.add_body('D', (x[0], x[1]))
  z.add_body('E', (0, x[2]))

  z.add_body('X', (-x[0], -x[1]))
  z.add_body('Y', (x[0], -x[1]))
  z.add_body('Z', (0, -x[2]))

  z.add_flexure('A', (-2*cm, 0), 'B', (-2*cm, 0))
  z.add_flexure('A', (2*cm, 0), 'B', (2*cm, 0))

  z.add_flexure('F', (-cm, 0), 'B', (-cm, 0))
  z.add_flexure('F', (cm, 0), 'B', (cm, 0))

  z.add_flexure('E', (0, 0), 'C', (0, 0))
  z.add_flexure('E', (0, 0), 'D', (0, 0))
  z.add_flexure('C', (0, 0), 'F', (0, 0))
  z.add_flexure('D', (0, 0), 'F', (0, 0))

  z.add_flexure('Z', (0, 0), 'X', (0, 0))
  z.add_flexure('Z', (0, 0), 'Y', (0, 0))
  z.add_flexure('X', (0, 0), 'F', (0, 0))
  z.add_flexure('Y', (0, 0), 'F', (0, 0))

  if retl: z.show()
  z.move('E', (0, x[3]))
  z.move('Z', (0, -x[3]))
  if retl:
    z.solve_pose('BFCDXY', A, E, I)
    z.show()

  l = []
  for i in range(nsamples):
    z.move('F', (i*mm, 0))
    z.solve_pose('BCDXY', A, E, I)
    l.append(z.energy(A, E, I))

  s = 0
  for i in range(len(l)-1):
    s += (l[i] - l[i+1])**2

  if retl: return l
  print(s)
  return s

sol = so.minimize(f, np.array([cm, cm, -2*cm, -15*mm]), bounds=((5*mm, 2*cm), (-2*cm, -5*mm), (-2*cm, -5*mm), (-2*cm, -5*mm)))
print(sol.x)

l = f(sol.x, True)
print(l)
mp.plot(l)
mp.ylabel('V (J)')
mp.xlabel('x (mm)')
mp.show()
