from prbm import PRBM
from constants import mm, cm
from math import sin, cos, pi

import scipy.optimize as so
import numpy as np
import matplotlib.pyplot as mp

t = 1*mm
h = 2*mm
E = 1650e6
A = t*h
I = h*t**3/12

def f(x, return_l=False):
  z = PRBM(2)
  z.add_body('A', (0, 0))
  z.add_body('B', (0, 2*cm))

  z.add_body('C', (-x[0], x[1]))
  z.add_body('D', (x[0], x[1]))
  z.add_body('E', (0, x[2]))

  z.add_flexure('A', (-2*cm, 0), 'B', (-2*cm, 0))
  z.add_flexure('A', (2*cm, 0), 'B', (2*cm, 0))
  z.add_flexure('E', (0, 0), 'C', (0, 0))
  z.add_flexure('E', (0, 0), 'D', (0, 0))
  z.add_flexure('C', (0, 0), 'B', (0, 0))
  z.add_flexure('D', (0, 0), 'B', (0, 0))

  if return_l:
    z.show()

  z.move('E', (0, x[3]))

  if return_l:
    z.solve_pose('CD', A, E, I)
    z.show()

  l = []

  if return_l:
    ran = range(-10, 11)
  else:
    ran = range(0, 11, 5)
  for i in ran:
    z.move('B', (.85*2*cm*sin(i*pi/180), .85*2*cm*cos(i*pi/180)+.15*2*cm))
    z.solve_pose('CD', A, E, I)
    l.append(z.energy(A, E, I))


  s = 0
  for i in range(len(l)-1):
    s += (l[i+1] - l[i])**2
  if return_l: return s, l
  print(s)
  return s

sol = so.minimize(f, np.array([cm, 3*cm, 4*cm, 35*mm]), bounds=((5*mm, 15*mm), (25*mm, 35*mm), (5*mm, 4*cm), (5*mm, 4*cm)))
print(sol.x)

s, l = f(sol.x, True)
print(s, np.array(l))
mp.plot(l)
mp.show()