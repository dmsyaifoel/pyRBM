from prbm import PRBM
from constants import mm, cm
from math import sin, cos, pi

p = PRBM(3)

p.add_body('A', (0, 0, 0))
p.add_body('B', (0, 0, 15*mm))
p.add_body('C', (0, 0, 3*cm))

t = mm
A = pi*t**2
E = 1650e6
I = pi*t**4/2

n = 3
r = 14*mm

for i in range(n):
  p.add_flexure('A', (r*cos(i/n*2*pi), r*sin(i/n*2*pi), 0),
                'B', (r*cos((i + 1)/n*2*pi), r*sin((i + 1)/n*2*pi), 0))
  p.add_flexure('C', (r*cos((i - .5)/n*2*pi), r*sin((i - .5)/n*2*pi), 0),
                'B', (r*cos((i + .5)/n*2*pi), r*sin((i + .5)/n*2*pi), 0))

p.add_force('C', (3, 2, -1))

p.solve_pose('BC', A, E, I)

p.print(A, E, I)
p.show()

