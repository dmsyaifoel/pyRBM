from prbm import PRBM
from math import sin, cos, pi

mm = 1e-3 # m
cm = 1e-2 # m

q = PRBM(3)

t = mm
A = pi*t**2
E = 1650e6 # Pa
I = pi*t**4/2

q.add_body('A', (0, 0, 0))
q.add_body('B', (0, 0, 0))

n = 3
ra = 2*cm
rb = cm

for i in range(n):
  q.add_flexure('A', (ra*cos(i/n*2*pi), ra*sin(i/n*2*pi), 0),
                'B', (rb*cos(i/n*2*pi), rb*sin(i/n*2*pi), -3*cm))

q.add_force('B', (1, 2, 3))

q.solve_pose('B', A, E, I)

q.show()
