from prbm import PRBM
from constants import mm, cm

t = 1*mm
h = 2*mm
E = 1650e6
A = t*h
I = h*t**3/12

p = PRBM(2)

p.add_body('A', (0, 0))
p.add_body('B', (0, 3*cm))
p.add_body('C', (5*cm, 0))
p.add_body('D', (9*cm, 10*cm))
p.add_flexure('A', (4*cm, 0), 'B', (4*cm, 0))
p.add_flexure('A', (-4*cm, 0), 'B', (-4*cm, 0))
p.add_flexure('A', (4*cm, 0), 'B', (4*cm, 0))
p.add_flexure('A', (4*cm, 0), 'C', (0, 0))
p.add_flexure('B', (4*cm, 0), 'C', (0, 3*cm))
p.add_flexure('A', (6*cm, 0), 'D', (-3*cm, 0))
p.add_flexure('A', (12*cm, 0), 'D', (3*cm, 0))
p.add_flexure('C', (0, 10*cm), 'D', (-3*cm, 0))
p.show()

p.add_force('B', (1, 0))
p.solve_pose('BCD', A, E, I)
p.show()