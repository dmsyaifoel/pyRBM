from prbm import PRBM
from constants import mm, cm

t = 1*mm
h = 3*mm
A = t*h
E = 1650e6
I = h*t**3/12

l = PRBM(2)

l.add_body('A', (0, 0))
l.add_body('B', (0, 2*cm))
l.add_body('C', (0, 0))

l.add_flexure('A', (-2*cm, 0), 'B', (-2*cm, 0))
l.add_flexure('A', (2*cm, 0), 'B', (2*cm, 0))
l.add_flexure('B', (-cm, 0), 'C', (-cm, 0))
l.add_flexure('B', (cm, 0), 'C', (cm, 0))

l.show()

l.add_force('C', (1, 2), (0, 0))
l.solve_pose('BC', A, E, I)

l.show()

print(l.bodies['C'].position)