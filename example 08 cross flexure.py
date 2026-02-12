from prbm import PRBM

mm = 1e-3 # m
cm = 1e-2 # m

t = 1*mm
h = 3*mm
E = 1650e6 # Pa
A = t*h
I = h*t**3/12

p = PRBM(2)

p.add_body('A', (0, 0))
p.add_body('B', (6*cm, 0))
p.add_flexure('A', (0, -3*cm), 'B', (0, 3*cm))
p.add_flexure('A', (0, 3*cm), 'B', (0, -3*cm))

p.show()

p.add_force('B', (0, .1))
p.solve_pose('B', A, E, I)

p.show()
