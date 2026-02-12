from prbm import PRBM

mm = 1e-3 # m
cm = 1e-2 # m

t = 1*mm
h = 3*mm
E = 1650e6 # Pa
A = t*h
I = h*t**3/12

g = PRBM(2)

g.add_body('A', (2*cm, 0))

g.add_body('B', (0, 3*cm))
g.add_body('C', (0, 6*cm))
g.add_body('D', (0, 9*cm))
g.add_body('E', (0, 12*cm))


g.add_body('F', (3*cm, 3*cm))
g.add_body('G', (2*cm, 6*cm))
g.add_body('H', (1*cm, 9*cm))


g.add_flexure('A', (-2*cm, 0), 'B', (0, 0))
g.add_flexure('B', (0, 0), 'C', (0, 0))
g.add_flexure('C', (0, 0), 'D', (0, 0))
g.add_flexure('D', (0, 0), 'E', (0, 0))

g.add_flexure('A', (2*cm, 0), 'F', (0, 0))
g.add_flexure('F', (0, 0), 'G', (0, 0))
g.add_flexure('G', (0, 0), 'H', (0, 0))

g.add_flexure('B', (0, 0), 'F', (0, 0))
g.add_flexure('C', (0, 0), 'G', (0, 0))
g.add_flexure('D', (0, 0), 'H', (0, 0))
g.add_flexure('E', (0, 0), 'H', (0, 0))
g.show()

g.add_force('B', (1, 0))
g.add_force('C', (2, 0))
g.add_force('D', (2, 0))
g.add_force('E', (1, 0))
g.solve_pose('BCDEFGH', A, E, I)
g.show()
