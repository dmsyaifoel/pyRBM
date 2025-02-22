from prbm import *

beamsize = 1e-3 # m
properties['I'] = beamsize**4/12
properties['A'] = beamsize**2

l = PRBM(properties)

l.add_body('A', (0, 0))
l.add_body('B', (0, 20e-3))
l.add_body('C', (0, 0))

l.add_flexure('A', (-20e-3, 0), 'B', (-20e-3, 0))
l.add_flexure('A', (20e-3, 0), 'B', (20e-3, 0))
l.add_flexure('B', (-10e-3, 0), 'C', (-10e-3, 0))
l.add_flexure('B', (10e-3, 0), 'C', (10e-3, 0))

# Solve pose
l.move('C', (.005, 0))
l.solve_pose(('B'))
l.show()

# Solve reaction
l.solve_reaction('C', ('B'))

# Solve pose given the same forces (should yield the same pose)
l.moveback()
l.clear_forces()
l.add_force('C', l.reaction[:2])
l.add_torque('C', (l.reaction[2]))
l.solve_pose(('B', 'C'))
l.show()
for body in l.bodies:
  print(body.position)

