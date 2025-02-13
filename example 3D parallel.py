from prbm import PRBM

enabletorsion = True
finite_difference_delta = .5e-3
optimizer = 'minimize'
gamma = .85
kappa_theta = 2.65

# https://www.oceanz.eu/app/uploads/2021/01/Oceanz_3D_printing_PA12.pdf

properties = {'E': 1650e6, # Pa
              'G': 600e6} # Pa, this is an estimate

beamradius = 1e-3 # m
properties['I'] = np.pi*beamradius**4/2
properties['A'] = np.pi*beamradius**2
properties['J'] = np.pi*beamradius**4/2

p = PRBM(properties)

p.add_body('A', (0, 0, 0))
p.add_body('B', (0, 0, 15e-3))
p.add_body('C', (0, 0, 30e-3))

n = 3
r = 15e-3

for i in range(n):
  p.add_flexure('A', (r*np.cos(i/n*2*np.pi), r*np.sin(i/n*2*np.pi), 0),
                'B', (r*np.cos((i + 1)/n*2*np.pi), r*np.sin((i + 1)/n*2*np.pi), 0))
  p.add_flexure('C', (r*np.cos((i - .5)/n*2*np.pi), r*np.sin((i - .5)/n*2*np.pi), 0),
                'B', (r*np.cos((i + .5)/n*2*np.pi), r*np.sin((i + .5)/n*2*np.pi), 0))

# Solve pose
p.move('C', (.001, .002, .025))
p.solve_pose(('B'))
p.show()
for body in p.bodies:
  print(body.position, body.angles)

# Solve reaction
p.solve_reaction('C', ('B'))

# Solve pose for same reaction, should yield same pose
p.moveback()
p.clear_forces()
p.add_force('C', p.reaction[:3])
p.add_torque('C', p.reaction[3:])
p.solve_pose(('B', 'C'))
p.show()
for body in p.bodies:
  print(body.position, body.angles)
