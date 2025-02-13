from prbm import *

beamradius = 1e-3 # m
properties['I'] = np.pi*beamradius**4/2
properties['A'] = np.pi*beamradius**2
properties['J'] = np.pi*beamradius**4/2


q = PRBM(properties)

q.add_body('A', (0, 0, 0))
q.add_body('B', (0, 0, 0))

n = 3
ra = 20e-3
rb = 10e-3

for i in range(n):
  q.add_flexure('A', (ra*np.cos(i/n*2*np.pi), ra*np.sin(i/n*2*np.pi), 0),
                'B', (rb*np.cos(i/n*2*np.pi), rb*np.sin(i/n*2*np.pi), -.03))

q.show()
q.add_force('B', (2, 1, -5))
q.solve_pose(('B'))
q.show()
