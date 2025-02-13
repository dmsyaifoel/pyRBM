from prbm import PRBM

enabletorsion = True
finite_difference_delta = .5e-3
optimizer = 'minimize'
gamma = .85
kappa_theta = 2.65

# https://www.oceanz.eu/app/uploads/2021/01/Oceanz_3D_printing_PA12.pdf
properties = {'E': 1650e6, # Pa
              'G': 600e6} # Pa, this is an estimate

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
