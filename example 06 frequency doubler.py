from prbm import PRBM

mm = 1e-3 # m
cm = 1e-2 # m

import numpy as np

t = 1*mm
h = 3*mm
A = t*h
E = 1650e6 # Pa
I = h*t**3/12

s = ((5*mm)**2 + (15*mm)**2)**.5
ds = s - 15*mm
ks = s/.85
ms = ks*.15/2
l = 15*mm
kl = l/.85
ml = kl*.15/2
w = 5*mm

d = PRBM(2)

d.add_body('C', (0, 0))
d.add_body('D', (-cm, 0))

d.add_body('B', (0, -s))

xa0 = np.array((0, -2*s))
d.add_body('A', xa0)

d.add_body('E', (l + ml + 15*mm, 3*cm))
xf0 = np.array((l + ml + 15*mm, 3*cm + ks))
d.add_body('F', xf0)

d.add_flexure('A', (-cm, -ms), 'B', (-cm, ms))
d.add_flexure('A', (cm, -ms), 'B', (cm, ms))

d.add_flexure('B', (0, -ms), 'C', (0, ms))

d.add_flexure('D', (0, 5*mm), 'C', (-5*mm, 0))
d.add_flexure('D', (0, -5*mm), 'C', (-5*mm, 0))

d.add_flexure('C', (-ml, 3*cm), 'E', (-15*mm, 0))

d.add_flexure('E', (-cm, -ms), 'F', (-cm, ms))
d.add_flexure('E', (cm, -ms), 'F', (cm, ms))

d.show()

d.move('A', xa0 + np.array((0, 2*ds)))
d.move('F', xf0 + np.array((-5*mm, -ds)))

d.solve_pose('BCE', A, E, I)

d.show('random')
