import numpy as np
import matplotlib.pyplot as mp
from constants import mm
from prbm import PRBM

cols = mp.rcParams['axes.prop_cycle'].by_key()['color']

scale = 1.4
beamradius = mm*scale # m
E = 1650e6
I = np.pi*beamradius**4/2
A = np.pi*beamradius**2

p = PRBM(3)

initialpos = np.array((0, 0, 29*mm*scale))

p.add_body('A', (0, 0, mm*scale))
p.add_body('B', (0, 0, 15*mm*scale))
p.add_body('C', initialpos)

n = 3
r = 14*mm*scale

for i in range(n):
  p.add_flexure('A', (r*np.cos(i/n*2*np.pi), r*np.sin(i/n*2*np.pi), 0),
                'B', (r*np.cos((i + 1)/n*2*np.pi), r*np.sin((i + 1)/n*2*np.pi), 0))
  p.add_flexure('C', (r*np.cos((i - .5)/n*2*np.pi), r*np.sin((i - .5)/n*2*np.pi), 0),
                'B', (r*np.cos((i + .5)/n*2*np.pi), r*np.sin((i + .5)/n*2*np.pi), 0))

# p.show()

p.move('C', initialpos - np.array([0, 0, 5*mm]))
p.solve_pose('B', A, E, I)
e1 = p.energy(A, E, I)

p.move('C', initialpos - np.array([0, 0, 4.9*mm]))
p.solve_pose('B', A, E, I)
e2 = p.energy(A, E, I)

f = abs(e1 - e2)/(.1*mm)
print(f)

data = np.genfromtxt('large test 4.csv', delimiter =',', skip_header=1)

n = 10
dat_ = [0]*18 + [f*((i+1)/n) for i in range(n)] + [f]*19 + [f - f*(i+1)/(n-1) for i in range(n-1)] + [0]*5

mp.plot(dat_[5:], label='Simulation')

n_ = 325
mp.plot(np.arange(len(data.T[i])-n_)/62.5, -data.T[2][n_:], '--', label='Measurements', color=cols[0])

mp.xlabel('Time (s)')
mp.ylabel('Force (N)')
mp.legend()
mp.show()