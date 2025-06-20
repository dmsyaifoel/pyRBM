import numpy as np
import matplotlib.pyplot as mp
from constants import mm
from prbm import PRBM
cols = mp.rcParams['axes.prop_cycle'].by_key()['color']

scale = 1.2
beamradius = mm*scale # m
I = np.pi*beamradius**4/2
A = np.pi*beamradius**2
E = 1650e6

p = PRBM(3)

initialpos = np.array((0, 0, 29*mm*scale))

p.add_body('A', (0, 0, mm*scale))
p.add_body('B', (0, 0, 15*mm*scale))
p.add_body('C', initialpos)

n = 3
r = 14*mm*scale

for i in range(n):
  p.add_flexure('A', (r*np.cos((i + 1)/n*2*np.pi), r*np.sin((i + 1)/n*2*np.pi), 0),
                'B', (r*np.cos((i + 0)/n*2*np.pi), r*np.sin((i + 0)/n*2*np.pi), 0))
  p.add_flexure('C', (r*np.cos((i + .5)/n*2*np.pi), r*np.sin((i + .5)/n*2*np.pi), 0),
                'B', (r*np.cos((i - .5)/n*2*np.pi), r*np.sin((i - .5)/n*2*np.pi), 0))

d = 2.5*mm

dz = d
dd = .01*mm

p.move('C', initialpos - np.array([0, 0, dz]))
p.solve_pose('B', A, E, I)
e1 = p.energy(A, E, I)

p.move('C', initialpos - np.array([0, 0, dz + dd]))
p.solve_pose('B', A, E, I)
e2 = p.energy(A, E, I)

f1 = abs(e1 - e2)/dd
print(f1)

dx = d/2**.5
dy = d/2**.5

p.move('C', initialpos - np.array([dx, dy, dz]))
p.solve_pose('B', A, E, I)
e1 = p.energy(A, E, I)

p.move('C', initialpos - np.array([dx - dd, dy, dz]))
p.solve_pose('B', A, E, I)
e2 = p.energy(A, E, I)
f2 = abs(e1 - e2)/dd

p.move('C', initialpos - np.array([dx, dy - dd, dz]))
p.solve_pose('B', A, E, I)
e3 = p.energy(A, E, I)
f3 = abs(e1 - e3)/dd

p.move('C', initialpos - np.array([dx, dy, dz - dd]))
p.solve_pose('B', A, E, I)
e4 = p.energy(A, E, I)
f4 = abs(e1 - e4)/dd

data = np.genfromtxt('medium test 11.csv', delimiter = ',', skip_header=1)


n1=4
z = 10*[0] + [i*f1/n1 for i in range(n1)] + 2*[f1]+ [f1 + i*(f4-f1)/n1 for i in range(n1)] + 3*[f4] + [f4 - i*(f4-f1)/n1 for i in range(n1)] + 2*[f1] + [f1 - i*f1/n1 for i in range(n1)] + 5*[0]
# z2 = np.ones_like(z)*f4
y = 16*[0] + [i*f2/n1 for i in range(n1)] + 3*[f2] + [f2 - i*f2/n1 for i in range(n1)] + 10*[0]
x = 16*[0] + [i*f3/n1 for i in range(n1)] + 3*[f3] + [f3 - i*f3/n1 for i in range(n1)] + 10*[0]

w__ = d*1e3/2**.5

z_ = 10*[0] + [i*2.5/n1 for i in range(n1)] + 15*[2.5] + [2.5 - i*2.5/n1 for i in range(n1)] + 5*[0]
y_ = 16*[0] + [i*w__/n1 for i in range(n1)] + 3*[w__] + [w__ - i*w__/n1 for i in range(n1)] + 10*[0]
x_ = 16*[0] + [i*w__/n1 for i in range(n1)] + 3*[w__] + [w__ - i*w__/n1 for i in range(n1)] + 10*[0]


fig, (ax1, ax2) = mp.subplots(
    2, 1,  # 2 rows, 1 column
    sharex=True,
    gridspec_kw={'height_ratios': [1, 10]},  # top 1x, bottom 5x
    figsize=(8, 6)
)

# Top: Displacement
ax1.plot(x_[5:], color=cols[0], label='$\Delta x=\Delta y$')
# ax1.plot(y_[5:], color=cols[1], label='$\Delta y$')
ax1.plot(z_[5:], color=cols[2], label='$\Delta z$')
ax1.set_ylabel('Displ. (mm)')

# Bottom: Force
ax2.plot(x[5:], color=cols[0], label='$F_x$ simulated')
ax2.plot(y[5:], color=cols[1], label='$F_y$ simulated')
ax2.plot(z[5:], color=cols[2], label='$F_z$ simulated')
# ax2.plot(z2[5:], color=cols[2], label='$F_z$ simulated')

ax2.plot((data.T[0]-data.T[1][0])[5:], 'x', label='$F_x$ measured')
ax2.plot((data.T[1]-data.T[1][0])[5:], 'x', label='$F_y$ measured')
ax2.plot((-data.T[2]+data.T[2][0])[5:], 'x', label='$F_z$ measured')
ax2.set_ylabel('Force (N)')
ax2.set_xlabel('Time (s)')

# Legend outside, for force (bottom) axis
ax1.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
ax2.legend(loc='upper left', bbox_to_anchor=(1.01, 1))

fig.tight_layout()
mp.show()

