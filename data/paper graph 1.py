n__ = 9
n___ = 12

from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, rcParams
import numpy as np

cols = rcParams['axes.prop_cycle'].by_key()['color']


data12 = np.genfromtxt('medium test 1.csv', delimiter=' , ')
sim12 = np.load('medium.npy')

sim1 = np.load('small.npy')
data1b = np.genfromtxt('small test 2.csv', delimiter=' , ')

data14 = np.genfromtxt('large test 1.csv', delimiter=' , ')
sim14 = np.load('large.npy')

data1h = np.genfromtxt('hollow test 1.csv', delimiter=' , ')

plot(np.arange(n__)/2,-sim1.T[2][:n__], color=cols[0], label='Scale 1.0 simulation')
plot(np.arange(n__)/2,-data1b.T[2][:n__]+data1b.T[2][0], 'x', color=cols[0], label='Scale 1.0 measurements')
plot(np.arange(n__)/2,-data1h.T[2][:n__]+data1h.T[2][0], 'o', mfc='none', color=cols[0], label='Scale 1.0 hollow measurements')



plot(np.arange(n___)/2, -sim12.T[2][:n___], color=cols[1], label='Scale 1.2 simulation')
plot(np.arange(n___)/2, -data12.T[2][:n___]+data12.T[2][0], 'x', color=cols[1], label='Scale 1.2 measurements')

plot(np.arange(len(sim14.T[2]))/2,-sim14.T[2], color=cols[2], label='Scale 1.4 simulation')
plot(np.arange(len(data14.T[2]))/2,-data14.T[2]+data14.T[2][0], 'x', color=cols[2], label='Scale 1.4 measurements')

print(-data1b.T[2][:n__]/-sim1.T[2][:n__])
print(-data12.T[2][:n__]/-sim12.T[2][:n__])
print(-data14.T[2][:n__]/-sim14.T[2][:n__])

l = legend()
l.get_frame().set_linewidth(0.0)
l.get_frame().set_facecolor('none')
xlabel('$u$ (mm)')
ylabel('$F$ (N)')
show()