from fem import Fem

mm = 1e-3
cm = 1e-2

m, n = 3, 3
dx, dy = cm, cm
E = 1650e6
h = 2*mm
t = .5*mm
A = h*t
I = h*t**3/12

f = Fem(m, n, dx, dy, E, h, t)
f.fix('x<1mm')
# f.force('1cm<x,.5mm<y', (1, 1))
# f.show(1)
print(f.energy(A, E, I))