from body import Body
from flexure import Flexure

class Topopt:
    def add_flexure(self, bodyA, bodyB):
        flexure = Flexure(bodyA, (0, 0), bodyB, (0, 0), self.E, self.h)
        self.flexures.append(flexure)

        bodyA.flexures.append(flexure)
        bodyA.letters.append('A')
        bodyB.flexures.append(flexure)
        bodyB.letters.append('B')

    def parse_rules(self, s: str):
        s = s.replace(' ', '').split(',')
        k = [i.split('<') for i in s]

        xmin = -1e9
        xmax = 1e9
        ymin = -1e9
        ymax = 1e9

        for l in k:
            if len(l) == 2:
                if l[0] == 'x': xmax = float(l[1])
                if l[1] == 'x': xmin = float(l[0])
            else:
                if l[1] == 'x':
                    xmin = float(l[0])
                    xmax = float(l[2])

            if len(l) == 2:
                if l[0] == 'y': ymax = float(l[1])
                if l[1] == 'y': ymin = float(l[0])
            else:
                if l[1] == 'y':
                    ymin = float(l[0])
                    ymax = float(l[2])

        return ((xmin, xmax), (ymin, ymax))


    def __init__(self, m, n, dx, dy, E, h, t0, fixrules, forcerules, backend):
        self.bodies = []
        for i in range(m):
            for j in range(n):
                self.bodies.append(Body(None, (i * dx, j * dy), 2))



        for i in range(m):
            for j in range(n - 1):
                bodyA = self.bodies[i*n + j]
                bodyB = self.bodies[i*n + j + 1]
                add_flexure(bodyA, bodyB)

        for i in range(m - 1):
            for j in range(n):
                bodyA = self.bodies[i*n + j]
                bodyB = self.bodies[(i + 1)*n + j]
                add_flexure(bodyA, bodyB)

        for i in range(m - 1):
            for j in range(n - 1):
                bodyA = self.bodies[i*n + j]
                bodyB = self.bodies[(i + 1)*n + j + 1]
                add_flexure(bodyA, bodyB)





def clear(self):
    t.fixed = [0 for i in range(len(t.bodies))]
    t.forced = [0 for i in range(len(t.bodies))]





def movet(self, x):
    i = 0
    for j in range(len(self.fixed)):
        if self.fixed[j] == 0:
            self.bodies[j].move((x[i], x[i + 1]), x[i + 2])
            i += 3


def getx(self):
    x = []
    nbodies = len(self.bodies)
    for i in range(nbodies):
        if self.fixed[i] == 0:
            a, b = self.bodies[i].position
            c = self.bodies[i].angles
            x = x + [a, b, c]
    return np.array(x)


def solve_pose_t(self):
    def optfun(x):
        self.movet(x)
        return self.energy()

    self.solt = so.minimize(optfun, self.getx(), method='Nelder-Mead')
    self.movet(self.solt.x)
    print(self.solt)


def forcet(self, s, vector):
    s = s.replace(' ', '').split(',')
    k = [i.split('<') for i in s]

    xmin = -1e9
    xmax = 1e9
    ymin = -1e9
    ymax = 1e9

    for l in k:
        if len(l) == 2:
            if l[0] == 'x': xmax = float(l[1])
            if l[1] == 'x': xmin = float(l[0])
        else:
            if l[1] == 'x':
                xmin = float(l[0])
                xmax = float(l[2])

        if len(l) == 2:
            if l[0] == 'y': ymax = float(l[1])
            if l[1] == 'y': ymin = float(l[0])
        else:
            if l[1] == 'y':
                ymin = float(l[0])
                ymax = float(l[2])

    for i in range(len(t.bodies)):
        pos = t.bodies[i].position
        if xmin < pos[0] < xmax and ymin < pos[1] < ymax:
            t.forced[i] = 1
            self.bodies[i].forces.append(((0, 0), np.array(vector)))


m, n = 3, 3

nbodies = m*n
nflexures = m*(n-1) + (m-1)*n + (m-1)*(n-1)

dim = 2
dim2 = 3

bodies = []
flexures = []

t = symvec('t', nflexures)

A = [h*t_ for t_ in t]
I = [h*t_**3/12 for t_ in t]

for i in range(m):
  for j in range(n):
    bodies.append(Body((i*cm, j*cm)))

def add_flexure(bodyA, bodyB, A, I):
  flexure = Flexure(bodyA, bodyB, A, I)

  flexures.append(flexure)

  bodyA.flexures.append(flexure)
  bodyA.letters.append('A')
  bodyB.flexures.append(flexure)
  bodyB.letters.append('B')

k = 0

for i in range(m):
  for j in range(n - 1):
    bodyA = bodies[i*n + j]
    bodyB = bodies[i*n + j + 1]
    add_flexure(bodyA, bodyB, A[k], I[k])
    k += 1

for i in range(m - 1):
  for j in range(n):
    bodyA = bodies[i*n + j]
    bodyB = bodies[(i + 1)*n + j]
    add_flexure(bodyA, bodyB, A[k], I[k])
    k += 1

for i in range(m - 1):
  for j in range(n - 1):
    bodyA = bodies[i*n + j]
    bodyB = bodies[(i + 1)*n + j + 1]
    add_flexure(bodyA, bodyB, A[k], I[k])
    k += 1