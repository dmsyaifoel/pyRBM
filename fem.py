from body import Body
from flexure import Flexure

class Fem:
  '''
  Class containing useful functions for setting up a finite element model based on PRBM.

  These functions are just for convenience: it might be easier to set up something manual.
  '''
  def __init__(self, m, n, dx, dy, E, h):
    self.bodies = []
    self.flexures = []
    self.h = h
    self.E = E

    for i in range(m):
      for j in range(n):
        self.bodies.append(Body((i*dx, j*dy)))

    for i in range(m):
      for j in range(n - 1):
        bodyA = self.bodies[i*n + j]
        bodyB = self.bodies[i*n + j + 1]
        self.add_flexure(bodyA, bodyB)

    for i in range(m - 1):
      for j in range(n):
        bodyA = self.bodies[i*n + j]
        bodyB = self.bodies[(i + 1)*n + j]
        self.add_flexure(bodyA, bodyB)

    for i in range(m - 1):
      for j in range(n - 1):
        bodyA = self.bodies[i*n + j]
        bodyB = self.bodies[(i + 1)*n + j + 1]
        self.add_flexure(bodyA, bodyB)

    self.clear()

  def clear(self):
    self.fixed = []
    self.forced = []
    for body in self.bodies:
      body.force = array((0, 0))
    self.targeted = []

  def show(self, showflexures=False, tlist=None, tmin=None, tmax=None):
    import numpy as np
    import matplotlib.pyplot as mp
    cols = mp.rcParams['axes.prop_cycle'].by_key()['color']
    points = np.array([[body.position[0], body.position[1]] for body in self.bodies]).T
    mp.scatter(points[0], points[1], color=cols[0])

    if showflexures:
      for i, flexure in enumerate(self.flexures):
        a, b, c, d = flexure.line()
        if tlist is None:
          mp.plot((a[0], b[0], c[0], d[0]), (a[1], b[1], c[1], d[1]), color=cols[0])
        else:
          mp.plot((a[0], b[0], c[0], d[0]), (a[1], b[1], c[1], d[1]), color=cols[0], alpha=(tlist[i] - tmin)/(tmax - tmin))


    for i in self.fixed:
      pos = self.bodies[i].position
      mp.plot(pos[0], pos[1], 'x', color=cols[1])

    for i in self.forced:
      pos = self.bodies[i].position
      mp.plot(pos[0], pos[1], 's', color=cols[2])

    for i in self.targeted:
      pos = self.bodies[i].position
      mp.plot(pos[0], pos[1], '^', color=cols[3])

    mp.axis('equal')
    mp.show()

  def add_flexure(self, bodyA, bodyB):
    flexure = Flexure(bodyA, bodyB)
    self.flexures.append(flexure)
    bodyA.flexures.append(flexure)
    bodyA.which.append(True)
    bodyB.flexures.append(flexure)
    bodyB.which.append(False)

  def parse_rules(self, rules):
    rules = rules.replace(' ', '').replace('*','').replace('mm','e-3').replace('cm', 'e-2').split(';')
    rulelist = []
    for rule in rules:
      s = rule.split(',')
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

      rulelist.append(((xmin, xmax), (ymin, ymax)))
    return rulelist

  def set_attribute(self, attribute, rules, vector=None):
    rulelist = self.parse_rules(rules)
    for i, body in enumerate(self.bodies):
      for rule in rulelist:
        applies = True
        for j in range(2):
          if not (rule[j][0] < body.position[j] < rule[j][1]):
            applies = False
        if applies:
          if 'fix' in attribute:
            self.fixed.append(i)
          if 'move' in attribute:
            body.move(body.position0 + matrix(vector))
          if 'force' in attribute:
            self.forced.append(i)
            body.force = matrix(vector)
          if 'target' in attribute:
            self.targeted.append(i)
            body.target = body.position0 + matrix(vector)

  def ndof(self):
    return 3*(len(self.bodies) - len(self.fixed))

  def nflex(self):
    return len(self.flexures)

  def energy(self, tlist):
    fen = sum([flexure.energy(self.h*tlist[i], self.E, self.h*tlist[i]**3/12) for i, flexure in enumerate(self.flexures)])
    ben = sum([body.energy() for body in self.bodies])
    return fen + ben

  # def targetloss(self):
  #   loss = 0
  #   for i in self.targeted:
  #     loss += dist(self.bodies[i].position, self.bodies[i].target)**2
  #   return loss