from body import Body
from flexure import Flexure

class Fem:
  def __init__(self, m, n, dx, dy, E, h, t):
    self.bodies = []
    self.flexures = []

    for i in range(m):
      for j in range(n):
        self.bodies.append(Body((i * dx, j * dy)))

    for i in range(m):
      for j in range(n - 1):
        bodyA = self.bodies[i * n + j]
        bodyB = self.bodies[i * n + j + 1]
        self.add_flexure(bodyA, bodyB)

    for i in range(m - 1):
      for j in range(n):
        bodyA = self.bodies[i * n + j]
        bodyB = self.bodies[(i + 1) * n + j]
        self.add_flexure(bodyA, bodyB)

    for i in range(m - 1):
      for j in range(n - 1):
        bodyA = self.bodies[i * n + j]
        bodyB = self.bodies[(i + 1) * n + j + 1]
        self.add_flexure(bodyA, bodyB)

    self.clear()

  def clear(self):
    self.fixed = len(self.bodies)*[False]
    self.forced = len(self.bodies)*[False]
    for body in self.bodies:
      body.clear_force()

  def show(self, showflexures=False):
    import numpy as np
    import matplotlib.pyplot as mp
    points = np.array([body.position for body in self.bodies])
    mp.scatter(points.T[0], points.T[1], color='blue')

    if showflexures:
      for flexure in self.flexures:
        line = np.array(flexure.line).T
        mp.plot(line[0], line[1], color='blue')

    for i, isfixed in enumerate(self.fixed):
      if isfixed:
        pos = self.bodies[i].position
        mp.plot(pos[0], pos[1], 'x', color='red')

    for i, isforced in enumerate(self.forced):
      if isforced:
        pos = self.bodies[i].position
        mp.plot(pos[0], pos[1], 's', color='green')

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
    rules = rules.replace(' ', '').replace('*',' ').replace('mm','e-3').replace('cm', 'e-2').split(';')
    ruleslist = []
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

      ruleslist.append(((xmin, xmax), (ymin, ymax)))
    return ruleslist

  def fix(self, rules, move=False):
    rulelist = self.parse_rules(rules)
    for i, body in enumerate(self.bodies):
      applies = True
      for rule in rulelist:
        for j in range(2):
          if not (rule[j][0] < body.position[j] < rule[j][1]):
            applies = False
        if applies:
          self.fixed[i] = True
          if move:
            body.move(move)

  def force(self, rules, force):
    rulelist = self.parse_rules(rules)
    for i, body in enumerate(self.bodies):
      applies = True
      for rule in rulelist:
        for j in range(2):
          if not (rule[j][0] < body.position[j] < rule[j][1]):
            applies = False
        if applies:
          self.forced[i] = True
          body.force = force

  def energy(self, A, E, I):
    return sum([body.energy for body in self.bodies]) + sum([flexure.energy(A, E, I) for flexure in self.flexures])