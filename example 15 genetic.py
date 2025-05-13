from prbm import PRBM
from constants import mm, cm
from random import random

tests_per_sample = 5

def f(x):
  p = PRBM(2)
  p.add_body('A', (0, 0))
  p.add_body('B', (0, 3*cm))
  p.add_body('C', (0, 6*cm))
  p.add_body('D', (0, 9*cm))

  p.add_body('E', (0, -cm))
  p.add_body('F', (0, -2*cm))
  p.add_body('G', (0, -3*cm))

  for i in x:
    if i['type'] == 'b':
      p.add_body(i['name'])
    elif i['type'] == 'f':
      p.add_flexure(i['A'], i['pA'], i['B'], i['pB'])

  for j in range(tests_per_sample):
    p.move('B', (random(), random()))
    p.move('C', (random(), random()))
    p.move('D', (random(), random()))

  p.solve_pose('EFG')

  p.show()

f(0)

