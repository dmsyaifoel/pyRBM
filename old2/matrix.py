class Matrix:
  def __init__(self, l):
    if isinstance(l, Matrix):
      l = l.l
    assert isinstance(l, (list, tuple))
    if not isinstance(l[0], (list, tuple)):
      l = [[i] for i in l]
    self.l = l
    self.dim = (len(l), len(l[0]))

  def __repr__(self):
    return str(self.l)

  def __getitem__(self, i):
    if self.dim[0] == 1:
      return self.l[i][0]
    if self.dim[1] == 1:
      return self.l[0][i]
    return self.l[i]

  def __add__(self, other):
    # assert self.dim == other.dim
    return Matrix([[self.l[i][j] + other.l[i][j] for j in range(self.dim[1])] for i in range(self.dim[0])])

  def __mul__(self, other):
    return Matrix([[other*self.l[i][j]for j in range(self.dim[1])] for i in range(self.dim[0])])

  def __rmul__(self, other):
    return self*other

  def __sub__(self, other):
    return self + ((-1)*other)

  def __truediv__(self, other):
    return self*(other**(-1))

  def dot(self, other):
    return (self.T@other)[0]

  @property
  def norm(self):
    return sum([sum([self.l[i][j]**2 for j in range(self.dim[1])]) for i in range(self.dim[0])])**.5

  def normalize(self):
    return self.__rmul__(1/self.norm)

  @property
  def T(self):
    return Matrix([[self.l[i][j] for i in range(self.dim[0])] for j in range(self.dim[1])])

  def __matmul__(self, other):
    assert self.dim[1] == other.dim[0]
    return Matrix([[sum([self.l[i][k]*other.l[k][j] for k in range(self.dim[1])]) for j in range(other.dim[1])]for i in range(self.dim[0])])