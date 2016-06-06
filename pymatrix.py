import itertools

class Matrix:
  def __init__(s, row, column, modulo=0, data=None):
    if data is None:
      s.data = [[0] * column for x in xrange(row)]
    else:
      s.data = data
    s.row = row
    s.column = column
    s.modulo = modulo

  def __getitem__(s, idx):
    if isinstance(idx, tuple):
      if len(idx) == 1:
        return s.data[idx[0]]
      elif len(idx) == 2:
        return s.data[idx[0]][idx[1]]
    return s.data[idx]

  def __setitem__(s, idx, v):
    if s.modulo != 0:
      v %= s.modulo
    if isinstance(idx, tuple):
      if len(idx) == 1:
        s.data[idx[0]] = v
      elif len(idx) == 2:
        s.data[idx[0]][idx[1]] = v
    else:
      raise RuntimeError("Invalid Argument")

  def __str__(s):
    res = ""
    value_str = [map(str, x) for x in s.data]
    max_len = max([max(map(len, x)) for x in value_str])
    value_str = ["[" + " ".join(map(lambda x: x.ljust(y + 1), x)).rstrip() + "]"
                 for x, y in zip(value_str, itertools.repeat(max_len))]
    return res + "\n" + "\n".join(value_str)

  def __repr__(s):
    return "Matrix(%d, %d, %d, %s)" % (s.row, s.column, s.modulo, s.data)

  def __add__(s, B):
    assert isinstance(B, Matrix)
    assert s.has_same_type(B)
    C = Matrix(s.row, s.column, s.modulo)
    for i in xrange(s.row):
      for j in xrange(s.column):
        C[i, j] = s[i, j] + B[i, j]
    return C

  def __sub__(s, B):
    return s + (-1*B)

  def __rsub__(s, A):
    return A + (-1 * s)

  def __mul__(s, rhs):
    A = s
    B = rhs
    return s.mult(A, B)

  def __rmul__(s, rhs):
    A = rhs
    B = s
    return s.mult(A, B)

  def __div__(s, x):
    C = Matrix(s.row, s.column)
    for i in xrange(s.row):
      for j in xrange(s.column):
        C[i, j] = s[i, j] / x
    return C

  def __pow__(s, n):
    A = s
    B = Matrix(A.row, A.column, A.modulo)
    for i in xrange(A.row):
      B[i, i] = 1
    b = map(int, bin(n)[2:])[::-1]
    for x in b:
      if x:
        B = B * A
      A = A * A
    return B

  def __invert__(s):
    return s.invert()

  def __len__(s):
    return s.row * s.column

  def __eq__(s, B):
    if s.has_same_type(B):
      for i in xrange(s.row):
        for j in xrange(s.column):
          if s[i, j] != B[i, j]:
            return False
      return True
    return False

  def mult(s, A, B):
    if isinstance(A, Matrix):
      if isinstance(B, Matrix):
        C = Matrix(A.row, B.column, s.modulo)
        for i in xrange(A.row):
          for k in xrange(B.row):
            for j in xrange(B.column):
              C[i, j] = C[i, j] + A[i, k] * B[k, j]
        return C
      else:
        C = Matrix(A.row, A.column, s.modulo)
        for i in xrange(A.row):
          for j in xrange(A.column):
            C[i, j] = A[i, j] * B
        return C
    else:
      C = Matrix(B.row, B.column, s.modulo)
      for i in xrange(B.row):
        for j in xrange(B.column):
          C[i, j] = A * B[i, j]
      return C

  def can_mult(s, B):
    if isinstance(B, Matrix):
      return s.row == B.column
    else:
      return True

  def can_invert(s):
    return s.determinant() != 0 and s.is_square()

  def has_same_type(s, B):
    if isinstance(B, Matrix):
      return s.row == B.row and s.column == B.column
    else:
      return False

  def is_square(s):
    return s.row == s.column

  def transpose(s):
    C = Matrix(s.column, s.row, s.modulo)
    for i in xrange(s.column):
      for j in xrange(s.row):
        C[i, j] = s[j, i]
    return C

  def determinant(s):
    assert s.is_square()
    k = s.row
    if k == 1:
      return s[0, 0]
    elif k == 2:
      return s[0, 0] * s[1, 1] - s[0][1] * s[1, 0]
    else:
      for i in xrange(k):
        if s[i, i] == 0:
          h = 0
          for h in xrange(k + 1):
            if h == k:
              break
            if s[i, h] != 0:
              break
          if h == k:
            raise ArithmeticError("Can't calculate determinant of matrix")
          for j in xrange(k):
            s[j, i] += s[j, h]
      for i in xrange(k):
        for j in xrange(k):
          if i < j:
            if s[i, i] == 0:
              raise ArithmeticError("Can't calculate determinant of matrix")
            t = s[i, j] / s[i, i]
            for l in xrange(k):
              s[l, j] -= s[l, i] * t
      det = reduce(lambda x,y:x*y, [s[i, i] for i in xrange(k)], 1)
      return det

  def cofactor(s, I, J):
    C = Matrix(s.row - 1, s.column - 1, s.modulo)
    i = 0
    for x in xrange(s.row):
      if x == I:
        continue
      j = 0
      for y in xrange(s.column):
        if y == J:
          continue
        C[i, j] = s[x, y]
        j += 1
      i += 1
    sign = (-1) ** (I + J)
    return C.determinant() * sign

  def adjugate(s):
    C = Matrix(s.row, s.column, s.modulo)
    for i in xrange(s.row):
      for j in xrange(s.column):
        C[i, j] = s.cofactor(i, j)
    return C

  def invert(s):
    return s.adjugate() / s.determinant()

  @classmethod
  def I(s, n):
    C = Matrix(n, n)
    for i in xrange(n):
      C[i, i] = 1
    return C

  @classmethod
  def O(s, n):
    C = Matrix(n, n)
    return C

def main():
  mat = Matrix(2, 2)
  mat[0, 0] = 1
  mat[0, 1] = 1
  mat[1, 0] = 1
  mat[1, 1] = 0
  r = (mat ** 65537)[0][1]
  print r

  m = mat ** 1024
  print ~m * m == Matrix.I(2)

if __name__ == "__main__":
  main()
