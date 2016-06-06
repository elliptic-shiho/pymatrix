from timeit import timeit
import copy
import Matrix

N = 4096

mat = Matrix.Matrix(2, 2)
mat[0, 0] = 1
mat[0, 1] = 1
mat[1, 0] = 1
mat[1, 1] = 0
r = (mat ** N)[0][1]
print r

def test():
  X = copy.deepcopy(mat)
  for x in xrange(N - 1):
    X *= mat

def main():
  d = timeit("test()",
             "from test import test", number=10)
  print "%.2f ms/pass" % (d * 1000)

if __name__ == "__main__":
  main()
