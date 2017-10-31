from sympy.core.symbol import Symbol
from sympy import simplify
from sympy.core.numbers import Float


class matrix:
    """
    A class for matrices, quite a good one at that if I may add
    """
    # usually nxm matrix n = no lines, m number of columns
    n = None
    m = None
    A = None
    type = None

    # Built in functions
    def __init__(self, *args):
        """
        If only input argument is an integer then the return is the identity matrix with
        dim = input int.
        input parameters are either a variable composed of an array of arrays, each array
        within the array must be of the same length. With two input parameters they must
        be integers first being the number of lines and the next being the number of
        columns, it will create a matrix of zeros
        """
        if len(args) == 1 and isinstance(args[0], int):
            self.n = args[0]
            self.m = args[0]
            self.A = [[0 for i in range(args[0])] for j in range(args[0])]
            self.type = 'n'
            for i in range(args[0]):
                self.A[i][i] = 1
        elif len(args) == 1 and isinstance(args[0], list) and not isinstance(args[0][0], list):
            self.n = 1
            self.m = len(args[0])
            self.A = [args[0]]
        elif len(args) == 1 and not isinstance(args[0],int):
            assert(isinstance(args[0], list))
            assert(isinstance(args[0][0], list))
            if not isinstance(args[0][0][0], float) and not isinstance(args[0][0][0], int): self.type = 's'
            else: self.type = 'n'
            self.n = len(args[0])
            self.m = len(args[0][0])
            self.A = args[0]
        elif len(args) == 2:
            n = args[0]
            m = args[1]
            self.type = 'n'
            self.n = n
            self.m = m
            self.A = [[0 for i in range(m)] for j in range(n)]

    def __add__(self, other):
        assert(other.n == self.n)
        assert(other.m == self.m)
        B = matrix(other.n, other.m)
        for i in range(B.n):
            for j in range(B.m):
                B.A[i][j] = self.A[i][j] + other.A[i][j]
        return B

    def __sub__(self, other):
        assert (other.n == self.n)
        assert (other.m == self.m)
        B = matrix(other.n, other.m)
        for i in range(B.n):
            for j in range(B.m):
                B.A[i][j] = self.A[i][j] - other.A[i][j]
        return B

    def __mul__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            A = self.copy()
            for i in range(A.n):
                for j in range(A.m):
                    A[i][j] = other*A[i][j]
            return A
        else:
            assert(self.m == other.n)
            mat1 = self.A
            mat2 = other.A
            mat3 = matrix(self.n, other.m)
            for i in range(len(mat1)):
                if isinstance(mat2[i], list):
                    for j in range(len(mat2[i])):
                        temp = 0
                        for k in range(len(mat1[i])):
                            temp += mat1[i][k] * mat2[k][j]
                        mat3.A[i][j] = temp
            if mat3.n == 1 and mat3.m == 1:
                mat3 = mat3[0][0]
            return mat3

    def __rmul__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            A = self.copy()
            for i in range(A.n):
                for j in range(A.m):
                    A[i][j] = other*A[i][j]
            return A

    def __repr__(self):
        str = "-\n"
        for i in range(self.n):
            for j in range(self.m):
                if j == 0:
                    str += "|"
                str += "{}".format(self.A[i][j])
                if j != self.m - 1:
                    str += "\t"
            str += "\n"
            if i != self.n - 1:
                str += "|\n"
        str += "-"
        return str

    def __str__(self):
        try:
            str = "-\n"
            for i in range(self.n):
                for j in range(self.m):
                    if j == 0:
                        str += "|"
                    str += "{}".format(self.A[i][j])
                    if j != self.m-1:
                        str+="\t"
                str += "\n"
                if i != self.n-1:
                    str += "|\n"
            str += "-"
            return str
        except IndexError:
            str = "["
            for i in range(len(self.A[0])):
                str += "{}".format(self.A[0][i])
                if i != len(self.A[0])-1:
                    str += "\t"
            str += "]"
            return str

    def __pow__(self, power, modulo=None):
        if not isinstance(power, int):
            raise SystemError("Matrix powers must be integers greater than or equal to -1")
        if power == -1:
            return self.inverse2()
        if power < -1:
            raise SystemError("negative powers not defined for matrices, -1 only implies the matrix's inverse")
        else:
            temp = matrix([[j for j in self.A[i]] for i in range(len(self.A))])
            B = matrix([[j for j in self.A[i]] for i in range(len(self.A))])
            for i in range(power-1):
                B *= temp
            return B

    def __truediv__(self, other):
        return self**-1*other

    def __eq__(self, other):
        if self.m != other.m or self.n != other.n:
            return False
        for i in range(self.n):
            for j in range(self.m):
                if self.A[i][j] != other.A[i][j]:
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, item):
        return self.A[item]

    # Matrix operations
    def inverse(self):
        """
        Matrix must be square!
        :return: the inverse of the matrix, if it exists
        """
        mat = self.copy()
        I = matrix(mat.n)
        mat = mat.A
        I = I.A
        for i in range(len(mat)):
            I[i] = self.divideRow(I[i], mat[i][i])
            mat[i] = self.divideRow(mat[i], mat[i][i])
            for j in range(i + 1, len(mat)):
                if mat[j][i] != 0:
                    I[j] = self.divideRow(I[j], mat[j][i])
                    mat[j] = self.divideRow(mat[j], mat[j][i])
                    I[j] = self.subRow(I[i], I[j])
                    mat[j] = self.subRow(mat[i], mat[j])
        for i in reversed(range(len(mat))):
            I[i] = self.divideRow(I[i], mat[i][i])
            mat[i] = self.divideRow(mat[i], mat[i][i])
            for j in reversed(range(len(mat[j]) - 1)):
                if mat[j][i] != 0 and i != j:
                    I[j] = self.divideRow(I[j], mat[j][i])
                    mat[j] = self.divideRow(mat[j], mat[j][i])
                    I[j] = self.subRow(I[i], I[j])
                    mat[j] = self.subRow(mat[i], mat[j])
        return matrix(I)

    def inverse2(self):
        A = self.copy().A

        n = len(A)
        m = len(A[0])
        I = matrix(n).A
        ord = []
        for i in range(n):
            j = 0
            for k in range(m):
                if A[i][k] != 0:
                    j = k
                    ord.append(j)
                    break
            div = A[i][j]
            A[i] = self.divideRow(A[i], div)
            I[i] = self.divideRow(I[i], div)
            for k in range(i+1, n):
                if A[k][j] != 0:
                    I[k] = self.divideRow(I[k], A[k][j])
                    A[k] = self.divideRow(A[k], A[k][j])
                    I[k] = self.subRow(I[i], I[k])
                    A[k] = self.subRow(A[i], A[k])
        A = matrix(A)
        temp = A.copy().A
        Itemp = matrix(I).copy().A
        A = A.A
        for i in range(n):
            I[ord[i]] = [j for j in Itemp[i]]
            A[ord[i]] = [j for j in temp[i]]

        for i in reversed(range(len(A))):
            I[i] = self.divideRow(I[i], A[i][i])
            A[i] = self.divideRow(A[i], A[i][i])
            for j in reversed(range(len(A[j]) - 1)):
                if A[j][i] != 0 and i != j:
                    I[j] = self.divideRow(I[j], A[j][i])
                    A[j] = self.divideRow(A[j], A[j][i])
                    I[j] = self.subRow(I[i], I[j])
                    A[j] = self.subRow(A[i], A[j])

        return matrix(I)

    def reduce(self):
        """
        returns the row reduced echelon form of the matrix
        """
        A = self.copy().A

        for i in range(min(len(A),len(A[0]))):
            if A[i][i] != 0:
                A[i] = self.divideRow(A[i], A[i][i])
                for j in range(i+1, len(A)):
                    if A[j][i] != 0:
                        A[j] = self.divideRow(A[j],A[j][i])
                        A[j] = self.subRow(A[i],A[j])
        for i in reversed(range(min(len(A[0]),len(A)))):
            if A[i][i] != 0:
                A[i] = self.divideRow(A[i],A[i][i])
                for j in reversed(range(i)):
                    if A[j][i] != 0:
                        A[j] = self.divideRow(A[j], A[j][i])
                        A[j] = self.subRow(A[i], A[j])
        A = matrix(A).smooth()
        return A

    def t(self):
        """
        Transposes a matrix
        :return: the matrix's transpose
        """
        n = self.n
        m = self.m
        A = matrix(m,n)
        try:
            for i in range(n):
                for j in range(m):
                    A[j][i] = self[i][j]
        except IndexError:
            return matrix([[i] for i in self.A[0]])
        return A

    def det(self):
        """
        Calculates the determinant of a matrix
        """
        if not self.isSquare():
            raise SystemError("Determinants can only be found for square matrices")
        if self.size() == [2, 2]:
            return self[0][0]*self[1][1] - self[0][1]*self[1][0]
        else:
            det = 0
            for i in range(self.n):
                A = self.copy()
                A = A.delRow(0)
                A = A.delCol(i)
                if isinstance(self[0][i], float) or isinstance(self[0][i], int):
                    if abs(self[0][i]) > 10**-5:
                        if i%2 == 0:
                            det += self[0][i]*A.det()
                        else:
                            det -= self[0][i]*A.det()
                else:
                    if i%2 == 0:
                        det += self[0][i]*A.det()
                    else:
                        det -= self[0][i]*A.det()
            return det

    def adj(self):
        A = self.copy()
        for i in range(A.n):
            for j in range(A.m):
                temp = A.delCol(j).delRow(i)
                mult = 1
                if i%2 == 0:
                    if j%2 != 0:
                        mult = -1
                if i%2 == 1:
                    if j%2 == 0:
                        mult = -1
                A[i][j] = mult*temp.det()
        return A

    def solve(self, b):
        """
        If  Ax=b, where A is the matrix in question then this
        functions solves for x
        :param b: the solution vector
        :return: x = A^-1b
        """
        return self**-1*b

    def diff(self):
        A = self.copy()
        for i in range(self.n):
            for j in range(self.m):
                A[i][j] = self[i][j].diff()
        return A

    # Matrix information
    def rank(self):
        A = self.reduce()
        rank = 0
        for i in A:
            if sum([abs(j) for j in i],0) > 0:
                rank += 1
        return rank

    def size(self):
        return [self.n, self.m]

    # Boolean return functions
    def linind(self):
        """
        Boolean, checks if matrix is linearly independent
        """
        return self.rank() == self.n

    def isSquare(self):
        return self.n == self.m

    # Miscellaneous and other helper functions
    def switchRow(self, i, j):
        A = self.copy()
        temp = A.A[i]
        A.A[i] = A.A[j]
        A.A[j] = temp
        return A

    def switchCol(self, i, j):
        A = self.copy().t()
        A = A.switchRow(i,j)
        return A.t()

    def delRow(self,i):
        A = self.copy()
        return matrix(A.A[0:i]+A.A[i+1:])

    def delCol(self,i):
        A = self.copy().t()
        A = A.delRow(i)
        return A.t()

    def smooth(self, div=1000):
        """
        Most of the time, numbers that are supposed to be zero are extremely close, but not
        quite, to make those numbers 0 we use this function
        :param div: accuracy-> the higher the more accurate
        :return: A more realistic matrix
        """
        A = self.copy()
        ave = sum([abs(sum(i,0)) for i in A.A], 0)/(A.n*A.m)
        for i in range(len(A.A)):
            for j in range(len(A.A[i])):
                if abs(A.A[i][j]) < ave/div:
                    A.A[i][j] = 0
        return A

    def copy(self):
        return matrix([[j for j in self.A[i]] for i in range(len(self.A))])

    def divideRow(self, row, num):
        for i in range(len(row)):
            row[i] = row[i] / num
        return row

    def multiplyRow(self, row, num):
        for i in range(len(row)):
            row[i] = row[i] * num
        return row

    def subRow(self, row1, row2):
        for i in range(len(row1)):
            row2[i] = row2[i] - row1[i]
        return row2

    def addRow(self, row1, row2):
        for i in range(len(row1)):
            row2[i] = row2[i] + row1[i]
        return row2

    def appendCol(self, row):
        if isinstance(row, matrix):
            temp = []
            for i in row.A:
                temp.append(i[0])
            row = temp
        A = self.copy().t().A
        A.append(row)
        A = matrix(A).t()
        return A

    def appendRow(self, row):
        if isinstance(row, matrix):
            temp = []
            for i in row.A:
                temp.append(i[0])
            row = temp
        A = self.copy().A
        A.append(row)
        A = matrix(A)
        return A

    def sub(self,*args):
        A = self.copy()
        for i in range(self.n):
            for j in range(self.m):
                for k in range(len(args)):
                    if k == 0 or k%2 == 0:
                        A[i][j] = A[i][j].subs(args[k],args[k+1])
        for i in range(self.n):
            for j in range(self.m):
                A[i][j] = float(A[i][j])
        return A

    def simplify(self):
        A = self.copy()
        for i in range(A.n):
            for j in range(A.m):
                A[i][j] = self.pre(A[i][j])
        for i in range(self.n):
            for j in range(self.m):
                    A[i][j] = simplify(A[i][j])
        for i in range(A.n):
            for j in range(A.m):
                A[i][j] = self.pre(A[i][j])
        return A

    def pre(self, B):
        if isinstance(B, Float):
            if abs(B) - 1 < 0.0001:
                B = int(B)
            else:
                B = round(B, 3)
        else:
            try:
                for arg in B.args:
                    if isinstance(arg, Float):
                        if abs(arg) - 1 < 0.0001:
                            B = B.subs(arg, int(arg))
                        else:
                            B = B.subs(arg, round(arg, 3))
                    try:
                        if arg.args:
                            arg2 = self.pre(arg)
                            B = B.subs(arg, arg2)
                    except AttributeError:
                        pass
            except AttributeError:
                pass
        return B

    def index(self, k):
        for i in range(self.n):
            for j in range(self.m):
                if self[i][j] == k:
                    return [j, i]
        return None

    def rindex(self, k):
        for i in reversed(range(self.n)):
            for j in reversed(range(self.m)):
                if self[i][j] == k:
                    return [j, i]

    def mat2mat(self):
        str = "["
        for i in range(self.n):
            for j in range(self.m):
                str += " {}".format(self[i][j])
                if j != self.m-1:
                    str += ", "
                else:
                    str += "; "
        str += "]"
        return str

    def tsum(self):
        tot = 0
        for i in range(self.n):
            for j in range(self.m):
                tot += self[i][j]
        return tot
def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

if __name__ == "__main__":
    A = matrix(5)
    print(A)
    A = matrix([[10, 2, 1, 4], [1, 2, 2, 4], [1, 2, 4, 4], [1,2,4,3]])
    b = A.delRow(1)
    b = b.delCol(1)
    b = matrix([[1,2,3,4]]).t()



