# Note: I used the in class solution to rewrite GaussSeidel.py since I originally used numpy and needed to change it

import copy


def GaussSeidel(Aaug, x, Niter=15, epsilon=1e-5):
    """Gauss-Seidel iterative to a set of equations in an augmented matrix by  ensuring diagonal dominant,
    solving each row in sequence from x[0} to x[n-1], then repeating for n iterations(Niter) until max change
    is less than epsilon
    :param Aaug: an augmented matrix
    :param x: initial guess vector
    :param Niter: number of iterations to get correct x
    :param epsilon: precision for early escape
    :return: x solution vector"""

    AA = copy.deepcopy(Aaug)  # deepcopy Aaug so we dont alter it

    AA = DiagDominant(AA)  # ensure diagonal dominance

    n = len(x)
    for j in range(Niter):
        maxErr = 0
        for r in range(n):
            xOld = x[r]
            rhs = AA[r][n]
            for c in range(n):
                if c != r:
                    rhs -= AA[r][c] * x[c]
            x[r] = rhs / AA[r][r]
            maxErr = max(maxErr, abs(xOld - x[r]))
        if maxErr <= epsilon:
            break
    return x


def DiagDominant(A):
    """Function to make matrix Diagonal Dominant
    :param A: matrix"""

    AA = copy.deepcopy(A)
    rows = len(AA)
    for i in range(rows):
        c = abs(AA[i][i])
        for k in range(i + 1, rows):
            if abs(AA[k][i]) > c:
                row = AA.pop(k)  # pop and insert move rows in a matrix
                AA.insert(i, row)
                c = abs(AA[i][i])
    return AA


def separateAugmented(Aaug):
    """
    This function separates the last column from Aaug and returns a group of matrices A, b
    :param Aaug: the augmented matrix
    :return: (A,b)
    """
    A = copy.deepcopy(Aaug)
    b = []
    n = len(A[0]) - 1
    for r in A:
        b.append(r.pop(n))
    return (A, b)


def checkMatrixSoln(A, x, augmented=True):
    """double checks the matrix is true by multiplying in solution vector
    :param augmented: The augmented matrix
    :param x: The solution vector transpose (a row vector)
    :return: The b vector transpose (a row vector)
    """
    if augmented:  # If A is augmented, strip off the last column
        AA, b = separateAugmented(A)
    else:
        AA = A
    B = []  # the result of multiplying AA*x
    for r in AA:
        s = 0
        rCntr = 0
        for c in r:
            s += c * x[rCntr]
            rCntr += 1
        B.append(s)

    rounded_B = []
    for num in B:
        rounded_B.append(round(num, 1))
    return rounded_B


def matrixMult(A,B):
    """
    Multiplies matrix A by Matrix B.  Note:  you are responsible for making sure nxm*mxp gives nxp
    :param A: a matrix
    :param B: another matrix
    :return: the matrix with size nxp
    """
    ARows=len(A)
    BCols=len(B[0])
    C=[[0 for c in range(BCols)] for r in range(ARows)]
    for r in range(ARows):
        for c in range(BCols):
            C[r][c]=multVecs(A[r],getCol(B,c))
    return C


def getCol(A, c):
    """
    Gets a column from A matrix
    :param A: a matrix
    :param c: the index for column you want
    :return: the row vector corresponding to column c
    """
    vec=[]
    for r in A:
        vec.append(r[c])
    return vec


def multVecs(A,B):
    """
    simply multiplies the vectors.  Note:  they should both be row vectors of same length
    :param A: a row vector
    :param B: another row vector
    :return: the product of the multiplication, a scalar
    """
    s=0
    for a in range(len(A)):
        s+=A[a]*B[a]
    return s


def main():
    """main function to solve problem"""

    A1 = [[3, 1, -1, 2], [1, 4, 1, 12], [2, 1, 2, 10]]

    x1 = [0, 0, 0]

    xSoln1 = GaussSeidel(A1, x1, Niter=22)
    print("xsoln1=", [round(x, 4) for x in xSoln1])


if __name__ == "__main__":
    main()
