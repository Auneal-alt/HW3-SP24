# note attempt on importing doolittle failed so I made this file to only use cholesky
# my attempt to import doolittle is in testa.py which I didnt have time to fully figure out.
import math


def symmetry(A):
    """ Symmetry tests of the matrix comparing Matrix A to the transposed matrix A to see if they are the same"""

    transpose_A = [[row[i] for row in A] for i in range(len(A[0]))]
    return A == transpose_A


# Used chat gpt to help figure out how to do xtransposeAx>0
def positivedefinite(A):
    """Check if a matrix is positive definite"""
    n = len(A)
    for i in range(n):
        for j in range(i, n):
            A[i][j] = A[i][j] + A[j][i]
    for i in range(n):
        det = 1
        for j in range(i):
            det *= A[i - 1][j]
        if det <= 0:
            return False
    return True


# I followed the equations in the book felt almost easier to solve on paper before trying to do any code
def CholeskyMethod(A, b):
    """Solve the matrix using Cholesky method if it is symmetrical and positive definite"""

    n = len(A)
    L = [[0] * n for _ in range(n)]

    # Fill the lower triangle of L
    for i in range(n):
        s = 0  # Initialize s to 0
        for j in range(i):
            s += L[i][j] * L[i][j]
        L[i][i] = math.sqrt(max(A[i][i] - s, 0))

        for j in range(i + 1, n):
            s = sum(L[j][k] * L[i][k] for k in range(i))
            L[j][i] = (A[j][i] - s) / L[i][i]

    # Solve the lower triangle
    y = [0] * n
    for i in range(n):
        y[i] = b[i][0]  # Extract the values from the column vector
        for j in range(i):
            y[i] -= L[i][j] * y[j]
        y[i] = y[i] / L[i][i]

    # Solve the upper triangle
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= L[j][i] * x[j]
        x[i] = x[i] / L[i][i]

    return x, y


# Important: I couldnt get Doolittle to integrate with my code I worked on it in test.py but kept getting divide by zero
# and other errors, so to show that the Cholesky method would still atleast give an output I took it out.
def main():
    """Main function to solve Cholesky Method"""

    # Problem 1
    A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
    b1 = [[15], [-35], [94], [1]]
    Aaug1 = [[1, -1, 3, 2, 15], [-1, 5, -5, -2, -35], [3, -5, 19, 3, 94], [2, -2, 3, 21, 1]]

    resultsym1 = symmetry(A1)
    print("result of symmetry test is:", resultsym1)

    resultposdef1 = positivedefinite(A1)
    print("result of positivedefinite test is:", resultposdef1)

    MatrixSoln1, y1 = CholeskyMethod(A1, b1)
    # It shouldve used DOolittle but I couldnt figure out a way to import it, my attempt is in the testa.py file
    print("x1 =")
    for i in range(len(MatrixSoln1)):
        print(MatrixSoln1[i])

    # Problem 2
    A2 = [[4, 2, 4, 0], [2, 2, 3, 2], [4, 3, 6, 3], [0, 2, 3, 9]]
    b2 = [[20], [36], [60], [122]]

    resultsym2 = symmetry(A2)
    print("result of symmetry test is:", resultsym2)

    resultposdef2 = positivedefinite(A2)
    print("result of positive definite test is:", resultposdef2)

    MatrixSoln2, y2 = CholeskyMethod(A2, b2)

    print("x2 using Cholesky Method =")
    for i in range(len(MatrixSoln2)):
        print(MatrixSoln2[i])


if __name__ == "__main__":
    main()
