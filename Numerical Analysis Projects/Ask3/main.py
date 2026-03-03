import numpy as np
import math


def lu(A):
    n = A.shape[0]
    P = np.identity(n, int)
    L = np.zeros((n, n))
    for k in range(n-1):
        # Partial pivoting
        maxA = abs(A[k, k])
        imax = k
        for i in range(k, n):
            if abs(A[i, k]) > maxA:
                maxA = abs(A[i, k])
                imax = i
        A[[k, imax]] = A[[imax, k]]  # swap row k with row imax in A
        P[[k, imax]] = P[[imax, k]]  # swap row k with row imax in P
        L[[k, imax]] = L[[imax, k]]  # swap row k with row imax in L
        pivot = A[k, k]
        if pivot == 0:
            continue

        # LU decomposition
        for i in range(k + 1, n):
            c = A[i][k] / pivot
            L[i, k] = c
            for j in range(k, n):
                A[i, j] -= c * A[k, j]
    for i in range(n):
        L[i, i] += 1
    U = A
    return [L, U, P]


def fsub(L, b):
    for i in range(L.shape[0]):
        for j in range(i):
            b[i] -= L[i, j] * b[j]
    y = b
    return y


def bsub(U, y):
    n = U.shape[0]
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            y[i] -= U[i, j] * y[j]
        y[i] = y[i] / U[i, i]
    x = y
    return x


def gauss(A, b):
    L, U, P = lu(A)  # PA = LU
    print('PA = LU')
    print('P = \n', P)

    print('L = \n', L)

    print('U = \n', U, '\n')
    b = np.matmul(P, b)  # because pivot changes haven't been applied to b yet
    print('Ly = b')
    y = fsub(L, b)  # solve the system Ly = b

    print('y = \n', y, '\n')
    x = bsub(U, y)  # solve the system Ux = y
    print('Ux = y')
    print('x = \n', x, '\n')
    return x


def cholesky(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    for i in range(n):
        for k in range(i + 1):
            s = 0
            for j in range(k):
                s += L[i, j] * L[k, j]

            if i == k:  # Diagonal elements
                L[i, k] = math.sqrt(A[i, i] - s)
            else:
                L[i][k] = (1.0 / L[k, k] * (A[i, k] - s))
    print('L =')
    print(L, '\n')
    return L


"""
# An example with other arrays to test it's not only working for the given example
A = np.array([[0.02, 0.01, 0., 0.],
             [1., 2., 1., 0.],
             [0., 1., 2., 1.],
             [0., 0., 100., 200.]])
b = np.array([[0.02],
              [1.],
              [4.],
              [800.]])
"""
# A is the array of the example in the slides
A = np.array([[10., -7., 0.],
             [-3., 2., 6.],
             [5., -1., 5.]])
b = np.array([[7.],
              [4.],
              [6.]])

gauss(A, b)


# here are two matrices that are symmetric and positive-definite
"""
Arr = np.array([[6., 3., 4., 8.],
                [3., 6., 5., 1.],
                [4., 5., 10., 7.],
                [8., 1., 7., 25.]])
"""

Arr = np.array([[5., 1., 1.],
                [1., 5., 1.],
                [1., 1., 2.]])


L = cholesky(Arr)
# to check we can just implement a function that transposes L and then we check if A == L*L^T


def transposer(L):
    n = L.shape[0]
    m = L.shape[1]
    LT = np.zeros((n, m))
    for i in range(n):
        for j in range(n):
            LT[i, j] = L[j, i]
    print('LT =')
    print(LT, '\n')
    return LT


def matrix_multiplication(A, B):
    n = A.shape[0]
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(n):
                s += A[i, k] * B[k, j]
            C[i, j] = s
    return C


LT = transposer(L)
print('L*LT = ')
print(matrix_multiplication(L, LT), '\n')  # this of course is equal with the original matrix


# Gauss-Seidel

def copy(A1):
    n = A1.shape[0]
    A2 = np.zeros((n, 1))
    for i in range(n):
        A2[i] = A1[i]
    return A2


def inf_norm(Ar):
    n = Ar.shape[0]
    maxs = abs(Ar[0])
    for i in range(1, n):
        if abs(Ar[i]) > maxs:
            maxs = abs(Ar[i])
    return maxs


def Gauss_Seidel_10(A, b):
    n = A.shape[0]
    x = np.zeros((n, 1))
    xold = copy(x)
    k = 1
    for i in range(n):
        d = b[i, 0]
        for j in range(n):
            if (i - 1) == j or i == (j - 1):
                d -= A[i, j] * x[j, 0]
        x[i, 0] = d / A[i, i]
    diff = np.subtract(x, xold)
    while inf_norm(diff) > 0.00005:
        xold = copy(x)
        k += 1
        for i in range(n):
            d = b[i, 0]
            for j in range(n):
                if (i - 1) == j or i == (j - 1):
                    d -= A[i, j] * x[j, 0]
            x[i, 0] = d / A[i, i]
        diff = np.subtract(x, xold)
    print('k = ', k)
    print('x =')
    print(x)
    print("\n")
    return x


def Gauss_Seidel_10000(A, b):
    n = A.shape[0]
    x = np.zeros((n, 1))
    xold = copy(x)
    k = 1

    d = b[0, 0]
    d -= A[0, 1] * x[1, 0]
    x[0, 0] = d / A[0, 0]
    for i in range(1, 9999):
        d = b[i, 0]
        d -= A[i, i + 1] * x[i + 1, 0]
        d -= A[i, i - 1] * x[i - 1, 0]
        x[i, 0] = d / A[i, i]
    d = b[9999, 0]
    d -= A[9999, 9998] * x[9998, 0]
    x[9999, 0] = d / A[9999, 9999]
    diff = np.subtract(x, xold)
    while inf_norm(diff) > 0.00005:
        xold = copy(x)
        k += 1

        d = b[0, 0]
        d -= A[0, 1] * x[1, 0]
        x[0, 0] = d / A[0, 0]
        for i in range(1, 9999):
            d = b[i, 0]
            d -= A[i, i + 1] * x[i + 1, 0]
            d -= A[i, i - 1] * x[i - 1, 0]
            x[i, 0] = d / A[i, i]
        d = b[9999, 0]
        d -= A[9999, 9998] * x[9998, 0]
        x[9999, 0] = d / A[9999, 9999]
        diff = np.subtract(x, xold)
    print('k = ', k)
    print('x =')
    print(x)
    return x


# create the two arrays for n = 10
A1 = np.zeros((10, 10))

for i in range(10):
    for j in range(10):
        if i == j:
            A1[i, j] = 5
        elif (i - 1) == j or i == (j - 1):
            A1[i, j] = -2


b1 = np.ones((10, 1))
b1[0, 0] = b1[9, 0] = 3

""""""

A2 = np.zeros((10000, 10000))

for i in range(1, 9999):
    A2[i, i] = 5
    A2[i, i + 1] = -2
    A2[i, i - 1] = -2
A2[0, 0] = 5
A2[9999, 9999] = 5
A2[0, 1] = -2
A2[9999, 9998] = -2


b2 = np.ones((10000, 1))
b2[0, 0] = b2[9999, 0] = 3


Gauss_Seidel_10(A1, b1)
Gauss_Seidel_10000(A2, b2)
