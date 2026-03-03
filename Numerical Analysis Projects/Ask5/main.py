import numpy as np
import matplotlib.pyplot as plt
import math


A = np.array([0, math.pi/9, 2*math.pi/9, math.pi/3, 4*math.pi/9, 5*math.pi/9, 2*math.pi/3, 7*math.pi/9, 8*math.pi/9, math.pi])
sinA = np.zeros(10)

for h in range(10):
    sinA[h] = round(math.sin(A[h]), 5)


def p(a):
    Lx = np.ones((1, 10))
    p = 0
    for i in range(10):
        for j in range(10):
            if j != i:
                Lx[0, i] *= (a - A[j]) / (A[i] - A[j])
        p += sinA[i] * Lx[0, i]
    return p


def f(b):
    flag = 0
    if b < 0:  # for the oddness of sinx
        b = -b
        flag = 1
    b %= 2*math.pi  # for the periodicity of sinx
    if math.pi <= b <= 2*math.pi:  # for the [π, 2π]
        b -= math.pi
        if flag == 0:
            return -p(b)
        return p(b)
    if flag == 0:
        return p(b)
    return -p(b)


def inf_normf():
    Xerr = np.arange(-math.pi, math.pi, math.pi / 100)
    # Xerr[i] = -math.pi + (math.pi / 100) * i
    Err = np.zeros(200)
    Err[0] = abs(np.sin(Xerr[0]) - f(Xerr[0]))
    max = Err[0]
    maxi = 0
    for i in range(1, 200):
        Err[i] = abs(np.sin(Xerr[i]) - f(Xerr[i]))
        if Err[i] > max:
            max = Err[i]
            maxi = i
    print("Lagrange Interpolation: ")
    print("max error = ", max)
    print("maxi = ", maxi)
    return Err


Arr = inf_normf()
Xaxis = np.zeros(200)
for i in range(200):
    Xaxis[i] = -math.pi + (math.pi / 100) * i
plt.plot(Xaxis, Arr, 'r')
plt.ylabel('Error')
plt.xlabel('Xerr cells in order')
plt.title('Graph of the error of Lagrange Interpolation')
plt.show()

# Here we can check that the graph of f indeed is like sin's
x = np.arange(-math.pi, math.pi, 0.01)
plt.plot(x, list(map(f, x)))
plt.ylim(-1, 1)
plt.ylabel('f(x)')
plt.xlabel('x')
plt.title('Graph of the Lagrange interpolation of sin')
# plt.plot(x, list(map(math.sin, x)))
plt.show()


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
    # print('PA = LU')
    # print('P = \n', P)

    # print('L = \n', L)

    # print('U = \n', U, '\n')
    b = np.matmul(P, b)  # because pivot changes haven't been applied to b yet
    # print('Ly = b')
    y = fsub(L, b)  # solve the system Ly = b

    # print('y = \n', y, '\n')
    x = bsub(U, y)  # solve the system Ux = y
    # print('Ux = y')
    # print('x = \n', x, '\n')
    return x


def P(rad):
    A1 = np.zeros((10, 10))
    pi = math.pi
    for i in range(1, 9):
        A1[i, i] = (4 * pi) / 9
        A1[i, i + 1] = pi / 9
        A1[i, i - 1] = pi / 9
    A1[9, 9] = A1[0, 0] = 1
    A1[9, 8] = A1[0, 1] = 0

    r = np.zeros((10, 1))
    for i in range(1, 9):
        r[i, 0] = (9/pi) * ((sinA[i+1] - sinA[i]) - (sinA[i] - sinA[i-1]))
    r[0, 0] = 0
    r[9, 0] = 0
    r *= 3

    c = gauss(A1, r)
    o = 1
    for i in range(1, 10):
        if rad <= A[i]:
            o = i
            break

    b = (9/pi) * (sinA[o] - sinA[o-1]) - (pi/27) * (2*c[o-1, 0] + c[o, 0])

    d = (3/pi) * (c[o, 0] - c[o-1, 0])
    s = rad - A[o-1]
    return sinA[o-1] + s*b + (s ** 2)*c[o-1, 0] + (s ** 3)*d


def fspl(b):
    flag = 0
    if b <= 0:
        b = -b
        flag = 1
    b %= 2*math.pi
    if math.pi <= b <= 2*math.pi:
        b -= math.pi
        if flag == 0:
            return -P(b)
        return P(b)
    if flag == 0:
        return P(b)
    return -P(b)


def inf_spl():
    Xerr = np.arange(-math.pi, math.pi, math.pi / 100)
    # Xerr[i] = -math.pi + (math.pi / 100) * i
    Err1 = np.zeros(200)
    Err1[0] = abs(np.sin(Xerr[0]) - fspl(Xerr[0]))
    max1 = Err1[0]
    maxi1 = 0
    for i in range(1, 200):
        Err1[i] = abs(np.sin(Xerr[i]) - fspl(Xerr[i]))
        if Err1[i] > max1:
            max1 = Err1[i]
            maxi1 = i
    print("Cubic Splines: ")
    print("max error = ", max1)
    print("maxi = ", maxi1)
    return Err1


Arr1 = inf_spl()
Xaxis1 = np.zeros(200)
for i in range(200):
    Xaxis1[i] = -math.pi + (math.pi / 100) * i
plt.plot(Xaxis1, Arr1, 'r')
plt.ylabel('Error')
plt.xlabel('Xerr cells in order')
plt.title('Graph of the error of cubic splines')
plt.show()


x1 = np.arange(-math.pi, math.pi, 0.01)
# plt.plot(x1, list(map(math.sin, x1)))
x = np.arange(A[0], A[2], 0.01)
plt.plot(x1, list(map(fspl, x1)))
plt.ylim(-1, 1.5)
plt.ylabel('P(x)')
plt.xlabel('x')
plt.title('Graph of the Cubic Splines interpolation of sin')
plt.show()


def transposer(L):
    n = L.shape[0]
    m = L.shape[1]
    LT = np.zeros((m, n))
    temp = np.zeros(m)
    for i in range(n):
        for j in range(m):
            temp[j] = L[i, j]
        for j in range(m):
            LT[j, i] = temp[j]
    return LT


def matrix_multiplication(A, B):
    n = A.shape[0]
    m = A.shape[1]
    p = B.shape[1]
    C = np.zeros((n, p))
    for i in range(n):
        for j in range(p):
            s = 0
            for k in range(m):
                s += A[i, k] * B[k, j]
            C[i, j] = s
    return C


def least_squares(x):
    X = np.zeros((10, 4))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = sinA[i]
    for i in range(10):
        X[i] = [1, A[i], A[i] ** 2, A[i] ** 3]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2) + β[3, 0] * (x ** 3)


def fls(b):
    flag = 0
    if b <= 0:
        b = -b
        flag = 1
    b %= 2*math.pi
    if math.pi <= b <= 2*math.pi:
        b -= math.pi
        if flag == 0:
            return -least_squares(b)
        return least_squares(b)
    if flag == 0:
        return least_squares(b)
    return -least_squares(b)


def inf_normls():
    Xerr = np.arange(-math.pi, math.pi, math.pi / 100)
    # Xerr[i] = -math.pi + (math.pi / 100) * i
    Err1 = np.zeros(200)
    Err1[0] = abs(np.sin(Xerr[0]) - fls(Xerr[0]))
    max1 = Err1[0]
    maxi1 = 0
    for i in range(1, 200):
        Err1[i] = abs(np.sin(Xerr[i]) - fls(Xerr[i]))
        if Err1[i] > max1:
            max1 = Err1[i]
            maxi1 = i
    print("Least Squares: ")
    print("max error = ", max1)
    print("maxi = ", maxi1)
    return Err1


Arr1 = inf_normls()
Xaxis1 = np.zeros(200)
for i in range(200):
    Xaxis1[i] = -math.pi + (math.pi / 100) * i
plt.plot(Xaxis1, Arr1, 'r')
plt.ylabel('Error')
plt.xlabel('Xerr cells in order')
plt.title('Graph of the error of least squares')
plt.show()

# Here we can check that the graph of the least squares cubic is indeed like sin's
x1 = np.arange(-math.pi, math.pi, 0.01)
plt.plot(x1, list(map(fls, x1)))
plt.ylim(-1, 1)
plt.ylabel('least_squares(x)')
plt.xlabel('x')
plt.title('Graph of the least squares cubic of sin')
# plt.plot(x, list(map(math.sin, x)))
plt.show()
