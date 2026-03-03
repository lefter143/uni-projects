import numpy as np
import matplotlib.pyplot as plt


dates = np.array([2, 3, 4, 5, 6, 9, 10, 11, 12, 13])
stocksΔΕΗ = np.array([9.255, 9, 9.19, 9.09, 8.805, 8.65, 9.175, 9.47, 9.195, 9.24])
stocksΟΤΕ = np.array([13.75, 13.53, 13.49, 13.61, 13.27, 12.89, 13.42, 13.25, 13.24, 13.28])


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


def least_squares_quadraticΟΤΕ(x):
    X = np.zeros((10, 3))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΟΤΕ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2)


def least_squares_quadraticΔΕΗ(x):
    X = np.zeros((10, 3))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΔΕΗ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2)


def least_squares_cubicΟΤΕ(x):
    X = np.zeros((10, 4))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΟΤΕ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2, dates[i] ** 3]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2) + β[3, 0] * (x ** 3)


def least_squares_cubicΔΕΗ(x):
    X = np.zeros((10, 4))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΔΕΗ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2, dates[i] ** 3]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2) + β[3, 0] * (x ** 3)


def least_squares_quarticΟΤΕ(x):
    X = np.zeros((10, 5))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΟΤΕ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2, dates[i] ** 3, dates[i] ** 4]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2) + β[3, 0] * (x ** 3) + β[4, 0] * (x ** 4)


def least_squares_quarticΔΕΗ(x):
    X = np.zeros((10, 5))
    Y = np.zeros((10, 1))
    for i in range(10):
        Y[i, 0] = stocksΔΕΗ[i]
    for i in range(10):
        X[i] = [1, dates[i], dates[i] ** 2, dates[i] ** 3, dates[i] ** 4]
    XT = transposer(X)
    X = matrix_multiplication(XT, X)
    Y = matrix_multiplication(XT, Y)
    β = gauss(X, Y)
    return β[0, 0] + β[1, 0] * x + β[2, 0] * (x ** 2) + β[3, 0] * (x ** 3) + β[4, 0] * (x ** 4)


real_stocksΔΕΗ = np.array([9.22, 9.295, 9.34, 9.19, 9.06])  # here are the real stock closes for ΔΕΗ after (and
# including) the 16th
real_stocksΟΤΕ = np.array([13.63, 13.45, 13.55, 13.53, 13.42])  # here are the real stock closes for ΟΤΕ after (and
# including) the 16th
dates1 = np.array([2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20])  # here are the 10 days before the 16th(the
# day closest to my birthday) and the 5 following days


plt.plot(dates1, list(map(least_squares_quadraticΟΤΕ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΟΤΕ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_quadraticΟΤΕ(dates1[i]), 'ro')  # the points that we want to predict are plotted
    # with red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΟΤΕ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with quadratic polynomial of ΟΤΕ\'s stock closings')
plt.show()


plt.plot(dates1, list(map(least_squares_quadraticΔΕΗ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΔΕΗ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_quadraticΔΕΗ(dates1[i]), 'ro')  # the points that we want to predict are plotted
    # with red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΔΕΗ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with quadratic polynomial of ΔΕΗ\'s stock closings')
plt.show()

print("QUADRATIC: ")
print("ΔΕΗ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_quadraticΔΕΗ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΔΕΗ[i-10], 4))
print("\nΟΤΕ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_quadraticΟΤΕ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΟΤΕ[i - 10], 4))


""""""""""""


plt.plot(dates1, list(map(least_squares_cubicΟΤΕ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΟΤΕ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_cubicΟΤΕ(dates1[i]), 'ro')  # the points that we want to predict are plotted with
    # red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΟΤΕ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with cubic polynomial of ΟΤΕ\'s stock closings')
plt.show()


plt.plot(dates1, list(map(least_squares_cubicΔΕΗ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΔΕΗ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_cubicΔΕΗ(dates1[i]), 'ro')  # the points that we want to predict are plotted with
    # red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΔΕΗ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with cubic polynomial of ΔΕΗ\'s stock closings')
plt.show()

print("\nCUBIC: ")
print("ΔΕΗ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_cubicΔΕΗ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΔΕΗ[i-10], 4))
print("\nΟΤΕ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_cubicΟΤΕ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΟΤΕ[i - 10], 4))


""""""""""""


plt.plot(dates1, list(map(least_squares_quarticΟΤΕ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΟΤΕ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_quarticΟΤΕ(dates1[i]), 'ro')  # the points that we want to predict are plotted
    # with red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΟΤΕ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with quartic polynomial of ΟΤΕ\'s stock closings')
plt.show()


plt.plot(dates1, list(map(least_squares_quarticΔΕΗ, dates1)))
for i in range(10):
    plt.plot(dates1[i], stocksΔΕΗ[i], 'bo')  # the points that we already know are plotted with blue color
for i in range(10, 15):
    plt.plot(dates1[i], least_squares_quarticΔΕΗ(dates1[i]), 'ro')  # the points that we want to predict are plotted
    # with red color
for i in range(10, 15):
    plt.plot(dates1[i], real_stocksΔΕΗ[i-10], 'go')  # the real values of the stock closes after(and
# including) the 16th are plotted in green
plt.ylabel('least_squares(date)')
plt.xlabel('date')
plt.title('Graph of the Least Squares Approximation with quartic polynomial of ΔΕΗ\'s stock closings')
plt.show()

print("\nQUARTIC: ")
print("ΔΕΗ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_quarticΔΕΗ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΔΕΗ[i-10], 4))
print("\nΟΤΕ: ")
for i in range(10, 15):
    print("Prediction for", dates1[i], "/10/2023: ", round(least_squares_quarticΟΤΕ(dates1[i]), 4))
    print("Real value: ", round(real_stocksΟΤΕ[i - 10], 4))
