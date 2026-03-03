import numpy as np
import math


A = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
              [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0]])

q = 0.15
n = A.shape[0]
# create G
G = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += A[j, k]
        G[i, j] = q / n + (A[j, i] * (1 - q)) / nj

# proof that G is stochastic
check = 0
for i in range(n):
    s = 0
    for j in range(n):
        s += G[j, i]
    if s != 1:
        check = 1
if check == 0:
    print("The matrix G is stochastic \n")
# since the sums of all the columns of G are 1 then G is left stochastic


def copy(A1):
    n = A1.shape[0]
    A2 = np.zeros((n, 1))
    for i in range(n):
        A2[i, 0] = A1[i, 0]
    return A2


def inf_norm(Ar):
    n = Ar.shape[0]
    maxs = abs(Ar[0, 0])
    for i in range(1, n):
        if abs(Ar[i, 0]) > maxs:
            maxs = abs(Ar[i, 0])
    return maxs


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


def straightSelectionSort(arr):
    n = arr.shape[0]
    rankings = np.zeros((n, 1))
    for i in range(n):
        rankings[i, 0] = i + 1
    C = np.zeros((n, 1))
    C = copy(arr)

    for i in range(n):
        kmax = i
        xmax = C[i, 0]
        for j in range(i + 1, n):
            if xmax > C[j, 0]:
                kmax = j
                xmax = C[j, 0]
        C[kmax, 0] = C[i, 0]
        C[i, 0] = xmax
        temp = rankings[i, 0]
        rankings[i, 0] = rankings[kmax, 0]
        rankings[kmax, 0] = temp
    return rankings


def find_p(A):
    n = A.shape[0]
    k = 1
    p = np.ones((n, 1)) / n
    pold = copy(p)
    p = matrix_multiplication(A, p)
    diff = np.subtract(p, pold)
    while inf_norm(diff) > 5 * 1e-10:
        k += 1
        pold = copy(p)
        p = matrix_multiplication(A, p)
        diff = np.subtract(p, pold)
    print('k = ', k)
    return p


p = find_p(G)
s = 0
for i in range(15):
    s += p[i, 0]  # we can't have infinite precision after all
print('s = ', s)
print('p = ')
print(p)  # so we've proven that this p is the same one as the given one
ranks = straightSelectionSort(p)
print('ranks = ')
for i in range(n-1, 0, -1):
    print(ranks[i, 0])
print("\n")

# let's say that I want to increase page 9's ranking of importance from 4th (15 and 13, 10 and 11 have the
# same rank since p[15,0] == p[13, 0] and p[10, 0] == p[11, 0]). So, if I add 4 connections from high ranked pages to it
# and cut its connection with page 10 (the highest ranked page that it's connected to) I can raise its rank

# the new A then becomes:
Anew = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                 [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                 [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
                 [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0]])
# A[8, 9] is removed (therefore connection 9→10severed) , A[13, 8] added(therefore connection 14→9 added),
# A[9, 8] added(therefore connection 10→9 added), A[10, 8] added(therefore connection 11→9 added),
# A[14, 8] added(therefore connection 15→9 added)

# creation of new G
q = 0.15
n = Anew.shape[0]
Gnew = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += Anew[j, k]
        Gnew[i, j] = q / n + (Anew[j, i] * (1 - q)) / nj

pnew = find_p(Gnew)
print('pnew = ')
print(pnew)
newranks = straightSelectionSort(pnew)
print('newranks = ')
for i in range(n-1, 0, -1):
    print(newranks[i, 0])
print("\n")

# we can observe that page 9 now rises to 1st rank, and all the pages that it connects to have risen in rank (pages 5,
# 6, 10 are at the top now)

# Experiment to see what happens if we change the transition probability q to (α) q = 0.02 and (β) q = 0.6 in our
# changed graph

# creation of new G
q = 0.02  # (α)
n = Anew.shape[0]
Gnewa = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += Anew[j, k]
        Gnewa[i, j] = q / n + (Anew[j, i] * (1 - q)) / nj

pnewa = find_p(Gnewa)
print('pnewa = ')
print(pnewa)
newranksa = straightSelectionSort(pnewa)
print('newranksa = ')
for i in range(n-1, 0, -1):
    print(newranksa[i, 0])
print("\n")

# creation of new G
q = 0.6  # (β)
n = Anew.shape[0]
Gnewb = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += Anew[j, k]
        Gnewb[i, j] = q / n + (Anew[j, i] * (1 - q)) / nj

pnewb = find_p(Gnewb)
print('pnewb = ')
print(pnewb)
newranksb = straightSelectionSort(pnewb)
print('newranksb = ')
for i in range(n-1, 0, -1):
    print(newranksb[i, 0])


ni = np.zeros(15)
for i in range(15):
    for j in range(15):
        ni[i] += Anew[j, i]
    print(i + 1, ': ', ni[i])

# we can see that the effect that the transition probability q has, as it increases, is that it reduces the influence of
# pages which themselves have higher traffic coming to them. In other words, if q  increases, the number of connections
# coming to a page, becomes more important than the rank of the pages with which  it's connected to. This is pretty much
# obvious from the way that we compute the matrix G(namely q / n is the term that increases and A(j,i) * (1 - q) / nj
# is the one that decreases). Another effect is that it makes the ranks p of all the pages closer to all being equal
# (which again can be seen from the formula for the creation of G). The opposite is true if q decreases, that the rank
# of the pages which connect to a page is more important than the number of pages connecting to it


# changed a with 3 at A[7, 10](therefore connection 8→11 is made stronger), changed a with 3 at A[11, 10](therefore
# connection 12→11 is made stronger),
A11 = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 3, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0]])
# creation of new G
q = 0.15
n = A11.shape[0]
G11 = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += A11[j, k]
        G11[i, j] = q / n + (A11[j, i] * (1 - q)) / nj

p11 = find_p(G11)
print('p11 = ')
print(p11)
ranks11 = straightSelectionSort(p11)
print('ranks11 = ')
for i in range(n-1, 0, -1):
    print(ranks11[i, 0])
print("\n")
# we can see that this strategy works as 11 is now 2nd in rank as compared to 10 that stays in the same place as before


# A with deleted the 10th row and column
Ano10 = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                  [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                  [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                  [0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0]])


def straightSelectionSortno10(arr):
    n = arr.shape[0]
    rankings = np.zeros((n, 1))
    k = 0
    for i in range(n):
        if i != 9:
            k += 1
        else:
            k += 2
        rankings[i, 0] = k
    C = np.zeros((n, 1))
    C = copy(arr)

    for i in range(n):
        kmax = i
        xmax = C[i, 0]
        for j in range(i + 1, n):
            if xmax > C[j, 0]:
                kmax = j
                xmax = C[j, 0]
        C[kmax, 0] = C[i, 0]
        C[i, 0] = xmax
        temp = rankings[i, 0]
        rankings[i, 0] = rankings[kmax, 0]
        rankings[kmax, 0] = temp
    return rankings


# creation of new G
q = 0.15
n = Ano10.shape[0]
Gno10 = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        nj = 0
        for k in range(n):
            nj += Ano10[j, k]
        Gno10[i, j] = q / n + (Ano10[j, i] * (1 - q)) / nj

pno10 = find_p(Gno10)
print('pno10 = ')
print(pno10)
ranksno10 = straightSelectionSortno10(pno10)
print('ranksno10 = ')
for i in range(n - 1, 0, -1):
    print(ranksno10[i, 0])
print("\n")
# the rank of page 13 is lowered as its rank was directly affected by its connection with 10. So page 11 which
# isn't affected by the change as it has no connection from page 10 rises to 2nd rank, 12 rises rank as it isn't
# affected by 10 too. So generally the pages that aren't affected by the changes, are going to rise ranks and the ones
# that lost the connection coming from 10 are going to lower ranks
