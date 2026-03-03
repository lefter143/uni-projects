import numpy as np
import math


# π/2 as the interval so for 11 equally spaced points we get (π/2)/10=π/20 as the size of our intervals(we'll have 10
# of them)
A = np.array([0, math.pi/20, math.pi/10, 3*math.pi/20, math.pi/5, math.pi/4, 3*math.pi/10, 7*math.pi/20, 2*math.pi/5,
              9*math.pi/20, math.pi/2])

sinA = np.zeros(11)
for h in range(11):
    sinA[h] = math.sin(A[h])


def trapezoid_rule(x, f):
    N = x.shape[0]
    S = 0
    for i in range(1, N-1):
        S += f[i]
    return ((x[N-1] - x[0]) / (2*(N-1))) * (f[0] + f[N-1] + 2*S)


def Simpson_rule(x, f):
    N = x.shape[0]
    S1 = S2 = 0
    for i in range(1, N//2):
        S1 += f[2*i]
    for i in range(1, N//2+1):
        S2 += f[2*i-1]
    return ((x[N-1] - x[0]) / (3*(N-1))) * (f[0] + f[N-1] + 2*S1 + 4*S2)


I1 = trapezoid_rule(A, sinA)
print("The result with trapezoid rule: ", I1)
print("The arithmetic error with trapezoid rule: ", abs(1-I1))  # since if we compute it in the classic calculus way
# we get 1
N = A.shape[0] - 1
M = 1  # since we know that the 2nd derivative of sin(x) is -sin(x) then abs(-sin(x))=sin(x) and so in the
# interval [0, π/2] its maximum value is 1
e = (((math.pi/2) ** 3) * M) / (12*N)
print("The theoretical maximum error with trapezoid rule is: ", e)

""""""""""""

I2 = Simpson_rule(A, sinA)
print("\nThe result with Simpson rule: ", I2)
print("The arithmetic error with Simpson rule: ", abs(1-I2))  # since if we compute it in the classic calculus way
# we get 1
N = A.shape[0] - 1
M = 1  # the 4th derivative of sin(x) is sin(x) then abs(sin(x))=sin(x) and so in the
# interval [0, π/2] its maximum value is 1 (the same as before)
e = (((math.pi/2) ** 5) * M) / ((180*N) ** 4)
print("The theoretical maximum error with Simpson rule is: ", e)
