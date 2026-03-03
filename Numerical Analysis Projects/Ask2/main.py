import math
import random


def f(x):
    return 54 * math.pow(x, 6) + 45 * math.pow(x, 5) - 102 * math.pow(x, 4) - 69 * math.pow(x, 3) + 35 * math.pow(x, 2) + 16 * x - 4


def df(x):
    return 324 * math.pow(x, 5) + 225 * math.pow(x, 4) - 408 * math.pow(x, 3) - 207 * math.pow(x, 2) + 70 * x + 16


def d2f(x):
    return 1620 * math.pow(x, 4) + 900 * math.pow(x, 3) - 1224 * math.pow(x, 2) - 414 * x + 70


# Newton-Raphson


def newton_raphson(x0):
    print('Newton-Raphson for f:')
    k = 1
    xprev = x0
    x = xprev - f(xprev)/df(xprev)
    while abs(x - xprev) > abs(x) * 0.00001:
        xprev = x
        x = xprev - f(xprev)/df(xprev)
        k = k + 1
    print('k = ', k)
    return x


def modified_newton_raphson(x0):
    print('Modified Newton-Raphson for f:')
    k = 1
    xprev = x0
    x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
    while abs(x - xprev) > abs(x) * 0.00001:
        xprev = x
        x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
        k = k + 1
    print('k = ', k)
    return x


def newton_raphson_more_precision_for_multiple_roots(x0):
    print('Newton-Raphson with more precision for the multiple roots of f:')
    k = 1
    xprev = x0
    x = xprev - f(xprev)/df(xprev)
    while abs(x - xprev) > abs(x) * 0.000008:
        xprev = x
        x = xprev - f(xprev)/df(xprev)
        k = k + 1
    print('k = ', k)
    return x


def altered_modified_newton_raphson(x0):
    print('Altered modified Newton-Raphson for f:')
    k = 1
    xprev = x0
    x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
    while abs(x - xprev) > abs(x) * 0.00001:
        xprev = x
        x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
        k = k + 1
    k -= 1  # in this case it needs one iteration less since we return xprev
    print('k = ', k)
    return xprev


# Bisection
def bisection(a, b):
    print('Bisection for f:')
    n = int((math.log(b - a) - math.log(0.000005)) / math.log(2) + 1)
    print('k = ', n)
    for j in range(n):
        m = (a + b) / 2
        if (f(m) * f(a)) < 0:
            b = m
        elif (f(m) * f(a)) > 0:
            a = m
        else:
            return m
    return m


def bisection_rand(a, b):
    print('Randomized Bisection for f:')
    k = 1
    m = random.uniform(a, b)
    if (f(m) * f(a)) < 0:
        b = m
    elif (f(m) * f(a)) > 0:
        a = m
    else:
        print('k = ', k)
        return m
    while abs(b - a) > abs(b) * 0.00001:  # to make sure that we get the precision we want since we don't know beforehand the error rate because it's random
        k = k + 1
        m = random.uniform(a, b)
        if (f(m) * f(a)) < 0:
            b = m
        elif (f(m) * f(a)) > 0:
            a = m
        else:
            print('k = ', k)
            return m
    print('k = ', k)
    return m


def bisection_rand_test(a, b):
    k = 1
    m = random.uniform(a, b)
    if (f(m) * f(a)) < 0:
        b = m
    elif (f(m) * f(a)) > 0:
        a = m
    else:
        return k
    while abs(b - a) > 0.000001:  # to make sure that we get the precision we want since we don't know beforehand the error rate
        k = k + 1
        m = random.uniform(a, b)
        if (f(m) * f(a)) < 0:
            b = m
        elif (f(m) * f(a)) > 0:
            a = m
        else:
            return k
    return k


# Secant
def secant(x0, x1):
    print('Secant for f:')
    k = 1
    a = x0
    b = x1
    c = b - (f(b) * (b-a)) / (f(b) - f(a))
    while abs(c - b) > abs(c) * 0.00001:
        a = b
        b = c
        c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
        k = k + 1
    print('k = ', k)
    return c


def modified_secant(x0, x1, x2):
    print('Modified secant for f:')
    k = 1
    a = x0
    b = x1
    c = x2
    q = f(a) / f(b)
    r = f(c) / f(b)
    s = f(c) / f(a)
    d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
    while abs(d - c) > abs(d) * 0.0001:
        a = b
        b = c
        c = d
        q = f(a) / f(b)
        r = f(c) / f(b)
        s = f(c) / f(a)
        d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
        k = k + 1
    k -= 1  # because we actually do one more iteration than needed but we do it for the root x4 = 0.5 to be displayed properly
    print('k = ', k)
    return d


def secant_more_precision_for_multiple_roots(x0, x1):
    print('Secant with more precision for the multiple roots of f:')
    k = 1
    a = x0
    b = x1
    c = b - (f(b) * (b-a)) / (f(b) - f(a))
    while abs(c - b) > abs(c) * 0.000005:
        a = b
        b = c
        c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
        k = k + 1
    print('k = ', k)
    return c


def modified_secant_more_precision_for_multiple_roots(x0, x1, x2):
    print('Modified secant with more precision for the multiple roots of f:')
    k = 1
    a = x0
    b = x1
    c = x2
    q = f(a) / f(b)
    r = f(c) / f(b)
    s = f(c) / f(a)
    d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
    while abs(d - c) > abs(d) * 0.000005:
        a = b
        b = c
        c = d
        q = f(a) / f(b)
        r = f(c) / f(b)
        s = f(c) / f(a)
        d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
        k = k + 1
    print('k = ', k)
    return d


# Newton-Raphson
print('x1 = ', newton_raphson(-1.5), '\n')
print('x1 = ', modified_newton_raphson(-1.5), '\n')

print('f(x2)*d2f(x2) = ', abs(f(-0.6666666666666)*d2f(-0.666666666666)))  # since it's equal to 0 we have a multiple root so it will be
# slower and the rate in which the error gets smaller will decrease with each iteration more than the point of  our
# precision which is 0.00001 before finding the root  so we need to check for a smaller error than what should be
# standard for our wanted precision (I decreased it to  0.000005) since it goes too slow. So we call the
# newton_raphson_more_precision_for_multiple_roots function instead.  Note that this isn't needed for the modified
# newton-raphson since it seems to be faster in general than the normal  newton-raphson, probably because it uses more
# information, namely the d2f, and similarly to what happens to the  secant method that doesn't use the derivative and
# turns out to be slower than newton-raphson(1.62 convergence instad  of quadratic), the normal newton-raphson is
# slower than the modified one

print('df(x2) = ', df(-0.666666666666666666666666666666666666))  # zero derivative
print('x2 = ', newton_raphson_more_precision_for_multiple_roots(-0.5), '\n')
print('x2 = ', modified_newton_raphson(-0.5), '\n')

print('x3 = ', newton_raphson(0), '\n')
print('x3 = ', modified_newton_raphson(0), '\n')

# since the modified newton-raphson reaches the root too fast it has an error rate that's too high even after it
# reaches the root so it does one more iteration and the error plummets again since it's already at the wanted precision
# but then it has more precision than we want (it's technically closer to the actual root 0.5 than the previous
# iteration) but it is from the back, namely 49999999999999994 so in order to avoid this we can call a different
# modified newton-raphson that just returns xprev instead of x
print('x4 = ', newton_raphson(0.75), '\n')
print('x4 = ', altered_modified_newton_raphson(0.75), '\n')

print('x5 = ', newton_raphson(2), '\n')
print('x5 = ', modified_newton_raphson(2), '\n')


# Bisection
print('x1 = ', bisection(-1.5, -1), '\n')
print('x1 = ', bisection_rand(-1.5, -1), '\n')

print('x2 = ', bisection(-0.75, -0.5), '\n')
print('x2 = ', bisection_rand(-0.75, -0.5), '\n')

print('x3 = ', bisection(0, 0.25), '\n')
print('x3 = ', bisection_rand(0, 0.25), '\n')

print('x4 = ', bisection(0.3, 0.6), '\n')
print('x4 = ', bisection_rand(0.3, 0.6), '\n')

print('x5 = ', bisection(1, 1.5), '\n')
print('x5 = ', bisection_rand(1, 1.5), '\n')


# Try 20 times the random bisection to see if it converges to a number of iterations (doesn't matter which root we
# choose for this)
avg = 0
for i in range(20):
    p = bisection_rand_test(1, 1.5)
    avg = avg + p
avg = avg / 20
print('average = ', avg, '\n')
# we can see that the average number of iterations seems to be around 28 (this changes according to our wanted level of
# precision) so this confirms that it converges

# Secant
print('x1 = ', secant(-1.5, -1.2), '\n')
print('x1 = ', modified_secant(-1.5, -1.3, -1.2), '\n')

print('x2 = ', secant_more_precision_for_multiple_roots(-0.8, -0.6), '\n')
print('x2 = ', modified_secant_more_precision_for_multiple_roots(-0.8, -0.7, -0.6), '\n')

print('x3 = ', secant(0.1, 0.25), '\n')
print('x3 = ', modified_secant(0.1, 0.2, 0.25), '\n')

print('x4 = ', secant(0.3, 0.7), '\n')
print('x4 = ', modified_secant(0.3, 0.6, 0.7), '\n')

print('x5 = ', secant(1, 1.3), '\n')
print('x5 = ', modified_secant(1, 1.1, 1.3), '\n')

# Experimental speed comparison between the different methods


def newton_raphson_experiment(x0):
    print('Experiment for Newton-Raphson for f:')
    k = 1
    xprev = x0
    x = xprev - f(xprev)/df(xprev)
    while abs(x - xprev) > abs(x) * 0.00001:
        xprev = x
        x = xprev - f(xprev)/df(xprev)
        k = k + 1
        print('abs(f(x)/f(xprev)) = ', abs(f(x) / f(xprev)))
    print('k = ', k)
    return x


def modified_newton_raphson_experiment(x0):
    print('Experiment for modified Newton-Raphson for f:')
    k = 1
    xprev = x0
    x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
    while abs(x - xprev) > abs(x) * 0.00001:
        xprev = x
        x = xprev - 1 / ((df(xprev) / f(xprev)) - (d2f(xprev) / (2 * df(xprev))))
        k = k + 1
        print('abs(f(x)/f(xprev)) = ', abs(f(x) / f(xprev)))
    print('k = ', k)
    return x


def bisection_experiment(a, b):
    print('Experiment for bisection for f:')
    n = int((math.log(b - a) - math.log(0.000005)) / math.log(2) + 1)
    print('k = ', n)
    for j in range(n):
        m = (a + b) / 2
        if (f(m) * f(a)) < 0:
            b = m
        elif (f(m) * f(a)) > 0:
            a = m
        else:
            return m
    return m


def secant_experiment(x0, x1):
    print('Experiment for secant for f:')
    k = 1
    a = x0
    b = x1
    c = b - (f(b) * (b-a)) / (f(b) - f(a))
    while abs(c - b) > abs(c) * 0.00001:
        a = b
        b = c
        c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
        k = k + 1
        print('abs(f(c)/f(b)) = ', abs(f(c) / f(b)))
    print('k = ', k)
    return c


def modified_secant_experiment(x0, x1, x2):
    print('Experiment for modified secant for f:')
    k = 1
    a = x0
    b = x1
    c = x2
    q = f(a) / f(b)
    r = f(c) / f(b)
    s = f(c) / f(a)
    d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
    while abs(d - c) > abs(d) * 0.00001:
        a = b
        b = c
        c = d
        q = f(a) / f(b)
        r = f(c) / f(b)
        s = f(c) / f(a)
        d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ((q - 1) * (r - 1) * (s - 1))
        k = k + 1
        print('abs(f(d)/f(c)) = ', abs(f(d) / f(c)))
    print('k = ', k)
    return d


print('x5 = ', newton_raphson_experiment(2), '\n')
# we can see that this (almost) follows (at the last iterations of course) the quadratic rate: (10^-3, 10^-6)
print('x5 = ', modified_newton_raphson_experiment(2), '\n')
# we can see that this (almost) follows (at the last iterations of course) the cubic rate: (10^-4, 10^-12)


print('x5 = ', bisection_experiment(1, 1.5), '\n')
# we know this is linear (in fact we know the exact number of iterations k = 17)

avg1 = 0
for i in range(1000):
    p = bisection_rand_test(-1.5, -1)
    avg1 += p
avg1 = avg1 / 1000
print('average of bisection_rand x1 = ', avg1, '\n')

avg2 = 0
for i in range(1000):
    p = bisection_rand_test(-0.75, -0.5)
    avg2 += p
avg2 = avg2 / 1000
print('average of bisection_rand x2 = ', avg2, '\n')

avg3 = 0
for i in range(1000):
    p = bisection_rand_test(0, 0.25)
    avg3 += p
avg3 = avg3 / 1000
print('average of bisection_rand x3 = ', avg3, '\n')

avg4 = 0
for i in range(1000):
    p = bisection_rand_test(0.3, 0.6)
    avg4 += p
avg4 = avg4 / 1000
print('average of bisection_rand x4 = ', avg4, '\n')

avg5 = 0
for i in range(1000):
    p = bisection_rand_test(1, 1.5)
    avg5 += p
avg5 = avg5 / 1000
print('average of bisection_rand x5 = ', avg5, '\n')
# this is slower as the average number of iterations that we get for this root is 28, so we can conclude that it's on
# average slower than linear


print('x5 = ', secant_experiment(1, 1.3), '\n')
# we can see that this (almost) follows (at the last iterations of course) the quadratic rate: (10^-3, 10^-6) but it's
# a little less than that because the secant method has convergence of p ≈ 1.6
print('x5 = ', modified_secant_experiment(1, 1.1, 1.3), '\n')
# we can see that this (almost) follows (at the last iterations of course) the cubic rate: (10^-4, 10^-8) but it's
# a little less than that because the modified secant method, similarly to the normal secant, doesn't use f', f'' so
# it's a little slower than cubic convergence but more than quadratic
