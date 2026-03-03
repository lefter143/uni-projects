import math


def f(x):
    return (14 * x - 12) * math.exp(x - 2) - 7 * math.pow(x, 3) + 20 * math.pow(x, 2) - 26 * x + 12


def df(x):
    return (14 * x + 2) * math.exp(x - 2) - 21 * math.pow(x, 2) + 40 * x - 26


def d2f(x):
    return (14 * x + 16) * math.exp(x - 2) - 42 * x + 40


def g(x):
    return math.pow(x - 2, 3)


def dg(x):
    return 3 * math.pow(x - 2, 2)


def bisection(a, b):
    print('Bisection for f:')
    print('k = 19')
    # N = (ln(1.5)-ln(0.5*10^(-k))/ln(2) = 18.19 rounded up makes 19 iterations for precision of k=5 decimal places
    for i in range(19):
        m = (a + b) / 2
        if (f(m) * f(a)) < 0:
            b = m
        elif (f(m) * f(a)) > 0:
            a = m
        else:
            return m
    return m


def newton_raphson_g(x0):
    print('Newton-Raphson for g:')
    k = 1
    xprev = x0
    x = xprev - g(xprev)/dg(xprev)
    err = x - xprev
    while abs(x - xprev) > abs(x) * 0.000002:  # just reduce the precision or increase manually the number of iterations since it's a multiple root
        errprev = err
        xprev = x
        x = xprev - g(xprev)/dg(xprev)
        k = k + 1
        err = x - xprev
        print('abs(err / errprev)) = ', abs(err / errprev))  # it's nearly constant so we have linear convergence
    print('k = ', k)
    return x


def newton_raphson_x1(x0):
    print('Newton-Raphson for f(x1):')
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


def newton_raphson_x2(x0):
    print('Newton-Raphson for f(x2):')
    print('Check for the multiple root x2 = 2 not having quadratic convergence(if false for all the approximations approaching our root then proven)')
    k = 1
    xprev = x0
    x = xprev - f(xprev)/df(xprev)
    err = x - xprev
    while abs(x - xprev) > abs(x) * 0.000005:  # we need one more iteration to get our wanted precision
        errprev = err
        xprev = x
        x = xprev - f(xprev)/df(xprev)
        k = k + 1
        err = x - xprev
        print('abs(err / errprev)) = ', abs(err / errprev))  # it's nearly constant so we have nearly linear convergence
    errprev = err
    xprev = x
    x = xprev - f(xprev) / df(xprev)
    k = k + 1
    err = x - xprev
    print('abs(err / errprev)) = ', abs(err / errprev))
    print('k = ', k)
    return x


def secant_x1(x0, x1):
    k = 1
    a = x0
    b = x1
    c = b - (f(b) * (b-a)) / (f(b) - f(a))
    while abs(c - b) > abs(c) * 0.00001:
        a = b
        b = c
        c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
        k = k + 1
    print('Secant for f:')
    print('k = ', k)
    return c


def secant_x2(x0, x1):
    k = 1
    a = x0
    b = x1
    c = b - (f(b) * (b-a)) / (f(b) - f(a))
    while abs(c - b) > abs(c) * 0.000005:
        a = b
        b = c
        c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
        k = k + 1
    a = b
    b = c
    c = b - ((f(b) * (b - a)) / (f(b) - f(a)))
    k = k + 1
    print('Secant for f:')
    print('k = ', k)
    return c


# Bisection

# interval in which Bolzano's theorem must hold true, hence (0, 1.5) and (1.5, 3)
print('x1 = ', bisection(0, 1.5), '\n')
print('x2 = ', bisection(1.5, 3), '\n')


# Newton-Raphson

# choose an interval in which f'(x),f''(x)!=0 so (0, 1) but for the 2nd root, the root itself has f'(x),f''(x)=0
print('f(0)*d2f(0) = ', f(0)*d2f(0), 'f(3)*d2f(3) = ', f(3)*d2f(3), '\n')
# since f(0)*d2f(0)>0 we choose as x0=0 and for the other root since f(3)*d2f(3)>0 we choose x0=3
print('x1 = ', newton_raphson_x1(1), '\n')
print('x2 = ', newton_raphson_x2(3), '\n')
# we can see from abs(f(x)/f(xprev)) that its values for iterations of x1 follow: (10^−3, 10^-6) at the end which makes
# it have quadratic convergence unlike the root x2 which goes almost linearly

# the roots for which the newton raphson method does not have quadratic convergence are pretty clear. In this example
# we can see that the root 2 has f'(2)=0 (in other words we have at least a double root) and therefore with each step
# of the iteration the next solution will have a derivative which gets closer and closer to 0 together with f itself.
# But the derivative will always be larger as it tends to 0 than the function f that also tends to 0 since it has a
# lower degree (since differentiating will lower the degree by 1) and therefore the quotient f(x)/f'(x) will improve
# the guess at a rate worse and worse every time so we won't have quadratic convergence, the convergence will look more
# linear in this case. Experimentally this can be seen by the amount of iterations that are needed for the first root
# with k=5 as compared to the 2nd with k=29 and as another experiment we can see that newton-raphson_g needs 29
# iterations (with g(x)=(x-2)^3 which has triple root at 0). We can also examine
# the fact that the rate between the change of the errors when using newton-raphson for x2 is nearly constant (so we end
# up with more iterations) compared to x1 in which with each passing iteration the rate between the errors gets smaller
# end smaller, meaning that we are approaching it in at least a quadratic rate.

# f(xk)≈f′(x∗)(xk−x∗). Therefore, assuming f′(x∗)≠0, we get |f(xk+1)| / |f(xk)| ≈ |xk+1−x∗| / |xk−x∗|,
# and we'll observe quadratic convergence by looking at the ratios |f(xk+1)/f(xk)| in the roots where f'(x0)!=0
# and for the ones where f'(x0) == 0 we'll just examine the error rates to see if they improve by a near constant factor
print('Also df = ', df(2), '\n')
print('Experiment with function g: ')
print('x = ', newton_raphson_g(3), '\n')


# Secant

# choose two intervals in which bolzano's theorem holds true, (0, 1.5) and (1.5, 3)
print(secant_x1(0, 1))
print(secant_x2(1.5, 3))
