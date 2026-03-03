import math
import matplotlib.pyplot as plt
import numpy as np

def mA(x):
    if x > 15:
        y = 1 / (1 + math.pow(x - 15, -2))
    else:
        y = 0
    return y

def mB(x):
    return 1 / (1 + math.pow(x - 17, 4))

def notmA(x):
    return 1 - mA(x)

#T-norms
def mCTmin(x):
    return min(mA(x), mB(x))

def mCTalg_prod(x):
    return mA(x) * mB(x)

def mCTbounded_prod(x):  # or Lukasiewicz T-norm
    return max(0, mA(x) + mB(x) - 1)

def mCTdrastic_prod(x): #can't compare floats with precision so there's some deformation
    if abs(mB(x) - 1) < 0.0001 * max(abs(mB(x)), 1):
        y = mA(x)
    elif abs(mA(x) - 1) < 0.0001 * max(abs(mA(x)), 1):
        y = mB(x)
    elif mA(x) < 1 and mB(x) < 1:
        y = 0
    return y

#T-conorms
def mDTmax(x):
    return max(mA(x), mB(x))

def mDTalg_sum(x):
    return mA(x) + mB(x) - (mA(x) * mB(x))

def mDTbounded_sum(x):  # or Lukasiewicz T-conorm
    return min(1, mA(x) + mB(x))

def mDTdrastic_sum(x):
    if abs(mB(x)) < 0.0001:
        y = mA(x)
    elif abs(mA(x)) < 0.0001:
        y = mB(x)
    elif mA(x) > 0 and mB(x) > 0:
        y = 1
    return y

#NOT mA AND mB
def mETmin(x):
    return min(notmA(x), mB(x))

def mETalg_prod(x):
    return notmA(x) * mB(x)

def mETbounded_prod(x):  # or Lukasiewicz T-norm
    return max(0, notmA(x) + mB(x) - 1)

def mETdrastic_prod(x):
    if abs(mB(x) - 1) < 0.0001 * max(abs(mB(x)), 1):
        y = notmA(x)
    elif abs(notmA(x) - 1) < 0.0001 * max(abs(notmA(x)), 1):
        y = mB(x)
    elif notmA(x) < 1 and mB(x) < 1:
        y = 0
    return y


x_vals = np.linspace(10, 25, 500)

# Compute y values for each method
y_min = [mCTmin(x) for x in x_vals]
y_alg_prod = [mCTalg_prod(x) for x in x_vals]
y_bounded_prod = [mCTbounded_prod(x) for x in x_vals]
y_drastic_prod = [mCTdrastic_prod(x) for x in x_vals]

y_max = [mDTmax(x) for x in x_vals]
y_alg_sum = [mDTalg_sum(x) for x in x_vals]
y_bounded_sum = [mDTbounded_sum(x) for x in x_vals]
y_drastic_sum = [mDTdrastic_sum(x) for x in x_vals]

y_min1 = [mETmin(x) for x in x_vals]
y_alg_prod1 = [mETalg_prod(x) for x in x_vals]
y_bounded_prod1 = [mETbounded_prod(x) for x in x_vals]
y_drastic_prod1 = [mETdrastic_prod(x) for x in x_vals]

#(a)
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_min, label='Minimum (min)', linewidth=2)
plt.plot(x_vals, y_alg_prod, label='Algebraic Product', linewidth=2)
plt.plot(x_vals, y_bounded_prod, label='Bounded Product (Lukasiewicz)', linewidth=2)
plt.plot(x_vals, y_drastic_prod, label='Drastic Product', linewidth=2)

plt.title('(A)T-Norm Combination Methods')
plt.xlabel('x')
plt.ylabel('Membership Value')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#(b)
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_max, label='Maximum (max)', linewidth=2)
plt.plot(x_vals, y_alg_sum, label='Algebraic Sum', linewidth=2)
plt.plot(x_vals, y_bounded_sum, label='Bounded Sum (Lukasiewicz)', linewidth=2)
plt.plot(x_vals, y_drastic_sum, label='Drastic Sum', linewidth=2)

plt.title('(B)T-Conorm Combination Methods')
plt.xlabel('x')
plt.ylabel('Membership Value')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#(c)
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_min1, label='Minimum (min)', linewidth=2)
plt.plot(x_vals, y_alg_prod1, label='Algebraic Product', linewidth=2)
plt.plot(x_vals, y_bounded_prod1, label='Bounded Product (Lukasiewicz)', linewidth=2)
plt.plot(x_vals, y_drastic_prod1, label='Drastic Product', linewidth=2)

plt.title('(C)')
plt.xlabel('x')
plt.ylabel('Membership Value')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()