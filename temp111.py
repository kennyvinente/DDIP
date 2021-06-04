

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(y,t):
    k = 0.3
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()


#################### Sympy

import sympy

#We need to start the pretty-printer to get nicely typeset math
sympy.init_printing()

#In order to do symbolic calculations, we need to create a symbol
x = sympy.Symbol('x')
print(x)


#Sympy allows us to do many mathematical operations that would be tedious by hand. For instance, we can expand a polynomial:
polynomial = (2*x + 3)**4
polynomial.expand()

polynomial

(x**2 + 2*x + 1).factor()


# Sympy knows how to integrate and differentiate
eq = sympy.tan(sympy.log(x**2 + 1))
#print(eq.diff(x))
eq.diff(x)

polynomial.diff(x) # First derivative
polynomial.diff(x, 2) # Second derivative
polynomial.integrate(x) # indefinite integral - note no constant of integration is added
polynomial.integrate((x, 1, 2)) # Note that integrate takes one argument which is a tuple for the definite integral


nonlinear_expression = sympy.sin(x)
sympy.series(nonlinear_expression, x, 2, 7) # taylor expansion in terms of the x variable, around x=2, first order

#To remove the order term use .removeO()
temp = sympy.series(nonlinear_expression, x, 2, 2)
temp.removeO()

#You will also notice that SymPy’s default behaviour is to retain exact representations of certain numbers:
number = sympy.sqrt(2)*sympy.pi
number

#o convert the exact representations above to an approximate floating point representations, use one of these methods.
# sympy.N works with complicated expressions containing variables as well. float will return a normal Python float and
# is useful when interacting with non-sympy programs.
sympy.N(number*x)
float(number)


sympy.plot(x**2 + 1)
solutions = sympy.solve(2*x**2 + 2 - 4)
solutions
solutions[0]

#We can also use sympy.Eq to construct equations
equation = sympy.Eq(2*x**2 + 2, 4)
equation

#The roots function will give us the multiplicity of the roots as well.
sympy.roots(equation)

#We can also solve systems of equations by passing a list of equations to solve and asking for a list of variables to solve for
x, y = sympy.symbols('x, y')
sympy.solve([x + y - 2,
             x - y - 0], [x, y])

#This even works with symbolic variables in the equations
a, b, c = sympy.var('a, b, c')
solution = sympy.solve([a*x + b*y - 2,
                        a*x - b*y - c], [x, y])
solution


import sympy
sympy.init_printing()
D, d, H, h = sympy.symbols('D, d, H, h', real=True)
R = D/2
r = d/2
Hprime = sympy.sqrt(H**2 - (R - r)**2)  # Pythagoras

#The radius changes linearly from the small one to the large one:
radius = r + h/Hprime*(R - r)

A = sympy.pi*radius**2
V = sympy.integrate(A, (h, 0, h))

V
print(V)


Vsymb = sympy.symbols('V', real=True)
hV = sympy.solve(Vsymb - V, h)
print(hV)


# Sympy can solve some differential equations analytically:
h = sympy.Function('h')  # This is how to specify an unknown function in sympy
t = sympy.Symbol('t', positive=True)

Fin = 2
K = 1
A = 1

Fout = K*h(t)
Fout

#We use .diff() to take the derivative of a function.
de = h(t).diff(t) - 1/A*(Fin - Fout)
de



#Here we calculate the general solution. Notice this equation just satisfies the original differential equation when we
# plug it in, we don’t have specific values at points in time until we specify boundary conditions.

solution = sympy.dsolve(de)
solution

#We need a name for the constant of integration which Sympy created. Expressions are arranged as trees with the
# arguments as elements. We can navigate this tree to get the C1 element:

C1 = solution.rhs.args[1].args[0]

#We can find the value of the constant by using an initial value:
h0 = 1
constants = sympy.solve(solution.rhs.subs({t: 0}) - h0, C1)
constants

import matplotlib.pyplot as plt
sympy.plot(solution.rhs.subs({C1: constants[0]}), (t, 0, 10))


#When the boundary conditions of differential equations are specified at t=0, this is known as an Initial Value Problem
# or IVP. We can solve such problems numerically using scipy.integrate.solve_ivp.

import numpy, scipy.integrate
Fin = 2
def dhdt(t, h):
    """Function returning derivative of h - note it takes t and h as arguments"""
    Fout = K*h
    return 1/A*(Fin - Fout)

#solve_ivp will automatically determine the time steps to use, integrating between the two points in tspan:
tspan = (0, 10)
sol = scipy.integrate.solve_ivp(dhdt, tspan, [h0])

#We’ll need a smooth set of time points to evaluate the analytic solution
tsmooth = numpy.linspace(0, 10, 1000)
hanalytic = 2 - numpy.exp(-tsmooth)

plt.plot(sol.t, sol.y.T, 'o-', label='solve_ivp solution')
plt.plot(tsmooth, hanalytic, label='Analytic solution')
plt.legend()


sol = scipy.integrate.solve_ivp(dhdt, tspan, [h0], t_eval=tsmooth)
plt.plot(tsmooth, sol.y.T)
plt.plot(tsmooth, hanalytic, '--')


import scipy.integrate
def Fin(t):
    """ A step which starts at t=2 """
    if t < 2:
        return 1
    else:
        return 2

def dhdt(t, h):
    Fout = K*h
    return 1/A*(Fin(t) - Fout)

import matplotlib.pyplot as plt

t, s = sympy.symbols('t, s')
a = sympy.symbols('a', real=True, positive=True)

f = sympy.exp(-a*t)
f

tspan = (0, 10)

sol = scipy.integrate.solve_ivp(dhdt, tspan, [h0])
smoothsol = scipy.integrate.solve_ivp(dhdt, tspan, [h0], t_eval=tsmooth)

plt.plot(sol.t, sol.y.T, 'o')
plt.plot(smoothsol.t, smoothsol.y.T)


sol = scipy.integrate.solve_ivp(dhdt, tspan, [h0], max_step=0.1)
smoothsol = scipy.integrate.solve_ivp(dhdt, tspan, [h0], t_eval=tsmooth, max_step=0.1)

plt.plot(sol.t, sol.y.T, 'o')
plt.plot(smoothsol.t, smoothsol.y.T)


def odeintdhdt(h, t):
    """Odeint expects a function with the arguments reversed from solve_ivp"""
    return dhdt(t, h)

#he order in which the initial values and the input times are specified is different.
# Also, unlike solve_ivp, you always give odeint the times you want the results at.

odeinth = scipy.integrate.odeint(odeintdhdt, h0, tsmooth)
plt.plot(sol.t, sol.y.T, 'o')
plt.plot(tsmooth, odeinth)


a = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1
a

a==1
b = 0.125 + 0.125 + 0.125 + 0.125 + 0.125 + 0.125 + 0.125 + 0.125
b

b == 1


import decimal

a = decimal.Decimal('0.1')

a = 1/decimal.Decimal(10)

sum(a for i in range(10))

decimal.Decimal('1.0')

sum(0.1 for i in range(10))


import sympy
sympy.init_printing()

b = sympy.Rational('0.1')

b

b = 1
c = 10

a = b/c

type(a)


sympy.nsimplify(a)

x = sympy.Symbol('x')

expr = sympy.sqrt(3)*x
expr

sympy.simplify(expr**2)

sympy.simplify(sympy.N(expr, 3)**2)






sympy.integrate(f*sympy.exp(-s*t), (t, 0, sympy.oo))

sympy.laplace_transform(f, t, s)


from numpy import genfromtxt
my_data = genfromtxt('C:\\Users\\User\\Documents\\git\julia\\data_expAC.csv', delimiter=',')
















