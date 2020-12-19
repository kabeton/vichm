import numpy as np
import matplotlib.pyplot as plt

def rhs(t, coords):
    ans = np.zeros(4)
    r = (coords[0]**2 + coords[2]**2)**(1/2)
    ans[0] = coords[1]
    ans[1] = -1 * coords[0] * 39.9533e13 / r**(3)
    ans[2] = coords[3]
    ans[3] = -1 * coords[2] * 39.9533e13 / r**(3)
    return ans

def rk4(f,y_0,x): 
    crmar = 0
    T = 0; tgot = 0
    n = len(x)       
    y = np.zeros((n,len(y_0)))
    h = x[1]-x[0]
    y[0] = np.copy(y_0)
    for i in range(1,n):
        k_1 = f(x[i-1], y[i-1])
        k_2 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k_1)
        k_3 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k_2)
        k_4 = f(x[i-1] + h, y[i-1] + h*k_3)
        y[i] = y[i-1] + h*(k_1 + 2*k_2 + 2*k_3 + k_4)/6
        if ((abs(y[i][0] - y_0[0])) < 1 and abs(y[i][2]) < 10) and (not tgot):
            T = x[i]
            tgot += 1
        if (y[i][0]**2 + y[i][2]**2)**(1/2) < 6380:
            crmar += 1
    return y, crmar, T

rc = 1e4
gamma = 6.67e-11
u = 0
M = 5.99e24
vc = (gamma*M/rc)**(1/2)
init = np.array([rc, 0, 0, vc - u])

crmar = 0
t = np.linspace(0, 1, 1000)
while(not crmar):
    y, crmar, T = rk4(lambda t,c:rhs(t, c), init, t)
    u += 10
    init[-1] -= u
    x = [y[i][0] for i in range(len(y))]
    a = (abs(min(x)) + abs(max(x)))/2
    t2 = 2*3.1415926*a**(3/2)/(gamma * M)**(1/2)
    print(T / t2)
    d = [y[i][2] for i in range(len(y))]

print(u)

x = [y[i][0] for i in range(len(y))]
d = [y[i][2] for i in range(len(y))]

earth = plt.Circle((0, 0), 6380)
fig, ax = plt.subplots()
ax.add_artist(earth)
ax.set_aspect(1)
plt.scatter(x, d)
plt.grid()
plt.show()


