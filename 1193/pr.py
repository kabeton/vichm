import numpy as np
import matplotlib.pyplot as plt

def rhs(t, coords):
    ans = np.zeros(4)
    ans[0] = coords[1]
    ans[1] = t * (coords[0])**(1/2)
    ans[2] = coords[3]
    if coords[0] != 0:
        ans[3] = t*coords[2]/(2*(coords[0])**(1/2))
    else:
        ans[3] = 0
    return ans

def rk4(f, y_0, x): 
    n = len(x)       
    y = np.zeros((n,len(y_0)))
    h = x[1]-x[0]
    y[0] = np.copy(y_0)
    for i in range(1,n):
        k_1 = f(x[i-1],y[i-1])
        k_2 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k_1)
        k_3 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k_2)
        k_4 = f(x[i-1] + h, y[i-1] + h*k_3)
        y[i] = y[i-1] + h*(k_1 + 2*k_2 + 2*k_3 + k_4)/6
    return y

def shoot(f, x, y_0, alpha_0, beta, eps): 
  y = np.zeros((len(x), 4))
  alpha = alpha_0
  while(abs(y[-1][0] - beta) > eps):           
    y = rk4(lambda x,y: f(x,y), np.array([y_0, alpha, 0, 1]), x)
    alpha = alpha - (y[-1][0] - beta) / (y[-1][2])
  return y, alpha

x = np.linspace(0, 1, 1000)
alpha = 0
init = 0
end = 2
eps = 10e-6

y, alpha = shoot(lambda x,y:rhs(x, y), x, init, alpha, end, eps)
print(alpha)
d = [y[i][0] for i in range(len(y))]
plt.scatter(x, d)
plt.grid()
plt.show()
