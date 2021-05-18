import numpy as np

def fx(x):
    return np.array([1,1])

def G(x):
    return np.array([x[0]**2 + x[1]**2 - 2])

def Gx(x):
    return np.array([[2*x[0], 2*x[1]]])

def fxx(x):
    return np.zeros((2,2))

def Gxx(x,mu):
    return 2*np.eye(2)*mu

x = np.array([0.5,-0.5])
mu = np.array([1])
gamma = 10
Z = np.array([10])
zeta = 0.99995
sigma = 0.1

for i in range(15):
    Lx = np.array(fx(x) + mu.T*Gx(x))
    Lz = np.array(mu.T - gamma/Z)
    Lm = Gx(x).T + Z.T
    Lxx = fxx(x) + Gxx(x,mu)

    A = np.hstack([Lxx, np.zeros((2,1)), Gx(x).T])
    A = np.vstack((A,np.hstack((np.zeros((2)),np.array(mu),Z))))
    A = np.vstack((A,np.hstack((Gx(x),np.ones((1,1)),np.zeros((1,1))))))

    b1 = mu*Z - gamma
    b2 = G(x)+Z
    b3 = np.vstack((Lx.T,b1,b2))
    b = np.array([b3[0][0],b3[1][0],b3[2][0],b3[3][0]])

    y = np.linalg.solve(A,-b)
    if np.linalg.norm(y) <= 10**-8:
        break
    Dx = np.array([y[0],y[1]])
    Dz = np.array([y[2]])
    Dm = np.array([y[3]])

    ap = min(zeta*(-Z[0]/Dz[0]),1) if Dz[0] < 0 else 1
    ad = min(zeta*(-mu[0]/Dm[0]),1) if Dm[0] < 0 else 1

    x = x + ap*Dx
    Z = Z + ap*Dz
    mu = mu + ad*Dm
    gamma = sigma*Z[0]*mu[0]

    print([x, Z, mu, gamma, np.linalg.norm(y)])

for i in range(15):
    A = np.array([[2*mu[0],0,0,2*x[0]],
                  [0,2*mu[0],0,2*x[1]],
                  [0,0,mu[0],Z[0]],
                  [-2*x[0],-2*x[1],-1,0]])
    b = np.array([1+2*x[0]*mu[0],1+2*x[1]*mu[0],mu[0]*Z[0]-gamma,2-x[0]**2-x[1]**2-Z[0]])

    y = np.linalg.solve(A,-b)
    if np.linalg.norm(y) <= 10**-8:
        break
    Dx = np.array([y[0],y[1]])
    Dz = np.array([y[2]])
    Dm = np.array([y[3]])

    ap = min(zeta*(-Z[0]/Dz[0]),1) if Dz[0] < 0 else 1
    ad = min(zeta*(-mu[0]/Dm[0]),1) if Dm[0] < 0 else 1

    x = x + ap*Dx
    Z = Z + ap*Dz
    mu = mu + ad*Dm

    print(x)

    gamma = sigma*Z[0]*mu[0]