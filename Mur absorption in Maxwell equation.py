import numpy as np
import matplotlib.pyplot as plt
import math as m

#Parameters of simulation
dx=0.05
dt=0.025
T=20
L1=4
Nx=int((L1+0)/dt)+1
Nt=int((T+0)/dx)+1
imp0=377
x=np.linspace(0,L1,Nx)
t=np.linspace(0,T,Nt)
dx=x[1]-x[0]
dt=t[1]-t[0]
Sc=dt/dx
print(Sc,'Curant namber')
plt.ion()
Hy=np.zeros(Nx)
Ez=np.zeros(Nx)
#Parameters of normal distribution (Source)
zigma=0.6
delay=3
#parameters of medium
eps=4
vp=1/(np.sqrt(eps))
for tn in t:
    Hy[-1] = Hy[-2]+(vp*dt-dx)*(Hy[-2]+(Ez[Nx-1]-Ez[Nx-2])*Sc/(imp0)-Hy[-1])/(vp*dt+dx)
    Hy[:-1] += (Ez[1:] - Ez[:-1])*Sc / (imp0)
    Ez[0] = Ez[1]+(vp*dt-dx)*(Ez[1]+(Hy[1]-Hy[0])* Sc*(imp0)/eps-Ez[0])/(vp*dt+dx)
    Ez[1:] += (Hy[1:] - Hy[:-1]) * Sc*imp0/eps
    # source
    Ez[int(Nx/2)-1] += 10**(-3) * m.exp(-((tn-delay)/zigma)**2)
#ploting in time
    plt.plot(x, Ez)
    plt.plot(x, Hy*imp0)
    plt.pause(0.0001)
    plt.clf()
#Last frame
plt.ioff()
plt.plot(x, Ez)
plt.plot(x, Hy*imp0)
plt.show()