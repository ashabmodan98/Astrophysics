import numpy as np
import matplotlib.pyplot as plt


n=float(input('Enter the polytropic index n (0<=n<=5) : '))
Rstar=float(input('Enter the radius of the star (in Solar radius) '))
Mstar=float(input('Enter the mass of the star (in Solar mass) '))

if 0<=n<0.5:m=0
if 0.5<=n<1.:m=2
if 1.<=n<1.5:m=4
if 1.5<=n<2.:m=5
if 2.<=n<2.5:m=6
if 2.5<=n<3.:m=8
if 3.0<=n<=5.0:m=9

Mstar=2
Rstar=1
Pi=3.141592
G=6.67e-11
k_B=1.38e-23
m_H=1.67e-27
mu=0.62
Msun=1.989e30
Rsun=7.0e8
Lsun=3.86e26


k=10**m
K=k**(n/(n+1))

Rstar=Rstar*Rsun
Mstar=Mstar*Msun  

Nmax=int(Rstar/1000.)
P = np.empty(Nmax)
M = np.empty(Nmax)
R = np.empty(Nmax)

P_d = lambda p,m,r : G*m*(p**(n/(n+1)))/(K* r*r)
M_d = lambda p,r : 4.0*Pi*r*r*(p**(n/(n+1)))/K
 
R[0]=Rstar
P[0]=1.0
M[0]=Mstar


dR=1000.0
i=0
for i in range(Nmax-1):

    p,m,r =P[i],M[i],R[i]
    
    M[i+1]=M[i] - dR *M_d(p,r)
    P[i+1]=P[i] + dR *P_d(p,m,r)
    R[i+1]=R[i] - dR
    #if(M[i+1]>M[i]  or R[i+1]<1.0 or P[i+1]<P[i]): break
    if(M[i+1]<1e-8  or R[i+1]<1e-8 ): break

R=R/Rsun
"""
R=[R[j] for j in range(i)]
P=[P[j] for j in range(i)]
M=[M[j] for j in range(i)]
"""
plt.plot(R[:i],P[:i],linewidth=3)
#plt.xlim(-.2,2)
#plt.ylim(-1e5,1e15)
#plt.scatter(R,M)
#plt.show

plt.xlabel('Radius ($R_{sun}$)',fontsize=12)
plt.ylabel('Pressure ($N/m^2$)',fontsize=12)
plt.title('Pressure profile',fontsize=12)
plt.text(1,1,'n='+str(n),fontsize=16,fontweight='bold')
plt.show()
