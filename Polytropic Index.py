import numpy as np
import matplotlib.pyplot as plt



n=float(input('Enter the polytropic index n : '))
Rstar=float(input('Enter the radius of the star (in Solar radius) '))
Mstar=float(input('Enter the mass of the star (in Solar mass) '))



if (n<0 or n>5):
    print("invalid n")
    exit()

Pi=3.141592
G=6.67e-11
k_B=1.38e-23
m_H=1.67e-27
mu=0.62
Msun=1.989e30
Rsun=7.0e8
Lsun=3.86e26

Rstar=Rstar*Rsun
Mstar=Mstar*Msun  
   
Nmax=100000

x=np.empty(Nmax)
y=np.empty(Nmax)
z=np.empty(Nmax)

# here x = ξ, 
#y = D(ξ), 
#z = dD /d ξ


y[0]=1.0
z[0]=0.0
x[0]=1.0e-6


dx = 0.001

i=0    
while y[i] > 0 and i < Nmax-1: 
    i=i+1   
    z[i]=z[i-1] - ( 2.0*z[i-1]/x[i-1] + y[i-1]**n )*dx
    y[i] = y[i-1] + z[i]*dx 
    x[i] = x[i-1] + dx
    ylast = y[i-1] 
    zlast = z[i-1]

if n>0:      
    x_1 = x[i] 
    z_1 = zlast
    
    a = -x_1/3.*1.0/z_1
    b = 1.0/((n+1)*x_1*(-z_1))
    c = 1.0/(4*Pi*(n+1.)*(-z_1)**2)
    rho_c = a*3.0*Mstar/(4.*Pi*Rstar**3)
    #T_c = b*G*Mstar*mu*m_H/(Rstar*k_B)
    P_c = c*G*Mstar**2/Rstar**4
    
    K = (4.*Pi)**(1./n)*G/(n+1)*x_1**(-(n+1.)/n)*(-z_1)**((1.-n)/n)*Mstar**((n-1.)/n)*Rstar**((3.-n)/n)
    alpha=np.sqrt((n+1.)*K*rho_c**(1./n-1.)/(4.0*Pi*G))
    x2 = np.empty(i)
    y2 = np.empty(i)
    rho = np.empty(i)
    Pres = np.empty(i)
    #Temp = np.empty(i)
    #L = np.empty(i)
    #Ltot = np.zeros(i) 
    
    radius = np.empty(i)
    
    for j in range (0,i): 
        x2[j]=x[j] 
    
        y2[j]=y[j] 
        rho[j]=rho_c*y[j]**n 
        Pres[j]=P_c*y[j]**(n+1) 
    #    Temp[j]=T_c*y[j] 
        radius[j]=alpha*x[j]/Rsun  
    #    L[j]=4.*Pi*(radius[j]*Rsun)**2*(dx*alpha)*2.6E-37*rho[j]**2*Temp[j]**4.
        
    #for j in range (1,i-1): 
    #    Ltot[j]=Ltot[j-1]+L[j]/Lsun


fig,ax=plt.subplots(figsize=(10,4))
    
plt.subplot(121)
plt.plot (x2,y2,linewidth=3)
plt.title (r'$D_n$ vs $\xi$',fontsize=12)
plt.xlabel(r'$\xi$',fontsize=12)
plt.ylabel(r'$D_n$',fontsize=12)


plt.subplots_adjust(hspace=.5)
"""
plt.subplot(212)
plt.plot(radius,rho,lw=3)
plt.xlabel('Radius ($R_{sun}$)',fontsize=12)
plt.ylabel('Density ($kg/m^3$)',fontsize=12)
plt.title('Density profile',fontsize=12)

plt.subplots_adjust(wspace=.5)
"""
if (n>0):
    plt.subplot(122)
    plt.plot (radius,Pres,linewidth=3)
    plt.xlabel('Radius ($R_{sun}$)',fontsize=12)
    plt.ylabel('Pressure ($N/m^2$)',fontsize=12)
    plt.title('Pressure profile',fontsize=12)
    #plt.ylim(-1e5,1e13)

fig.text(0.38,0.95,'polytropic index n = '+str(n),fontweight='bold',fontsize=16)
plt.show()


def le_soln(n):

    x=np.empty(Nmax)
    y=np.empty(Nmax)
    z=np.empty(Nmax)
    
    y[0]=1.0
    z[0]=0.0
    x[0]=1.0e-6
    
    dx = 0.001
    
    i=0    
    while y[i] > 0 and i < Nmax-1: 
        i=i+1   
        z[i]=z[i-1] - ( 2.0*z[i-1]/x[i-1] + y[i-1]**n )*dx
        y[i] = y[i-1] + z[i]*dx 
        x[i] = x[i-1] + dx
        
    x2 = [x[j] for j in range(i)]
    y2 = [y[j] for j in range(i)]

    plt.plot (x2,y2,label='n='+str(n))

le_soln(0)
le_soln(1)
le_soln(2)
le_soln(3)
le_soln(4)
le_soln(5)

plt.title (r'$D_n$ v s $\xi$')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$D_n$')
plt.xlim(-1,15)
plt.ylim(0,1.02)
plt.legend(loc='upper right')
plt.show()


"""
plt.plot (radius,Temp,'b.')
plt.xlabel('Radius (Rsun)')
plt.ylabel('Temperature (K)')
plt.title('Temperature versus radius')
plt.text(1,1,'n = '+str(n)+' polytrope',color='blue',fontsize=16)
plt.show()

plt.plot (radius,Ltot,'b.')
plt.xlabel('Radius (Rsun)')
plt.ylabel('Luminosity (Lsun)')
plt.title('Luminosity versus radius')
plt.text(1,0,'n = '+str(n)+' polytrope',color='blue',fontsize=16)
plt.show()
"""
