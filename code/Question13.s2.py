#Question 13 : Calcul de l'espérance E(exp(-uY)) par la méthode de Monte-Carlo
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.stats as st
import scipy.special as ss
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

d=float(1)
nu=0.5
#nu=1.5
n=1000
def estDansA(X,a) : #teste si X appartient à A, l'ensemble A dépend du paramètre a=alpha
    if(a>1) :
        lb=(a-1/(6*a))/(a-1) #lambda dans la méthode CF
        if ( (2/(a-1))*np.log(X[0])-np.log(lb*X[1]/X[0])+lb*X[1]/X[0]-1<=0 ) :
            return 1
        return 0
    else :
        lb=(a+np.e)/np.e #lambda dans la méthode AD
        if ( lb*X[0]<=1 ):
            if( X[1]<=np.exp(-np.power(lb*X[0],1/a)) ) :
                return 1
            else :
                return 0
        else:
            if( X[1]<=np.power(-np.log(lb*(1-X[0])/a),a-1) ) :
                return 1
            else :
                return 0

def X_T(a):  # fonction qui donne une valeur aléatoire de X_T
    X=[random(),random()]
    while(estDansA(X,a)!=1):
        X=[random(),random()]
    return X

def Gamma_CF_AD(a): #Simulation de la loi Gamma(alpha,1/2)
    X=X_T(a)
    if(a>1) :
        lb=(a-1/(6*a))/(a-1)
        return 2*(a-1)*lb*X[1]/X[0]
    else :
        lb=(a+np.e)/np.e
        if ( lb*X[0]<=1 ) :
            return 2*np.power(lb*X[0],1/a)
        else:
            return -2*np.log(lb*(1-X[0])/a)

def Chi2(): #Simulation de chi^2 par conditionnement avc N une loi de poisson de paramètre d/2
    N=np.random.poisson(d/2)
    return Gamma_CF_AD(N+nu/2)
    
def EsperanceMC(u,n):
    E=0
    for i in range(n) :
        E+=np.exp(-u*Chi2())
    return E/n

def EsperanceTH(u) :
    return np.power(1/(2*u+1),nu/2)*np.exp(d/(4*u+2)-d/2)
    

X=np.linspace(0,1,100)
Y=[EsperanceMC(x,n) for x in X] #(Y=EsperanceMC(X,n) garde apparemment les mêmes valeurs de chi2() calculées durant le premier appel de cette
                                  #fonction et les utilise pour le calcul de tous les Y[i])
Z=EsperanceTH(X)
Ypls=Y+1.96*(EsperanceTH(2*X)-EsperanceTH(X)**2)/np.sqrt(n) #Y+1.96*\sigma(chi2)/racine(n)
Ymns=Y-1.96*(EsperanceTH(2*X)-EsperanceTH(X)**2)/np.sqrt(n)
plt.plot(X,Y,color='r',label=r'$\mathbb{E}(e^{-uY})$ par Monte-Carlo')
plt.fill_between(X,Ypls,Ymns,where=Ypls >= Ymns,facecolor='green',alpha=0.3,interpolate=True) # la bande constituée par les intervalles de confiance
plt.plot(X,Z,color='blue',label=r'$\mathbb{E}(e^{-uY})$ th\'eorique')
plt.axis('equal')
plt.xlabel('u')
plt.ylabel(r'$\mathbb{E}(e^{-uY})$')
plt.title(r'Calcul de l\textquoteright esp\'erance $\mathbb{E}(e^{-uY})$ avec $n=$'+str(n))
plt.legend()
plt.show()

    
    