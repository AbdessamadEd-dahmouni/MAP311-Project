#Question 13 : simulation chi² par le conditionnement avec une loi de poisson 
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.stats as st
import scipy.special as ss
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

nbs=10000
d=float(1)
nu=1.5

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

#def Chi2():    #Simulation de chi^2 par la méthode de laquestion 15
#    return (np.random.normal(0.0,1.0,1)[0]+np.sqrt(d))**2+Gamma_CF_AD((nu-1)/2)

debut=st.ncx2.ppf(0.01,nu,d)
fin=st.ncx2.ppf(0.99, nu, d)
X=np.linspace(debut,fin,100)
plt.plot(X, st.ncx2.pdf(X, nu, d),'r-', label=r'Densit\'e de $\chi^2(\nu,d)$')
Y=[Chi2() for i in range(nbs)]
nbs, bins, patches = plt.hist(Y, X, normed=1, facecolor='yellowgreen', alpha=0.5,label=r'L\textquoteright histogramme de $\chi^2(\nu,d)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$p_{\nu,d}(x)$')
plt.title(r'Simulation et densit\'e de la loi $\chi^2(\nu,d)$')
plt.axis([debut,fin,0,0.5])
plt.legend()
plt.show()



