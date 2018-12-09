#Qestion 10 : Simulation de la loi Gamma(1/2,1) par la méthode d'Ahrens et Dieter
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.stats as st
import scipy.special as ss
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

nbs=10000 #nombre de simulations
a=0.5 #alpha
lb=(a+np.e)/np.e #lambda dans la question 9
def estDansA(X) :    #Fonction qui teste si X appartient à A
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
            
def X_T():  # fonction qui donne une valeur aléatoire de X_T
    X=[random(),random()]
    while(estDansA(X)!=1):
        X=[random(),random()]
    return X

def Gamma_AD():     #la loi Gamma(alpha,1)
    X=X_T()
    if ( lb*X[0]<=1 ) :
        return np.power(lb*X[0],1/a)
    else:
        return -np.log(lb*(1-X[0])/a)
debut=st.gamma.ppf(0.01,a)
fin=st.gamma.ppf(0.99, a)
num_bins=np.linspace(debut,fin,100)
Y=[Gamma_AD() for i in range(nbs)]
plt.figure(3)
nbs, bins, patches = plt.hist(Y, num_bins, normed=1, facecolor='deepskyblue', alpha=0.5,label=r'L\textquoteright histogramme de $\Gamma(\alpha,1)$')
plt.plot(num_bins,st.gamma.pdf(num_bins,a),color='r',label=r'La densit\'e de la loi $\Gamma(\alpha,1)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{1}{\Gamma(\alpha)}e^{-v}v^{\alpha-1}\mathbf{1}_{\{x>0\}}$')
plt.title(r'Simulation et densit\'e de la loi $\Gamma(\alpha,1)$ pour $\alpha=\frac{1}{2}$')
plt.xlim(debut,fin)
plt.ylim(0,2)
plt.legend()
plt.show()
