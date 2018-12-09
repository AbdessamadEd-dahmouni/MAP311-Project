#Question 7 : simuler Gamma(alpha,1) par la méthode de Cheng and Feast
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.stats as st
import scipy.special as ss
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
nbs=10000 #nombre de simulations
a=2.5 #alpha dans l'énoncé
lbCF=(a-1/(6*a))/(a-1) #lambda_CF dans l'énoncé
lb_s=np.power(((a+1)/(a-1)),(a+1)/2) #lambda soulignée dans l'énoncé

def estDansA(X,lb) :    #teste si X appartient à A
    if ( (2/(a-1))*np.log(X[0])-np.log(lb*X[1]/X[0])+lb*X[1]/X[0]-1<=0 ) :
        return 1
    return 0

def X_T(lb):  # fonction qui donne une valeur aléatoire de X_T
    X=[random(),random()]
    while(estDansA(X,lb)!=1):
        X=[random(),random()]
    return X
def Gamma_CF(lb): #Simulation de la loi Gamma(alpha,1)
    X=X_T(lb)
    return (a-1)*lb*X[1]/X[0]

debut=st.gamma.ppf(0,a)
fin=st.gamma.ppf(0.99,a)
num_bins=np.linspace(debut,fin,100)
Y1=st.gamma.pdf(num_bins,a) #densité de Gamma(5/2,1)

def Simulation(lb,nbs):
    Y2=[Gamma_CF(lb) for i in range(nbs)]
    nbs, bins, patches = plt.hist(Y2, num_bins, normed=1, facecolor='green', alpha=0.5,label=r'L\textquoteright histogramme de $\Gamma(\alpha,1)$')
    plt.plot(bins,Y1,color='r',label=r'La densit\'e de la loi $\Gamma(\alpha,1)$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\frac{1}{\Gamma(\alpha)}e^{-v}v^{\alpha-1}\mathbf{1}_{\{x>0\}}$')
    if lb==lb_s :
      plt.title(r'Simulation et densit\'e de la loi $\Gamma(\alpha,1)$ pour $\alpha=\frac{5}{2}$ et $\lambda=\underline{\lambda}$')
    if lb==lbCF :
      plt.title(r'Simulation et densit\'e de la loi $\Gamma(\alpha,1)$ pour $\alpha=\frac{5}{2}$ et $\lambda=\lambda^{CF}$')
    plt.legend()
    plt.xlim(debut,fin)
    return

plt.figure(1)
Simulation(lbCF,nbs)
plt.figure(2)
Simulation(lb_s,nbs)

plt.show()

