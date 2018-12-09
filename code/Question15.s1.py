#Question 15 : Calcul de l'espérance E(exp(-uY)) par la méthode de Monte-Carlo
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.stats as st
import scipy.special as ss
import time

d=float(5)       #paramètre de décentralité
rd=np.sqrt(d)
nu=7/float(3)    #nombre de degrés de liberté
nbs=100000         #nombre de simulations

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

def X_T(a):             # fonction qui donne une valeur aléatoire de X_T
    X=[random(),random()]
    while(estDansA(X,a)!=1):
        X=[random(),random()]
    return X

def Gamma_CF_AD(a):     #Simulation de la loi Gamma(alpha,1/2)
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

def Chi2():             #Simulation de chi^2 par conditionnement avc N une loi de poisson de paramètre d/2
    N=np.random.poisson(d/2)
    return Gamma_CF_AD(N+nu/2)

def Chi2Amelioree():    #Simulation de chi^2 par la méthode de la question 15
    return (np.random.normal(0,1,1)[0]+rd)**2+Gamma_CF_AD((nu-1)/2)

start_time=time.time()
Y=[Chi2() for i in range(nbs)]
print('le temps de calcul par la méthode de la question 13 est :',time.time()-start_time)
start_time=time.time()
Z=[Chi2Amelioree() for i in range(nbs)]
print('le temps de calcul par la méthode de la question 15 est :',time.time()-start_time)