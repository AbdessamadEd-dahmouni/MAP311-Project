#Question 10 : Temps de calcul nécessaire à simuler 10^6 simulations
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.special as ss
import time
nbs=1000000 #nombre de simulations
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

start_time=time.time()
Y=[Gamma_AD() for i in range(nbs)]
print('le temps de calcul est :',time.time()-start_time)