#Question 10 : Intervalle de confiance de la probabilité de rejet par la méthode de Monte-Carlo
import numpy as np
import matplotlib.pyplot as plt
from random import *
import scipy.special as ss
nbs=10000 #nombre de simulations
a=0.5 #alpha
lb=(a+np.e)/np.e #lambda dans la question 9
pRejet=1-a*ss.gamma(a)/lb
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

n=0
X=[random(),random()]
for i in range(nbs):
    if estDansA(X)==0 :
        n+=1
    X=[random(),random()]
print("l'intervalle de confiance pour la probabilité de rejet est ["+str(n/nbs-1/np.sqrt(nbs))+","+str(n/nbs+1/np.sqrt(nbs))+"]")
print("la valeur théorique est : ",pRejet)