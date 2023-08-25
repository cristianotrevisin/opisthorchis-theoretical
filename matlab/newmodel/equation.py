import numpy as np
from scipy.integrate import quad
import math

  


def integrand(s,x,k,theta):
    return x**(k-1)*math.e**(-s/theta)

def probabilite(x,L,lambda_var,k,theta):
    prob = 0;
    for l in range(1,L+1):
        for y in range(1,x+1):
            prob += math.e**(-lambda_var)*lambda_var**x/np.math.factorial(x) * quad(integrand, l-1/2, l+1/2, args=(x,k,theta))[0]/(math.gamma(theta)*theta**k)
    return prob


    
theta = 1 #scale
k = 5 #shape 
lambda_var = 1
L = 16

results = np.zeros((30,))
valeurs = np.zeros((30,))

for i in range(30):
    results[i] = probabilite(i+1,L,lambda_var,k,theta)
    valeurs[i] = i+1

import matplotlib.pyplot as plt


print(results)



