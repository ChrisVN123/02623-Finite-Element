import numpy as np
from scipy.sparse import diags
import time
import matplotlib.pyplot as plt

f = 1
psi = 1
x = np.linspace(0,1,1000)
epsilons = [1,0.01,0.0001]
for eps in epsilons:
    ys = 1/psi*((1+(np.exp(psi/eps)-1)*(x)-np.exp(x*(psi/eps)))/(np.exp(psi/eps)-1))
    plt.plot(x,ys,label = f"Epsilon= {eps}")
