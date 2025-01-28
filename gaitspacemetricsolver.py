#!/usr/bin/env python
# coding: utf-8

# In[1]:


# adapted from Simon Wilshin's toroidalGeometry.py code, included as Supplementary Information in:
#  Wilshin, S., Shamble, P. S., Hovey, K. J., Harris, R., Spence, A. J., Hsieh, S. T. (2018). 
#  Limping following limb loss increases locomotor stability.
#  Journal of Experimental Biology, 221(18)jeb174268.

import sympy
import numpy as np
import math
p0,p1,p2,p3,p4,t0,t1,t2,t3,t4,t5,psi = sympy.symbols('p0 p1 p2 p3 p4 t0 t1 t2 t3 t4 t5 psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,p2-t2+t3,p3-t3+t4,p4-t4+t5,t0+t1+t2+t3+t4+t5-psi],t0,t1,t2,t3,t4,t5)
g5T = np.zeros((5,5))
eta = [p0,p1,p2,p3,p4]
for j in range(len(eta)):
    for k in range(len(eta)):
        g5T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g5T)      # 5 dimensional, 6 legs 


# In[2]:


p0,p1,p2,t0,t1,t2,t3,psi = sympy.symbols('p0 p1 p2 t0 t1 t2 t3 psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,p2-t2+t3,t0+t1+t2+t3-psi],t0,t1,t2,t3)
g3T = np.zeros((3,3))
eta = [p0,p1,p2]
for j in range(len(eta)):
    for k in range(len(eta)):
        g3T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g3T*4)             # 3 dimensional, 4 legs 


# In[ ]:





# In[3]:


p0,p1,p2,p3,p4,p5,t0,t1,t2,t3,t4,t5,t6,psi = sympy.symbols('p0 p1 p2 p3 p4 p5 t0 t1 t2 t3 t4 t5 t6 psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,p2-t2+t3,p3-t3+t4,p4-t4+t5,p5-t5+t6,t0+t1+t2+t3+t4+t5+t6-psi],t0,t1,t2,t3,t4,t5,t6)
g6T = np.zeros((6,6))
eta = [p0,p1,p2,p3,p4,p5]
for j in range(len(eta)):
    for k in range(len(eta)):
        g6T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g6T*7)         # 6 dimensional, 7 legs 


# In[4]:


p0,p1,t0,t1,t2,psi = sympy.symbols('p0 p1 t0 t1 t2  psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,t0+t1+t2-psi],t0,t1,t2)
g2T = np.zeros((2,2))
eta = [p0,p1]
for j in range(len(eta)):
    for k in range(len(eta)):
        g2T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g2T*3)         # 2 dimensional, 3 legs 


# In[5]:


p0,p1,p2,p3,t0,t1,t2,t3,t4,psi = sympy.symbols('p0 p1 p2 p3 t0 t1 t2 t3 t4 psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,p2-t2+t3,p3-t3+t4,t0+t1+t2+t3+t4-psi],t0,t1,t2,t3,t4)
g4T = np.zeros((4,4))
eta = [p0,p1,p2,p3]
for j in range(len(eta)):
    for k in range(len(eta)):
        g4T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g4T*5)        # 4 dimensional, 5 legs 


# In[6]:


p0,p1,p2,p3,p4,p5,p6,t0,t1,t2,t3,t4,t5,t6,t7,psi = sympy.symbols('p0 p1 p2 p3 p4 p5 p6 t0 t1 t2 t3 t4 t5 t6 t7 psi')
eqn = sympy.solvers.solve([p0-t0+t1,p1-t1+t2,p2-t2+t3,p3-t3+t4,p4-t4+t5,p5-t5+t6,p6-t6+t7,t0+t1+t2+t3+t4+t5+t6+t7-psi], \
                          t0,t1,t2,t3,t4,t5,t6,t7)
g7T = np.zeros((7,7))
eta = [p0,p1,p2,p3,p4,p5,p6]
for j in range(len(eta)):
    for k in range(len(eta)):
        g7T[j,k] = np.sum([sympy.diff(eqn[th],eta[j])*sympy.diff(eqn[th],eta[k]) for th in eqn])
        
print(g7T*8)        # 7 dimensional, 8 legs 

# and so on, changing the entries above for pn,tn, gNT, etc. for the correct number of legs


# In[ ]:




