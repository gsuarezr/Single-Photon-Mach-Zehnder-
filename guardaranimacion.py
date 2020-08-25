import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from sympy.plotting import plot3d,plot
from mpl_toolkits.mplot3d import Axes3D
from sympy.utilities.autowrap import ufuncify
from sympy.utilities.lambdify import lambdify
jones1= Matrix([0,1])
jones0= Matrix([1,0])

theta=Symbol('theta',real=True)
beta=Symbol('beta',real=True)
def N_mach(N,k):
    if k==1:
        lista=list()
        K=symarray('theta', N+1,real=True)
        for i in range(1,N+1):
            l=i
            BS=Matrix([[cos(K[i]),I*sin(K[i])],[I*sin(K[i]),cos(K[i])]])
            lista.append(BS)
    else:
        t=np.pi/(2)
        lista=Matrix([[cos(t/N),I*sin(t/N)],[I*sin(t/N),cos(t/N)]])
    return lista

def N_mach_zehnder(N,k):
    K=N_mach(N,k)
    jones1= Matrix([0,1])
    jones0= Matrix([1,0])
    M1=Matrix([[1,0],[0,1]])
    A=Matrix([[ 1,0],[0,beta*(cos(theta)+I*sin(theta))]])
    r=jones0
    if k==1:
        for i in range(0,N):
            if i==N-1:
                r=simplify(K[i]*r)
            else:
                r=simplify(M1*A*K[i]*r)
    else:
        r=K*(M1*A*K)**(N-1)*r
    return r
        

def N_Interferomemer(N,k):
    r=N_mach_zehnder(N,k)
    p1=jones0.T*r
    p1=p1[0]*conjugate(p1[0])
    p1=re(p1)
    p2=jones1.T*r
    p2=p2[0]*conjugate(p2[0])
    p2=re(p2)
    pabs=1-p1-p2
    K=symarray('theta', N+1,real=True)
    return pabs,p1,p2,K

def animacion(N):
    AA=[]
    BB=[]
    CC=[]
    x = np.linspace(-np.pi/2, 3*np.pi/2, 50)
    y = np.linspace(0, 1,50)
    X, Y = np.meshgrid(x, y)
    jones1= Matrix([0,1])
    jones0= Matrix([1,0])
    M1=Matrix([[1,0],[0,1]])
    A=Matrix([[ 1,0],[0,beta*(cos(theta)+I*sin(theta))]])
    r=jones0
    for i in range(2,N+1):
        K=N_mach(N,2)
        if i==N-1:
            r=simplify(K*r)
        else:
            r=simplify(M1*A*K*r)
        pd1=jones0.T*r
        pd1=pd1[0]*conjugate(pd1[0])
        pd1=re(pd1)
        pd2=jones1.T*r
        pd2=pd2[0]*conjugate(pd2[0])
        pd2=re(pd2)
        pabs=1-pd1-pd2
        pd1=simplify(pd1)
        f=ufuncify((theta,beta),pd1)
        Z= f(X, Y)
        pd2=simplify(pd2)
        f2=ufuncify((theta,beta),pd2)
        Z2= f2(X, Y)
        pabs=simplify(pabs)
        f3=ufuncify((theta,beta),pabs)
        Z3 = f3(X, Y)
        AA.append(Z)
        BB.append(Z2)
        CC.append(Z3)
        print(i)
        np.save('nuevoPD1hastan'+str(i)+'.npy',AA)
        np.save('nuevoPD2hastan'+str(i)+'.npy',BB)
        np.save('nuevoPabshastan'+str(i)+'.npy',CC)
    return AA,BB,CC

