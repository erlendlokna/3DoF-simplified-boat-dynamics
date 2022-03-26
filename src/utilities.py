import numpy as np

#Constant parametes 
sigma_0 = 1000              #water density [kg/m^2] 
sigma = 500                 #Ship density [kg/m^2] 
R = 10                      #Ship radius [m]
A_s = 1/2 * np.pi * R**2    #Ship cross section 
m = A_s * sigma             #Ship's mass per length unit (in z direcion)
A_0 = sigma*np.pi*R**2 / (2*sigma_0)
I_c = 1/2 * m * R**2 * (1-(32)/(9*np.pi**2)) #Ships moment of intertia      #/(9*np.pi**2) Legg til for riktig I_c
g = -981/100                #Gravitational acceleration9*np.pi**2))   
h = 0.42* R                 #Distanve from center of mass to physical senter
w_0 = np.sqrt(m*abs(g)*h/I_c)


def newton(f, df, x0, m = 0, tol=1.e-8, max_iter=30):
    ''' Solve f(x)=0 by Newtons method
        The output of each iteration is printed
        Input:
        f, df:   The function f and its derivate f'.
        x0:  Initial values
        tol: The tolerance
      Output:
        The root and the number of iterations
    '''
    x = x0
    #print(f"k ={0:3d}, x = {x:18.15f}, f(x) = {f(x):10.3e}")
    for k in range(max_iter):
        fx = f(x,m)
        if abs(fx) < tol:           # Aksepterer løsningne 
            break 
        x = x - fx/df(x)            # Newton-iterasjon
        #print(f"k ={k+1:3d}, x = {x:18.15f}, f(x) = {f(x):10.3e}")
    return x, k+1

def f3(x, m=0):
    sigma_m = m/(1/2 *np.pi*R**2)
    return x - np.sin(x) -np.pi*((sigma+sigma_m)/sigma_0)

def df3(x):
    return 1-np.cos(x)

def get_beta(m=0): #Hvor m=load mass. Brukes til oppgave 2
    beta, its = newton(f3, df3, 1 , m, tol=1.e-8, max_iter=30)
    return beta

def Y_MB(gamma=get_beta(0)):
    '''
    Input: Gamma (sector angle)
    Output: Center of gravity of repressed water
    '''
    return R * (4*(np.sin(gamma/2)**3))/(3*(gamma-np.sin(gamma)))

