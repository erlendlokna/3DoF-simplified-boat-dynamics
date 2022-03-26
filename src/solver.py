from ode import *
from methods import *
from utilities import *
from scipy.integrate import RK45


def solver(x0, xend, y0, loadmass = 0, theta_area=True, fence=False, h=0.01, method=RK4_step, env=ode.Enviroments(), verbose=True):
    """
    Løser dynamikk/ode til båten.
    Input:
        x0: starttid
        xend: slutttid
        y0: initial state
        mass_load: lastens masse
        theta_area: (bool) for om vi skal ta hensyn til areal som funksjon av theta
        fence: (bool) for om det skal være gjerde på båten
        h: stepsize til metoden
        method: step metode
        env: Enviroments, ulike fysiske konstanter.
        verbose: (bool) for om solver skal printe events til konsoll.
    
    output:
        x_num: liste med x-verdier som metoden har løst
        y_num: liste med y-verdier som metoden har løst
        tipped: (bool) for om båten har veltet
        load_overboard: (bool) for om lasten har falt av skipet.

    """

    if(len(y0) != 8): #om init state er mindre enn 8 legger vi til 0 på resten. Dette blir som en failsafe.
        for _ in range(8 - len(y0)):
            y0.append(0)
        
    tipped_print = False #disse endres til True når båten velter eller lasten faller over. SLik at vi kun printer 1 gang.
    load_print = False

    ode = Ode(loadmass, theta_area, env) #initialiserer ode

    #initialiserer
    x_num = np.array([x0])
    y_num = np.array([y0])
    xn = x0
    yn = y0
    
    while xn <= xend:
        if(not ode.tipped):
            xn, yn = method(ode.f, xn, yn, h) #step. Hiver inn dynamikk objektets f (/ode).
        else:
            if(not tipped_print): #båten har veltet, så printer det til konsoll.
                if(verbose): print("boat tipped at", xn, "s") #sørger for at det kun blir printet en gang.
                tipped_print = True

            xn = xn + h #gjør et tidssteg, men state blir konstant.
            yn = np.array([0, 0, yn[2], yn[3], 0, np.sign(yn[5])*np.pi / 2, 0, 0])

        if(abs(yn[7]) >= R and fence): #lasten beveger seg på gjerdet. Setter hastigheten til 0 og posisjonen til +-R
            yn[7] = np.sign(yn[7])*R
            yn[6] = 0
            
        if(abs(yn[7]) >= R and not fence): #Om det ikke er noe gjerde faller lasten av.
            yn[7] = 0
            yn[6] = 0
            ode.loadmass = 0 #setter lastmassen til 0¨

            if(not load_print):
                if(verbose): print("load fell overboard at ", xn)
                load_print = True
    

        y_num = np.concatenate((y_num, np.array([yn]))) #legger til ny state i historikken
        x_num = np.append(x_num,xn)
    
    if(verbose and not tipped_print): #printer viktige events:
        print("Boat did not tip during the simulation.")
    if(verbose and not load_print and ode.loadmass != 0):
        print("Load mass did not fall over board.")

    return  x_num, y_num, tipped_print, load_print # de to siste er booleans for om båten har veltet eller lasten har fallt over bordprint # de to siste er booleans for om båten har veltet eller lasten har fallt over bord

def odeint_solver(x0, xend, y0, loadmass = 0, theta_area=True, fence=False, method = RK45, max_step=0.01, env=ode.Enviroments(), verbose=True):

    if(len(y0) != 8): #om init state er mindre enn 8 legger vi til 0 på resten. Dette blir som en failsafe.
        for _ in range(8 - len(y0)):
            y0.append(0)
        
    tipped_print = False #disse endres til True når våren velter eller lasten faller over. SLik at vi kun printer 1 gang.
    load_print = False

    ode = Ode(loadmass, theta_area, env) #lager ode

    solution = method(ode.f, x0, y0, xend, max_step=max_step)  #solution holder på tid, state og har step funksjonen innebygd.

    #initialiserer
    x_num = np.array([x0])
    y_num = np.array([y0])
    xn = x0
    yn = y0
    
    while(solution.t < xend):
        if(not ode.tipped):
            solution.step() #gjør et step.
        else:
            if(not tipped_print): #båten har veltet, så printer det til konsoll.
                if(verbose): print("boat tipped at", xn, "s") #sørger for at det kun blir printet en gang.
                tipped_print = True
            solution.y = [0, 0, solution.y[2], solution.y[3], 0, np.sign(solution.y[5])*np.pi / 2, 0, 0]

        if(abs(solution.y[7]) >= R and fence): #lasten beveger seg på gjerdet. Setter hastigheten til 0 og posisjonen til +-R
            solution.y[7] = np.sign(solution.y[7])*R
            solution.y[6] = 0
            
        if(abs(solution.y[7]) >= R and not fence): #Om det ikke er noe gjerde faller lasten av.
            solution.y[7] = 0
            solution.y[6] = 0
            ode.loadmass = 0 #setter lastmassen til 0¨

            if(not load_print):
                if(verbose): print("load fell overboard at ", xn)
                load_print = True
    

        y_num = np.concatenate((y_num, np.array([solution.y]))) #legger til ny state i historikken
        x_num = np.append(x_num,solution.t)
    
    if(verbose and not tipped_print): #printer viktige events:
        print("Boat did not tip during the simulation.")
    if(verbose and not load_print and ode.loadmass != 0):
        print("Load mass did not fall over board.")


    return  x_num, y_num, tipped_print, load_print