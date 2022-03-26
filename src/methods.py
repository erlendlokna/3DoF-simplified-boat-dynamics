
def euler_step(f, t_i, w_i, h):
    return t_i + h, w_i + h * f(t_i, w_i)

def RK4_step(f, t_i, w_i, h):
    ''' 
    One step using Runge_kutta
    Intput: 
        w_i: [angle theta, frequency ohmega], 
        t_i: time
        h: step length
        f: function to be solved
    Return: next function value
    '''
    k1 = f(t_i,w_i)
    k2 = f(t_i + h/2, w_i + h*k1/2)
    k3 = f(t_i + h/2, w_i + h*k2/2)
    k4 = f(t_i+ h, w_i + h*k3)

    return t_i + h,  w_i + h/6 *(k1+ 2*k2 + 2*k3 + k4)
