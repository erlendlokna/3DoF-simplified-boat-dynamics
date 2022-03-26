import numpy as np


#container for holding physical constants:
class Enviroments:
    def __init__(self, k_f=60, mu_last=0.2, F_0 = 0, w_w =0):
        self.k_f = k_f #water friction coef.
        self.mu_last = mu_last #load friction coef.
        self.F_0 = F_0 #wind force amplitude
        self.w_w = w_w #wind force frequency


class Ode:
    def __init__(self, loadmass, theta_area, env):
        self.tipped = False

        self.loadmass = loadmass
        self.theta_area = theta_area
        self.env = env

    def f(self, t, s):
        """
        input:
            t: tid
            s: state, [x, y, vx, vy, omega, theta, vs, s]
        output:
            ds: den deriverte av state.
        """
        ds = np.zeros(len(s))
        
        #pakker ut state vektoren
        x_dot, y_dot = s[0], s[1]
        x, y = s[2], s[3]
        theta_dot = s[4]
        theta = s[5]
        s_Ldot = s[6]
        s_L = s[7]

        beta = get_beta(self.loadmass) #finner beta gitt lastens masse.

        delta_y = y - (R * np.cos(beta / 2) - h) # y - y_C0, tyngde punktets relative bevegelse.

        #sektor vinkel:
        gamma = 2 * np.arccos(np.cos(beta/2) - 4 / (3 * np.pi) * (1 - np.cos(theta)) + delta_y / R)

        #krefter på lasten 
        friction_x =  self.env.mu_last * self.loadmass * abs(g) * np.cos(theta) #friksjon på lasten
        gravity_x =  self.loadmass * abs(g) * np.sin(theta) #tyngdekraft på lasten
            
        F_sx = - np.sign(theta) * abs(friction_x) + np.sign(theta) * abs(gravity_x)
        F_sy = - self.loadmass * abs(g) * np.cos(theta) #y-komponent
        F_Ly = - self.loadmass * abs(g) * (np.cos(theta))**2 #kontakt kraft fra last i x retning
        F_Lx = self.loadmass * abs(g) * np.cos(theta) * np.sin(theta) # -=- i y retning.

        #areal
        if(self.theta_area): A = 0.5 * R**2 * (gamma - np.sin(gamma)) #fortrengt vann
        else: A = 0.5 * R**2 * (beta - np.sin(beta)) #ikke fortrengt vann

        #krefter på båten:
        F_B = A * sigma_0 * abs(g) #boyuansy force. Fancy french
        F_G = - (m + self.loadmass) * abs(g) #gravity
        f_small = - self.env.k_f * R * gamma * theta_dot #friksjons kraft fra vannet på båten
        F_w = self.env.F_0 * np.cos(self.env.w_w * t) #kraft fra vind

        #dreie momenter:
        tau_B = - F_B * h * np.sin(theta) #dreiemoment fra bouyancy kraft
        tau_f = f_small * (y - (R * (np.cos(gamma/2)) - 1))
        tau_w = F_w * y
        tau_L = F_sy * s_L

        #regner ut summen av kreftene på båten: 
        F_sum_x = f_small + F_w + F_Lx #newtons
        F_sum_y = F_G + F_B + F_Ly #newtons

        #regner ut summen av dreiemomenter på båten:
        tau_sum = tau_B + tau_f + tau_w + tau_L

        #lager ny state gitt disse summene:
        ds[0] = F_sum_x / m  #akselerasjon i x retning
        ds[1] = F_sum_y / m  #akselerasjon i y retning
        ds[2] = x_dot #x -> xdot
        ds[3] = y_dot #y -> ydot
        ds[4] = tau_sum / I_c  #theta_dot -> theta_dotdot
        ds[5] = theta_dot #theta -> theta_dot

        if(self.loadmass != 0): ds[6] = F_sx / self.loadmass 
        else: ds[6] = 0 #unngår å dele på 0 om loadmass = 0.

        ds[7] = s_Ldot

        self.tipped = (theta > (np.pi - gamma) / 2) #bool for om båten har tippet.

        return ds