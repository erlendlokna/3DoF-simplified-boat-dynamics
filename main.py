from src.ode import *
from src.graphic import *
from src.solver import *
from src.methods import *
from src.utilities import *


def main():
    #simulating a mass of 0.4 kg on the boat using RK45 method:

    xs_fence2, ys_fence2, _, _ = solver(0, 60 , y0, loadmass= 0.08 * m, theta_area = True, fence=True) #gjør en simulering med gjerde og med en last masse


if __name__ == "__main__":
    main()
