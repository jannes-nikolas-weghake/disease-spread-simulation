### model with immunity, mobility and death###


# imports

import numpy as np
import matplotlib
import time
import matplotlib.pyplot as plt

matplotlib.rcParams["text.usetex"] = True
plt.rc("font", family="serif")

mobility = False
vaccinations = False
infections = True
deaths = False


def simulateEuler(precision):
    start = time.time()

    # parameters
    time_period = 180  # in days
    # in days
    sigma = 6  # number of contacts per day
    T = 0.05  # infection probability

    N = 1  # total number of people in city

    N_calcs = int(time_period / precision)

    alpha = 1 / 25 * np.ones(N_calcs)  # recovery rate
    beta = 0.5 * np.ones(N_calcs)  # infections rate

    # parameter changes
    I = np.zeros(N_calcs)  # infected people
    R = np.zeros(N_calcs)  # recovered people
    S = np.ones(N_calcs)  # susceptible people

    I[0] = 1 * 10**-6

    S[0] = S[0] - I[0]

    #################################################################

    for t in range(0, N_calcs - 1):
        infe_people = beta[t] * S[t] * I[t] / N
        reco_people = alpha[t] * I[t]

        S[t + 1] = S[t] + precision * (-infe_people)

        I[t + 1] = I[t] + precision * (+infe_people - reco_people)

        R[t + 1] = R[t] + precision * (+reco_people)

    return (I, S, R, np.arange(0.0, time_period, precision), time.time() - start)
