### model with immunity, mobility h###

# no vaccination, no mobility


import numpy as np
import matplotlib
import time
import matplotlib.pyplot as plt

matplotlib.rcParams["text.usetex"] = True
plt.rc("font", family="serif")

vaccinations = False
infections = True


# precision = 0.02


def simulateRungeKutta(precision):
    start = time.time()

    # parameters
    time_period = 180  # in days
    # in days
    sigma = 6  # number of contacts per day
    T = 0.05  # infection probability

    N = 1  # total number of people in city

    N_calcs = int(time_period / precision)

    gamma = 1 / 25 * np.ones(N_calcs)  # recovery rate, once alpha
    beta = 0.5 * np.ones(N_calcs)  # infections rate

    alpha = 10 ** (-4) * np.ones(N_calcs)

    I = np.zeros(N_calcs)  # infected people
    R = np.zeros(N_calcs)  # recovered people
    S = np.ones(N_calcs)  # susceptible people

    # raw setup done

    I[0] = 1 * 10**-6  # patient zero

    # correction for patient zero
    S[0] = S[0] - I[0]

    #################################################################

    for t in range(0, N_calcs - 1):
        # calc of k for susceptible
        k_S_1 = -beta[t] * S[t] * I[t]
        k_I_1 = beta[t] * S[t] * I[t] - gamma[t] * I[t]

        k_S_2 = (
            -beta[t] * (S[t] + k_S_1 * precision / 2) * (I[t] + k_I_1 * precision / 2)
        )
        k_I_2 = beta[t] * (S[t] + k_S_1 * precision / 2) * (
            I[t] + k_I_1 * precision / 2
        ) - gamma[t] * (I[t] + k_I_1 * precision / 2)

        k_S_3 = (
            -beta[t] * (S[t] + k_S_2 * precision / 2) * (I[t] + k_I_2 * precision / 2)
        )
        k_I_3 = beta[t] * (S[t] + k_S_2 * precision / 2) * (
            I[t] + k_I_2 * precision / 2
        ) - gamma[t] * (I[t] + k_I_2 * precision / 2)

        k_S_4 = -beta[t] * (S[t] + k_S_3 * precision) * (I[t] + k_S_3 * precision)
        k_I_4 = beta[t] * (S[t] + k_S_3 * precision) * (
            I[t] + k_S_3 * precision
        ) - gamma[t] * (I[t] + k_I_3 * precision)

        S[t + 1] = S[t] + precision / 6 * (k_S_1 + 2 * k_S_2 + 2 * k_S_3 + k_S_4)

        I[t + 1] = I[t] + precision / 6 * (k_I_1 + 2 * k_I_2 + 2 * k_I_3 + k_I_4)

        R[t + 1] = N - S[t + 1] - I[t + 1]

    return (I, S, R, np.arange(0.0, time_period, precision), time.time() - start)
