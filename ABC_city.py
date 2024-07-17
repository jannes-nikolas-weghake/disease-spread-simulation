### model with immunity, mobility and death###


# imports

import numpy as np
import matplotlib
import time
import matplotlib.pyplot as plt

matplotlib.rcParams["text.usetex"] = True
plt.rc("font", family="serif")

mobility = True
vaccinations = True
infections = True
deaths = True

# precision = 0.02


def simulateEuler(precision, plot_result=False):
    start = time.time()

    # parameters
    time_period = 240  # in days
    # in days
    sigma = 6  # number of contacts per day
    T = 0.05  # infection probability
    N_cities = 3
    N = [0.250, 0.222, 0.096]  # total number of people in city

    N_calcs = int(time_period / precision)

    # print("Cities: ",N_cities,"| Time span: ",time_period,"| Mobility: ",mobility,"| Vaccinations: ",vaccinations,"| Infections: ",infections,"| Deaths: ",deaths)
    # print("Simulating...")

    alpha = 1 / 16 / precision * np.ones((N_cities, N_calcs))  # recovery rate
    beta = 0.54 / precision * np.ones((N_cities, N_calcs))  # infections rate
    gamma = np.ones((N_cities, N_calcs))  # vaccination rate
    delta = 5 * 1e-6 * np.ones((N_cities, N_calcs))  # death rate

    ### parameter changes
    beta[0] = beta[0] / 4
    beta[2] = beta[2] * 1.3

    gamma[0] = np.zeros(N_calcs)
    gamma[1] = np.ones(N_calcs) * 20000 / (N[1] * 10e6)
    gamma[1][0] = 0
    gamma[1][1] = 0
    gamma[1][2] = 0
    gamma[2] = 1000 / (N[2] * 10e6) * np.ones(N_calcs)
    gamma = gamma / precision

    # beta[1]=beta[1]*1/3
    # gamma[1]=gamma[1]*0.01

    # W = np.random.random((N_cities, N_cities))
    # W = (W+W.T)/2 * 10**-4

    W = np.array([[0, 1.1, 0.2], [1.1, 0, 1.9], [0.2, 1.9, 0]]) * 10e-5

    I = np.zeros((N_cities, N_calcs))  # infected people
    R = np.zeros((N_cities, N_calcs))  # recovered people
    S = np.ones((N_cities, N_calcs))  # susceptible people
    D = np.zeros((N_cities, N_calcs))  # dead people

    for i in range(N_cities):
        for j in range(N_calcs):
            S[i][j] = S[i][j] * N[i]

    ### raw setup done

    I[2][0] = 1e-6  # patient zero

    S[2][0] = S[2][0] - I[2][0]  # correction for patient zero

    ### changing parameters according to turned off features

    if not mobility:
        W = np.zeros((N_cities, N_cities))
    if not vaccinations:
        gamma = np.zeros((N_cities, N_calcs))
    if not infections:
        beta = np.zeros((N_cities, N_calcs))
    if not deaths:
        delta = 0

    #################################################################

    for t in range(0, N_calcs - 1):
        for n in range(0, N_cities):
            I_mob = 0
            S_mob = 0
            R_mob = 0

            for m in range(N_cities):
                if n != m:  # only valid for symmetric W
                    I_mob += W[n][m] * (I[m][t] - I[n][t])
                    S_mob += W[n][m] * (S[m][t] - S[n][t])
                    R_mob += W[n][m] * (R[m][t] - R[n][t])

            infe_people = beta[n][t] * S[n][t] * I[n][t] / N[n]
            reco_people = alpha[n][t] * I[n][t]
            vacc_people = gamma[n][t] * S[n][t]
            dead_people = 0

            if deaths:
                dead_people = delta[n][t] * I[n][t]

            S[n][t + 1] = S[n][t] + precision * (-infe_people - vacc_people) + S_mob

            I[n][t + 1] = (
                I[n][t] + precision * (+infe_people - reco_people - dead_people) + I_mob
            )

            R[n][t + 1] = R[n][t] + precision * (+vacc_people + reco_people) + R_mob

            D[n][t + 1] = D[n][t] + precision * dead_people

    # print(f"finished in {time.time()-start:.6e} sec.")

    if plot_result:
        ### Plotting

        x_time_axis = np.arange(0.0, time_period, precision)
        figure, axis = plt.subplots(N_cities, 1)

        for i in range(N_cities):
            axis[i].plot(x_time_axis, I[i])
            axis[i].plot(x_time_axis, S[i])
            axis[i].plot(x_time_axis, R[i])
            if deaths:
                axis[i].plot(x_time_axis, D[i])

            axis[i].set(ylabel=r"\textbf{people in Million}")
            # axis[i].set(ylim=[-0.01,1.01])

        axis[0].set_title(
            r"\textbf{Disease spreading in cities $A$, $B$ and $C$}", fontsize=16
        )

        if deaths:
            axis[0].legend(
                ["infected", "susceptible", "recovered", "dead"],
                prop={"size": 11},
                loc="upper right",
            )
        else:
            axis[0].legend(
                ["infected", "susceptible", "recovered"], prop={"size": 11}, loc="best"
            )

        axis[0].text(210, 0.1, r"City A", fontsize=14)
        axis[1].text(210, 0.096, r"City B", fontsize=14)
        axis[2].text(210, 0.04, r"City C", fontsize=14)

        axis[N_cities - 1].set(xlabel=r"\textbf{time in days}")

        plt.savefig("figures/ABC_cities_dead.svg", bbox_inches="tight")

        print(f"{D[0][-1]*1e6},{D[1][-1]*1e6},{D[2][-1]*1e6}")
        plt.show(block=False)
        plt.pause(17)
        plt.close()

    return (
        I[0],
        S[0],
        R[0],
        np.arange(0.0, time_period, precision),
        time.time() - start,
    )


simulateEuler(1, True)
