import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import SIR_Euler as Euler
import SIR_Runge_Kutta as Runge_Kutta

matplotlib.rcParams["text.usetex"] = True
plt.rc("font", family="serif")


exactness = 0.001

I_exact = Euler.simulateEuler(exactness)[0]


errorsE = []
errorsRK = []

for i in range(0, 9):
    tempEuler = Euler.simulateEuler(2**-i)[0]
    temprunge = Runge_Kutta.simulateRungeKutta(2**-i)[0]

    temptrueerrorEuler = []
    temptrueerrorRunge = []

    for j in range(np.shape(tempEuler)[0]):
        tempI_exact = I_exact[int(j / exactness * 2**-i)]

        temptrueerrorEuler.append(abs(tempI_exact - tempEuler[j]) / tempI_exact)
        temptrueerrorRunge.append(abs(tempI_exact - temprunge[j]) / tempI_exact)

    errorsE.append(max(temptrueerrorEuler))
    errorsRK.append(max(temptrueerrorRunge))


plt.figure("true_error")

plt.plot(range(0, 9), errorsE, label="Euler")
plt.plot(range(0, 9), errorsRK, label="Runge-Kutta", linestyle="dashed")

plt.title(r'\textbf{Maximal true error to the "exact" solution}', fontsize=16)
plt.xlabel(r"$\textbf{precision in days}$", fontsize=14)
plt.ylabel(r"$\textbf{Maximal relative error in \%}$", fontsize=14)
plt.xticks(
    range(0, 9),
    (r"1", r"$2^{-1}$", r"$2^{-2}$", r"$2^{-3}$", r"$2^{-4}$", r"$2^{-5}$", r"$2^{-6}$", r"$2^{-7}$", r"$2^{-8}$"),
)
plt.legend(prop={"size": 12}, loc="upper right")
axes = plt.gca()
plt.gca().set_aspect(abs((plt.gca().get_xlim()[0] - plt.gca().get_xlim()[1]) / (plt.gca().get_ylim()[0] - plt.gca().get_ylim()[1])) * 1 / 2)
plt.savefig("figures/true_errors.svg", bbox_inches="tight")


plt.show(block=False)
plt.pause(17)
plt.close("all")
