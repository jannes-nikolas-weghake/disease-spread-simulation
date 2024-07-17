from SIR_Euler import *
from SIR_Runge_Kutta import simulateRungeKutta
import matplotlib.pyplot as plt

matplotlib.rcParams["text.usetex"] = True
plt.rc("font", family="serif")


IsE = []
IsRK = []
xaxis = []


precision_steps = 9
algo_error_scalingE = 2  # Euler
algo_error_scalingRK = 5  # Runge kutta 4th order
error_scaling_factor = 5


figure, axis = plt.subplots(2, 1)

for i in range(precision_steps):
    IsE.append(simulateEuler(2**-i)[0])
    IsRK.append(simulateRungeKutta(2**-i)[0])

    xaxis.append(np.arange(0.0, 180, 2**-i))
    if i % 2 == 0 and i < 6:
        axis[0].plot(xaxis[i], IsE[i])
        axis[1].plot(xaxis[i], IsRK[i], linestyle="dashed")


legend1 = [r"$h = 1$", r"$h = 1/4$", r"$h = 1/16$", r"$h = 1/64$"]

legend2 = [r"$h = 1$", r"$h = 1/2$", r"$h = 1/4$", r"$h = 1/8$", r"$h = 1/16$"]


axis[0].set_title(r"\textbf{Euler method}", fontsize=16)
axis[1].set_title(r"\textbf{Runge-Kutta method}", fontsize=16)
axis[1].set(xlabel=r"$\textbf{time in days}$")
axis[0].set(ylabel=r"$\textbf{infected population}$")
axis[1].set(ylabel=r"$\textbf{infected population}$")
axis[0].set(xlim=[20, 100])
axis[1].set(xlim=[20, 100])
axis[0].legend(legend1, prop={"size": 12}, loc="upper right")
axis[1].legend(legend1, prop={"size": 12}, loc="upper right")
figure.tight_layout(pad=1.8)


plt.savefig("figures/Euler_RK_infections.svg", bbox_inches="tight")


###numerical error


errorsE = []
errorsRK = []
globalerrorsE = []
globalerrorsRK = []
for i in range(precision_steps - 1):
    tempE = []
    tempRK = []
    for j in range(0, len(IsE[i])):
        tempnumberE = abs(IsE[i][j] - IsE[i + 1][2 * j]) / (2**algo_error_scalingE - 1)
        tempE.append(error_scaling_factor * tempnumberE)
        tempnumberRK = abs(IsRK[i][j] - IsRK[i + 1][2 * j]) / (2**algo_error_scalingRK - 1)
        tempRK.append(error_scaling_factor * tempnumberRK)

    errorsE.append(tempE)
    errorsRK.append(tempRK)

for i in range(precision_steps - 1):
    tempE = 0
    tempRK = 0
    for j in range(0, len(IsE[i])):
        tempE += errorsE[i][j]
        tempRK += errorsRK[i][j]

    globalerrorsE.append(tempE * (2**-i) ** 2)
    globalerrorsRK.append(tempRK * (2**-i) ** 2)


plt.figure("global_errors")
plt.axhline(y=globalerrorsE[-1], color="g", linestyle="-")
plt.scatter(range(precision_steps - 1), globalerrorsE, marker="x", label=r"Euler method")
plt.axhline(y=globalerrorsRK[-1], color="g", linestyle="-")
plt.scatter(range(precision_steps - 1), globalerrorsRK, marker=".", label=r"Runge-Kutta method")


# plt.title(r"\textbf{Global truncation errors scaling with precision}", fontsize=16)
plt.xlabel(r"$\textbf{precision in days}$", fontsize=16)
plt.ylabel(r"$\textbf{global truncation error}$", fontsize=16)
plt.text(5, 0.4, r"$y\approx 3.4e-05$", fontsize=14)
plt.legend(prop={"size": 13}, loc="upper right")
plt.xticks(
    range(precision_steps - 1),
    (
        r"1",
        r"$2^{-1}$",
        r"$2^{-2}$",
        r"$2^{-3}$",
        r"$2^{-4}$",
        r"$2^{-5}$",
        r"$2^{-6}$",
        r"$2^{-7}$",
    ),
)


plt.gca().set_aspect(abs((plt.gca().get_xlim()[0] - plt.gca().get_xlim()[1]) / (plt.gca().get_ylim()[0] - plt.gca().get_ylim()[1])) * 1 / 2)

plt.savefig("figures/global_errors.svg", bbox_inches="tight")
# plt.close("global_errors")
plt.close("all")


figure, axis = plt.subplots(2, 1)


for i in range(0, precision_steps - 1):
    temp = np.arange(0, len(IsE[0]), 2 ** (-i))
    axis[0].plot(temp, errorsE[i])
    axis[1].plot(temp, errorsRK[i], linestyle="dashed")


axis[0].set_title(r"\textbf{local truncation errors}", fontsize=16)
axis[1].set(xlabel=r"$\textbf{time in days}$")
axis[0].set(ylabel=r"$\textbf{local truncation error}$")
axis[1].set(ylabel=r"$\textbf{local truncation error}$")
axis[0].legend(legend2, prop={"size": 12}, loc="upper right")
axis[1].legend(legend2, prop={"size": 12}, loc="upper right")
axis[0].set(ylim=[0, 0.44])
axis[1].set(ylim=[0, 0.44])
axis[0].set(xlim=[10, 80])
axis[1].set(xlim=[10, 80])
axis[0].text(12, 0.3, r"Euler", fontsize=14)
axis[1].text(12, 0.3, r"Runge-Kutta", fontsize=14)
plt.savefig("figures/local_num_error_euler.svg", bbox_inches="tight")

# plt.close("num_err")


plt.figure("inf+num_err")

plt.plot(xaxis[0], IsE[0], label=r"h = 1 day")
plt.fill_between(xaxis[0], IsE[0] - errorsE[0], IsE[0] + errorsE[0], alpha=0.5, label="error band")


plt.plot(xaxis[4], IsE[4], label=r"h = 1/16 day", linestyle="dotted")
plt.fill_between(xaxis[4], IsE[4] - errorsE[4], IsE[4] + errorsE[4], alpha=0.5, label="error band")


plt.title(r"infections with error bands (local truncation error scaled by 5)", fontsize=12)
plt.xlabel(r"$\textbf{time in days}$", fontsize=16)
plt.ylabel(r"$\textbf{people in millions}$", fontsize=16)
plt.legend(prop={"size": 12}, loc="upper right")

plt.gca().set_aspect(abs((plt.gca().get_xlim()[0] - plt.gca().get_xlim()[1]) / (plt.gca().get_ylim()[0] - plt.gca().get_ylim()[1])) * 1 / 2)

plt.savefig("figures/numerical_errorbands.svg", bbox_inches="tight")

plt.close("inf+num_err")


plt.show(block=False)
plt.pause(17)
plt.close("all")
