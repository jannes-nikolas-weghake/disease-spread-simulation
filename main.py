import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import SIR_Euler as Euler
import SIR_Runge_Kutta as Runge_Kutta


euler = Euler.simulateEuler(0.01)
runge_kutta = Runge_Kutta.simulateRungeKutta(0.01)

plt.figure
plt.plot(euler[3], euler[0])
plt.plot(runge_kutta[3], runge_kutta[0])

plt.show()
