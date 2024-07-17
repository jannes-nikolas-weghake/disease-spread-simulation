import numpy as np
import matplotlib
import time
import matplotlib.pyplot as plt

import SIR_Euler as Euler
import SIR_Runge_Kutta as Runge_Kutta

time_euler = []
time_runge_kutta = []

for i in range(1000):
    time_euler.append(      Euler.simulateEuler(2**-8)[-1])
    time_runge_kutta.append(Runge_Kutta.simulateRungeKutta(2**-8)[-1])
    if i%10==0:
        print(i)

time_euler =        np.array(time_euler)
time_runge_kutta =  np.array(time_runge_kutta)

print(np.mean(time_euler),'+-',np.std(time_euler))
print(np.mean(time_runge_kutta),'+-',np.std(time_runge_kutta))



plt.figure("euler_vs_rungekutta")

temp1 = Euler.simulateEuler(1/32)
temp2 = Runge_Kutta.simulateRungeKutta(1/32)
plt.plot(temp1[-2], temp1[0])
plt.plot(temp2[-2], temp2[0])


