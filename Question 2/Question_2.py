import math
import matplotlib.pyplot as plt

# Known values
k1 = 100/60
k2 = 600/60
k3 = 150/60
E0 = 1

# Define functions of ES and S
def f1(time, ES, S):
    function = k1 * E0 * S - (k2 + k3 + k1 * S) * ES
    return function

def f2(time, ES, S):
    function = - k1 * E0 * S + (k1 * S + k2) * ES
    return function

# Input the initial conditions
ES = [0]
S = [10]
h = 0.005 # define a step size
n = int(30/h) # define the time range to be 30s
t = []
# To create the time list
for i in range(n+1):
    t_add = i * h
    t.append(t_add)

# The Runge-Kutta method
def main():
    for i in range(n):
        m_1 = f1(t[i], ES[i], S[i])
        q_1 = f2(t[i], ES[i], S[i])
        m_2 = f1(t[i] + h/2, ES[i] + h/2 * m_1, S[i] + h/2 * q_1)
        q_2 = f2(t[i] + h/2, ES[i] + h/2 * m_1, S[i] + h/2 * q_1)
        m_3 = f1(t[i] + h/2, ES[i] + h/2 * m_2, S[i] + h/2 * q_2)
        q_3 = f2(t[i] + h/2, ES[i] + h/2 * m_2, S[i] + h/2 * q_2)
        m_4 = f1(t[i] + h, ES[i] + h * m_3, S[i] + h * q_3)
        q_4 = f2(t[i] + h, ES[i] + h * m_3, S[i] + h * q_3)

        ES_add = ES[i] + (m_1 + 2 * m_2 + 2 * m_3 + m_4) * h / 6
        S_add = S[i] + (q_1 + 2 * q_2 + 2 * q_3 + q_4) * h / 6
        ES.append(ES_add)
        S.append(S_add)

main()

E = []
P = []
Pv = [] # the velocity of enzyme reaction (rate of change of P)
ES_sum = 0
# To get the value of E, P, and P velocity
for i in ES:
    Ei = E0 - i
    P_velo = k3 * i
    ES_sum += i
    Pi = k3 * ES_sum * h
    E.append(Ei)
    Pv.append(P_velo)
    P.append(Pi)

# Output the result at 30s range
print('The concentration of E: '+ str(E[n]) + ' μM')
print('The concentration of S: '+ str(S[n]) + ' μM')
print('The concentration of ES: '+ str(ES[n]) + ' μM')
print('The concentration of P: '+ str(P[n]) + ' μM')

# Visualization 8.2: changes in concentration of four species
plt.figure()
plt.title('Changes in concentration of four species', family = 'Times New Roman', size = 20)
plt.xlabel('Time(s)', family = 'Times New Roman', size = 15)
plt.ylabel('Concentration(μM)', family = 'Times New Roman', size = 15)
plt.plot(t, E, t, S, t, ES, t, P)
plt.legend(labels=['E(enzyme)', 'S(substrate)', 'ES(intermediate)', 'P(product)'])
plt.show()

# Visualization 8.3: Plot the velocity as a function of the concentration of the substrate
plt.figure()
plt.xlabel('Concentration of Substrate(S)', family = 'Times New Roman', size = 15)
plt.ylabel('Velocity(V) ', family = 'Times New Roman', size = 15)
plt.plot(S, Pv)
plt.show()

# To find exact Vm in the Pv list
print('The maximum value of Velocity: ' + str(max(Pv)))
index_Vm = Pv.index(max(Pv)) # To find the location of Vm in the list
# To find the exact concentration of S at Vm
print('The concentration of S when velocity reaches Vmax: ' + str(S[int(index_Vm)]))