import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import cmath


#################################

# DEBUT DONNEES
L_V = np.array([25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1100, 1158])*10**(-3)
L_V_inc = np.ones(len(L_V))*10**(-3)*2

L_P = np.array([0, 0.01, 0.01, 0.02, 0.02, 0.03, 0.04, 0.05, 0.07, 0.12, 0.76, 6.44, 16.15, 25.14, 34, 42.7, 51.6, 60.5, 69, 74])*10**(-3)
L_P_inc = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])*10**(-3)
# FIN DONNEES

def modele(x, a, b):
	return a*x + b

opt, pcov = curve_fit(modele, L_V[10:], L_P[10:], sigma=L_V_inc[10:], absolute_sigma=True)
print("a = ", round(opt[0], 3), " +- ", round(np.sqrt(np.diag(pcov))[0], 3), " W*V**(-1)")
print("b = ", round(opt[1], 3), " +- ", round(np.sqrt(np.diag(pcov))[1], 3), " W")
L_V_fit = np.linspace(min(L_V[10:]), max(L_V[10:]), 1000)
Fact_Barres = 10


fig, ax = plt.subplots()
ax.errorbar(L_V[:10], L_P[:10], xerr=L_V_inc[:10]*Fact_Barres, yerr=L_P_inc[:10]*Fact_Barres, elinewidth=1, ecolor="red", fmt='None')
ax.scatter(L_V[:10], L_P[:10], label="Mode spontanné", s=10, c="red")

ax.errorbar(L_V[10:], L_P[10:], xerr=L_V_inc[10:]*Fact_Barres, yerr=L_P_inc[10:]*Fact_Barres, elinewidth=1, ecolor="blue", fmt='None')
ax.scatter(L_V[10:], L_P[10:], label="Mode laser", s=10, c="blue")

ax.plot(L_V_fit, opt[0]*L_V_fit + opt[1], label="fit : P = {}V + {}".format(round(opt[0], 3), round(opt[1], 3)))
ax.axvline(330*10**(-3), c="green", label="Séparations V=0.33 V", linestyle="dashdot", linewidth=1)
ax.set_xlabel("Tension (V)")
ax.set_ylabel("Puissance (W)")
plt.title("Les régimes de la diode laser")
plt.grid()
plt.legend()




#On connait la puissance de la lumière bleue par rapport au deltaV
#################################

# DEBUT DONNEES
L_V = np.array([150, 200, 250, 315, 320, 340, 350, 400, 500, 600, 680, 700, 900])*10**(-3)
L_V_inc = np.ones(len(L_V))*10**(-3)*2

L_Pmax = np.array([0.02, 0.03, 0.04, 0.09, 0.1, 0.14, 0.5, 5.05, 12.69, 19.76, 17.85, 25.56, 37.7])*10**(-3)
L_Pmin = np.array([0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.07, 0.13, 0.19, 0.19, 0.26, 0.42])*10**(-3)

L_Pmin_inc = np.ones(len(L_Pmax))*0.01*10**(-3)*0
L_Pmax_inc = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.04, 0.04, 0.04, 0.04, 0.1])*10**(-3)
# FIN DONNEES

L_R = (L_Pmax-L_Pmin)/(L_Pmax+L_Pmin)
print(L_R)
#Propagation incertitudes
delta_d = L_Pmin_inc + L_Pmax_inc
L_R_inc = np.abs(L_R*np.sqrt( ((delta_d)/(L_Pmax-L_Pmin))**2 + ((delta_d)/(L_Pmax+L_Pmin))**2 ))

k = 5
Fact_Barres = 0.5
fig, ax = plt.subplots()
plt.plot(L_V, L_R, linewidth=1, color="black")
ax.errorbar(L_V[:k], L_R[:k], xerr=L_V_inc[:k]*Fact_Barres, yerr=L_R_inc[:k]*Fact_Barres, elinewidth=1, ecolor="red", fmt='None')
ax.scatter(L_V[:k], L_R[:k], label="Mode spontanné", s=10, c="red")

ax.errorbar(L_V[k:], L_R[k:], xerr=L_V_inc[k:]*Fact_Barres, yerr=L_R_inc[k:]*Fact_Barres, elinewidth=1, ecolor="blue", fmt='None')
ax.scatter(L_V[k:], L_R[k:], label="Mode laser", s=10, c="blue")
ax.axvline(330*10**(-3), c="green", label="Séparations V=0.33 V", linestyle="dashdot", linewidth=1)
ax.set_xlabel("Tension (V)")
ax.set_ylabel("Facteur R")
plt.title("Polarisation du laser")
plt.grid()
plt.legend()






########################################################################
# DEBUT DONNEES
Cb = (4.10*10**(-6))/(837*10**(-3))

Cb_inc = Cb*np.sqrt(((0.02*10**(-6))/(4.10*10**(-6)))**2 + ((10*10**(-3))/(837*10**(-3)))**2)
print(Cb, Cb_inc)





L_P = np.array([3.75, 3.66, 3.02, 2.27, 1.420, 0.841, 0.3, 0.1, 0.029, 0.027, 0.039, 0.099, 0.329])*10**(-6)
L_P_inc = np.array([0.01, 0.01, 0.01, 0.01, 0.002, 0.002, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001])*10**(-6)

L_Delta_V = np.array([812, 775, 650, 495, 300, 170, 62.6, 25, 0.5, 0.0001, 0.5, 20, 62.7])*10**(-3)
L_Delta_V_inc = np.ones(len(L_Delta_V))*10*10**(-3)

L_P_Delta = Cb*L_Delta_V
L_P_Delta_inc = L_P_Delta*np.sqrt(((Cb_inc)/(Cb))**2 + ((L_Delta_V_inc)/(L_Delta_V))**2)


L_alpha = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])*np.pi/180
L_alpha_inc = np.ones(len(L_alpha))*np.pi/180


P0 = 160*10**(-3)
P0_inc = 1*10**(-3)

# FIN DONNEES

L_P = L_P[:-3]
L_P_inc = L_P_inc[:-3]
L_Delta_V = L_Delta_V[:-3]
L_alpha = L_alpha[:-3]
L_alpha_inc = L_alpha_inc[:-3]
L_P_Delta = L_P_Delta[:-3]
L_P_Delta_inc = L_P_Delta_inc[:-3]

L_P_IR = P0*np.cos(2*L_alpha)**2
print("B  ", np.sqrt(L_P_Delta))
L_P_IR_inc = np.sqrt((np.cos(2*L_alpha)**4)*(P0_inc**2) + (L_alpha_inc**2)*(4*P0*np.cos(2*L_alpha)*np.sin(2*L_alpha))**2)





def modele(x, a, b):
	return a*x + b

opt, pcov = curve_fit(modele, L_P_IR, np.sqrt(L_P), sigma=L_P_IR_inc, absolute_sigma=True)
a, b = round(opt[0], 5), round(opt[1], 5)
print(a, b)
print(np.sqrt(np.diag(pcov))[0])


opt2, pcov2 = curve_fit(modele, L_P_IR, np.sqrt(L_P_Delta), sigma=L_P_IR_inc, absolute_sigma=True)
a2, b2 = round(opt2[0], 5), round(opt2[1], 5)
print(a2, b2)
print(np.sqrt(np.diag(pcov2))[0])


mu0 = 4*np.pi*10**(-7)
epsi = 8.85*10**(-12)
w = 2.2*10**(15)
n = 2.3
l = 10**(-3)
S = 5*10**(-9)


print("AAAA1111", a*np.sqrt(0.5*((S*n**3)/((w*l)**2)) * (epsi/mu0)**(3/2)), np.sqrt(np.diag(pcov))[0]*np.sqrt(0.5*((S*n**3)/((w*l)**2)) * (epsi/mu0)**(3/2)))
print("AAA2222  ", a2*np.sqrt(0.5*((S*n**3)/((w*l)**2)) * (epsi/mu0)**(3/2)), np.sqrt(np.diag(pcov2))[0]*np.sqrt(0.5*((S*n**3)/((w*l)**2)) * (epsi/mu0)**(3/2)))
Fact_Barres = 1

fig, ax = plt.subplots()

ax.errorbar(L_P_IR, np.sqrt(L_P), xerr=L_P_IR_inc*Fact_Barres, yerr=L_P_inc*Fact_Barres, elinewidth=1, ecolor="blue", fmt='None')
ax.scatter(L_P_IR, np.sqrt(L_P), label="Mesure Puissance-mètre", s=10, c="blue")
ax.plot(np.linspace(min(L_P_IR), max(L_P_IR), 1000), a*np.linspace(min(L_P_IR), max(L_P_IR), 1000) + b, linestyle="--", linewidth=1, color="blue", label="fit : sqrt(P(2w)) = {}*P(w) + {}".format(a, b))

ax.errorbar(L_P_IR, np.sqrt(L_P_Delta), xerr=L_P_IR_inc*Fact_Barres*0.75, yerr=L_P_Delta_inc*Fact_Barres, elinewidth=1, ecolor="red", fmt='None')
ax.scatter(L_P_IR, np.sqrt(L_P_Delta), label="Mesure Oscilloscope", s=10, c="red")
ax.plot(np.linspace(min(L_P_IR), max(L_P_IR), 1000), a2*np.linspace(min(L_P_IR), max(L_P_IR), 1000) + b2, linestyle="--", linewidth=1, color="red", label="fit : sqrt(P(2w)) = {}*P(w) + {}".format(a2, b2))

ax.set_xlabel("P(w) en W")
ax.set_ylabel("sqrt(P(2w)) en W**1/2")
plt.title("Coefficient non linéaire")
plt.grid()
plt.legend()


fig, ax = plt.subplots()
ax.errorbar(L_alpha, L_P, xerr=L_alpha_inc*Fact_Barres, yerr=L_P_inc*Fact_Barres, elinewidth=1, ecolor="blue", fmt='None')
ax.scatter(L_alpha, L_P, label="Mesure Puissance-mètre", s=10, c="blue")

ax.set_xlabel("alpha-alpha0 en rad")
ax.set_ylabel("P(2w) en W")
plt.title("Coefficient non linéaire")
plt.grid()
plt.legend()






plt.show()
