import matplotlib.pyplot as plt

Fa = 0.01
Fb = 0.99
Ra = 1.1
Rb = 1.0
t = 120

Fat = []
Fbt = []

for i in range(1,t+1):
	Fa = 1/(1+(Fb*(Rb/Ra)**i)/Fa)
	Fb = 1 - Fa
	Fat.append(Fa)
	Fbt.append(Fb)
	
plt.plot(Fbt)
plt.plot(Fat)
plt.show()

