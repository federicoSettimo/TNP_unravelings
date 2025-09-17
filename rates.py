import matplotlib.pyplot as plt
import numpy as np

filein = open("tmax.txt")
tmax = float(filein.readline())

def read_file_columns(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            numbers = line.split()
            for i, number in enumerate(numbers):
                if len(data) <= i:
                    data.append([])
                data[i].append(float(number))
    return data

rates = read_file_columns("rates_H.txt")

t = np.linspace(0, tmax, len(rates[0]))

fig, ax = plt.subplots(1,2, figsize=(26,5), sharex=True, sharey=False)
fontsize = 18

for i in range(3):
    ax[0].plot(t, rates[i])
ax[0].axhline(0, color='black', linewidth=1.)

ax[0].set_title(r"Eigenvalues of the Choi matrix", fontsize=fontsize-2)

ax[1].plot(t, rates[4])
ax[1].axhline(1, color='black', linewidth=1.)
ax[1].set_title(r'Maximum of the Bloch vector', fontsize=fontsize-2)

for i in range(2):
    ax[i].set_xlabel(r"$t$", fontsize=fontsize-3)
    ax[i].tick_params(axis='both', which='major', labelsize=fontsize-3)

#plt.suptitle(r"P and CP conditions for $\Lambda_{t+dt,t}$", fontsize=fontsize)

plt.savefig("Choi_Schro.png", dpi=300)
plt.show()