import matplotlib.pyplot as plt

processes = [1, 2, 3, 4, 5, 6, 7, 8]
time_1 = [0.027, 0.033, 0.036, 0.033, 0.041, 0.038, 0.038, 0.042]
time_2 = [0.034, 0.036, 0.036, 0.044, 0.044, 0.046, 0.047, 0.047]
time_3 = [0.059, 0.056, 0.055, 0.067, 0.070, 0.065, 0.067, 0.064]
time_4 = [0.171, 0.170, 0.154, 0.155, 0.143, 0.135, 0.137, 0.134]
time_5 = [0.606, 0.567, 0.465, 0.518, 0.456, 0.441, 0.458, 0.411]
time_6 = [2.350, 1.798, 1.798, 1.729, 1.721, 1.639, 1.622, 1.566]
time_7 = [9.480, 8.582, 7.835, 7.722, 6.817, 6.419, 6.238, 6.303]

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.axis([min(processes), max(processes), 0, max(time_7) * 1.1 ])

plt.plot(processes, time_7, label="24577x18433 (I=7)")
plt.plot(processes, time_6, label="12289x9217 (I=6)")
plt.plot(processes, time_5, label="6145x4609 (I=5)")
plt.plot(processes, time_4, label="3073x2305 (I=4)")
plt.plot(processes, time_3, label="1537x1153 (I=3)")
plt.plot(processes, time_2, label="769x577 (I=2)")
plt.plot(processes, time_1, label="385x289 (I=1)")

handles, labels = ax.get_legend_handles_labels()
# Legend on the right of the figure
# lgd = ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# Legend under the figure
lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
plt.title("Computation time given N processes and I iterations")
plt.xlabel('Number of processes N')
plt.ylabel('Computation time (seconds)')
plt.show()

fig.savefig('benchmarks', bbox_extra_artists=(lgd,), bbox_inches='tight')
