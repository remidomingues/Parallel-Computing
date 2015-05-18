import matplotlib.pyplot as plt

processes = [1, 2, 3, 4, 5, 6, 7, 8]

# With gathering data
time_10 = [2143, 1573, 1253, 1126, 1053, 936, 760, 606]
time_11 = [8710, 6350, 5090, 4390, 3840, 3683, 3493, 3206]
time_12 = [35100, 25256, 20336, 18370, 17573, 15146, 13183, 12990]

# Without gathering data
# time_10 = [1880, 1130, 883, 720, 626, 513, 316, 326]
# time_11 = [7663, 4606, 3450, 2740, 2413, 2036, 1256, 1213]
# time_12 = [32186, 18773, 14030, 10820, 9473, 8136, 5036, 4746]

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.axis([min(processes), max(processes), 0, max(time_12) * 1.05 ])

plt.plot(processes, time_12, label="28673x28673 (I=12)")
plt.plot(processes, time_11, label="14337x14337 (I=11)")
plt.plot(processes, time_10, label="7169x7169 (I=10)")

handles, labels = ax.get_legend_handles_labels()
# Legend on the right of the figure
# lgd = ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# Legend under the figure
lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1))
plt.title("Computation time given N processes and I iterations")
plt.xlabel('Number of processes N')
plt.ylabel('Computation time (ms)')
plt.show()

fig.savefig('benchmarks', bbox_extra_artists=(lgd,), bbox_inches='tight')
