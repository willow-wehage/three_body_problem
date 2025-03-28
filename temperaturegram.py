import matplotlib.pyplot as plt
import numpy as np

file_path = '/home/willow/Python/three_star_csv_files/system_states_20250224_T18_18_51_2.csv'

data = np.genfromtxt(
        file_path, delimiter=',', names=True, dtype=None
)

plt.plot(data["time_s"]/86400, data["T_K"], color='orange')
#plt.ylim(0, 325)

#plt.xscale('log')
plt.title('Time vs Temperature')
plt.xlabel('Time in Days')
plt.ylabel('Temperature in Kelvins')

plt.savefig("/home/willow/Python/Plots/temperature_plot_zoomed_out.png", dpi=600)
plt.show()