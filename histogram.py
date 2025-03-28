import matplotlib.pyplot as plt
import json
import numpy as np

with open('/home/willow/Python/results/results.json', 'r') as f:
    d = json.load(f)
        
time_in_days = [i/86400 for i in d['all_periods']]
plt.hist(time_in_days, bins=np.logspace(np.log10(1.0), np.log10(10000000.0)), color='lightgreen', edgecolor='green')

plt.yscale('log')
plt.xscale('log')
plt.title('Histogram of Habitable Periods')
plt.xlabel('Time in Days')
plt.ylabel('Number of Habitable Periods')

plt.savefig("/home/willow/Python/Plots/histogram_of_habitable_periods_2nd_attempt.png", dpi=600)

plt.show()
