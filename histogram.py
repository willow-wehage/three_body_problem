import matplotlib.pyplot as plt
import json

with open('/home/willow/Python/results/results.json', 'r') as f:
    d = json.load(f)
        
time_in_days = [i/86400 for i in d['all_periods']]
plt.hist(time_in_days, bins=50, color='lightgreen', edgecolor='green')

plt.yscale('log')
plt.title('Histogram of Habitable Periods')
plt.xlabel('Time in Days')
plt.ylabel('Count')

plt.show()
