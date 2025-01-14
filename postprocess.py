import os
import numpy as np
import json

folder_path = "/home/willow/Python/three_star_csv_files"

max_habitable_period = 0
file_with_max = None
i = 0
d = {}
d["results"] = []
d["all_periods"] = []

for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
     
    # Check if it's a file (not a directory)
    if os.path.isfile(file_path):
        
        habitable_periods = []
        habitable_bool = False
        start_time = 0
        end_time = 0
        
        data = np.genfromtxt(
            file_path, delimiter=',', names=True, dtype=None
        )
        for t, temp in zip(data["time_s"], data["T_K"]):
            if temp >= 273.15 and temp <= 310.928:
                if habitable_bool is True:
                    habitable_bool = True
                else: 
                    habitable_bool = True
                    start_time = t
            else:
                if habitable_bool is False:
                    habitable_bool = False
                else:
                    habitable_bool = False
                    end_time = t
                    habitable_periods.append(end_time - start_time)
                    d["all_periods"].append(end_time - start_time) 
        
        if habitable_periods:
            
            mean_time = np.mean(habitable_periods)
            max_time = np.max(habitable_periods)
            if max_time > max_habitable_period:
                max_habitable_period = max_time
                file_with_max = filename
        else:
            max_time = 0
            mean_time = 0
            
        d["results"].append(
            {
                "filename": filename, 
                "habitable_periods": len(habitable_periods),
                "mean": mean_time,
                "max": max_time
            }
        )
        print(f"[{i}] In file {filename} there were {len(habitable_periods)} habitable periods. Mean: {mean_time} seconds, Max: {max_time} seconds.")
        
        i+=1
        
        '''for i, period in enumerate(habitable_periods):
            print(f"    {i}: {period} s")'''
            
print (f"The longest period lasted for {max_habitable_period} seconds in file {file_with_max}.")

d["max_habitable_period"] = max_habitable_period
d["file_with_max_habitable_period"] = file_with_max

with open("/home/willow/Python/results/results.json", "w") as f:
    json.dump(d, f)