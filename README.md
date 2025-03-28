# three_body_problem

## Overview
This is code for a three-star system (with four bodies including the planet) that I used in my eighth grade science project. I used much of Patrick Youssef's code for the physics and plotting of the simulation itself (https://patrickyoussef.com/blog/nbody/), while referencing my dad for most of the other features that I implemented.

## How to use
`n_body_simulation.py` runs the simulation. Under the NBody class, you can change the distance of two stars from the first star in AU and the velocity in m/s when defining __init__. You can also change each star's luminosity, the albedo of the planet, and how close two bodies have to be in AU in order for a collision to be detected.

The last lines of code save the results of the simulation to csv files. Using `postprocess.py`, you can open the csv files and generate results (longest habitable period of each simulation, longest of all simulations, etc.) as a json file.

`histogram.py` takes the json you just retrieved and creates a histogram from it.

`temperaturegram.py` takes a specific csv file and creates a plot of the change in temperature.
