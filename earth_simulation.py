import copy # Used for creating non-linked copies of states
from math import ceil # Helpful for rounding off floats to conservative ceilings
from itertools import combinations # Produces the sets of combinations for pairs
import numpy as np # type: ignore # General matrix compute platform
import matplotlib.pyplot as plt # type: ignore # Used for plotting and analysis
import time 
import os
import datetime


def flux(luminosity, d):
        return luminosity / (4.0 * np.pi * d**2)
    
def rand_unit_vec():
    vec = np.random.normal(size=3)
    return vec / np.linalg.norm(vec)

sigma = 5.6866e-8 # Wm-2K-4

class NBody:

    def __init__(self, Xi=None, masses=None, G=6.67408e-11, distance_au=20.0, velocity_mps=1000.0, luminosity=[3.846e+26], albedo=.3, collision_tol_au=.01):

        # Init State and Masses
        # Set the initial state matrix
        if Xi is None:
            self.Xi = self.get_init_conditions_from_dist_au(distance_au, velocity_mps)
        else:
            self.Xi = Xi 

        self.masses = masses # Set the body masses array
        # Set gravitational constant
        # This really shouldn't change but it's fun for some special cases
        self.G = G
        self.history = None
        self.energies = None
        self.temperature = None
        self.luminosity = luminosity
        self.albedo = albedo
        self.collision_tol_m = collision_tol_au * 149597870700

    def get_energy(self, X):
        # Useful Variables
        N, D = X.shape # Get the number of bodies, and dimensionality
        D = D // 2 # Get the number of dimensions, are we 2d or 3d?
        R = X[:, :D] # Submatrix with all positions
        V = X[:, D:] # Submatrix with all velocities

        # Determine Kinetic Energy
        # 1/2 * mass * v ^ 2
        KE = 0
        for i in range(N):
            KE += 0.5 * self.masses[i] * np.linalg.norm(V[i]) ** 2

        # Determine Potential Energy
        # (-G * m_i * m_j) / r_ij
        PE = 0
        for body_i, body_j in self.pairs:
            r = np.linalg.norm(R[body_j] - R[body_i]) # Distance between bodies
            PE -= self.masses[body_i] * self.masses[body_j] / r
        PE *= self.G # Multiplying is expensive so I only do one at the end

        return KE, PE

    def get_state_deriv(self, X):
        """
        Inputs:
        - X: The current state to determine the state derivative for
        Returns:
         - Xdot: The corresponding state derivative for the input state
        """
        print(f"X.dtype: {X.dtype}")
         # Useful Variables
        N, D = X.shape # Get the number of bodies, and dimensionality
        D = D // 2 # Get the number of dimensions, are we 2d or 3d?
        R = X[:, :D] # Submatrix with all positions
        V = X[:, D:] # Submatrix with all velocities
        print ("Number of bodies ", N)
        print ("D: ", D)
        # Build Placeholder Structure
        Xdot = np.zeros_like(X) # Xdot is the same size as X
        
        print(f"Xdot.dtype: {Xdot.dtype}")

        Xdot[:, :D] = V # Fill in velocities from state 

        # Iterate Over Pairs and Fill Out Acceleration
        # self.pairs gets defined when we start a sim
        # body_i, body_j are the indices of the bodies
        for body_i, body_j in self.pairs:

            # Get vector from body_i => body_j and its magnitude
            r1, r2 = R[body_i], R[body_j] # Positions of body_i and body_j
            r_vec = r2 - r1 # Vector from body_i => body_j
            r = np.linalg.norm(r_vec) # Distance from body_i => body_j

            # Find Force from body_i => body_j
            F = self.G * self.masses[body_i] * self.masses[body_j] * r_vec / r**3
            a1 =  F / self.masses[body_i] # Compute acceleration for body_i
            a2 = -F / self.masses[body_j] # Compute acceleration for body_j

            print(f"a1: {a1} units")
            
            print(Xdot[body_i, D:])
            
            # Apply acceleration to body_i and body_j
            Xdot[body_i, D:] += a1
            Xdot[body_j, D:] += a2

        return Xdot

    def get_temperature(self, X):
        sum = 0
        d = []
        for i in range (1):
            d.append(np.sqrt((X[1][0] - X[i][0])**2 + (X[1][1] - X[i][1])**2 + (X[1][2] - X[i][2])**2))
        
        for d_i, l_i in zip(d, self.luminosity):
            sum += flux(l_i, d_i) 
        numerator = sum * (1.0 - self.albedo)
        denominator = 4.0 * sigma
        return (numerator / denominator)**0.25

    def detect_collision(self, X):
        
        for i in range (2):
            for j in range (i+1, 2):
                d = np.sqrt((X[i][0] - X[j][0])**2 + (X[i][1] - X[j][1])**2 + (X[i][2] - X[j][2])**2)
                if d < self.collision_tol_m:
                    return True
        return False
    
    def get_init_conditions_from_dist_au(self, distance_au, velocity):
        
        init_conditions = []
        sun_i = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        init_conditions.append(sun_i) 
        
        distance_m = distance_au * 149597870700.0      

        planet_i = [149597870700.0, 0.0, 0.0, 0.0, 29750.0, 0.0]
        init_conditions.append(planet_i)
        
        return (np.array(init_conditions))

    def rk4(self, X, dt, evaluate):
        """
        Inputs:
        - X: Current state of system
        - dt: Integration Timestep
        - evaluate: Function that will return the derivative for the state
        Returns:
        - X: Updated state one timestep later
        """
        # Calculate Terms
        k1 = evaluate(X)
        k2 = evaluate(X + 0.5*k1*dt)
        k3 = evaluate(X + 0.5*k2*dt)
        k4 = evaluate(X + k3*dt)

        # Update X
        X_prime = (1/6.)*(k1 + 2*k2 + 2*k3 + k4)
        return X + X_prime * dt


    def run_simulation(self, T, dt):
        """
        Inputs:
        - T: Total runtime of simulation
        - dt: Timestep for integration
        Returns:
        - history: Matrix of the history of states
        """
    
        # Check to ensure initial conditions and masses have been set
        assert self.Xi is not None
        assert self.masses is not None
    
    
        if self.detect_collision(self.Xi) is True:
            raise Exception("Error, bodies started in collision state. Refusing to run simulation.")

        # Setup Sim Params
        iters = ceil(T / dt) # Number of simulation iterations
    
        # Init History
        N, D = self.Xi.shape
        print(N, D, iters)
        self.history = np.zeros((iters+1, N, D))
        self.history[0] = self.Xi # First history is our initial conditions

        # Determine Force Pair Indexes
        self.pairs = list(combinations(range(N), 2))
        print ("self pairs", self.pairs)

        # Init Energies
        self.energies = np.zeros((iters+1, 3))
        KE, PE = self.get_energy(self.Xi)
        self.energies[0] = np.array([KE, PE, KE+PE])
        
        # Init Temperatures
        self.temperature = np.zeros((iters+1, 1))
        self.temperature[0] = np.array(self.get_temperature(self.Xi))
    
        # Run Simulation Iterations
        X = copy.deepcopy(self.Xi) # Copy as to not modify Xi
        print("starting simulation...")
        start = time.time()
        self.times = np.zeros((iters+1, 1))
        self.times[0] = 0.0
        for i in range(iters):
            X = self.rk4(X, dt, self.get_state_deriv) # Get new state
            self.history[i+1] = X # Store new state
            KE, PE = self.get_energy(X) # Get new state's energy 
            T = self.get_temperature(X)
            self.energies[i+1] = np.array([KE, PE, KE+PE]) # Store energy
            self.temperature[i+1] = np.array([T])
            self.times[i+1] = self.times[i] + dt
            if self.detect_collision(X) is True:
                raise Exception("Error, bodies in collision. Stopping simulation.")

            # print status every 1000 iterations
            if not i % 10000:
                print(f"\r{int(time.time() - start)}s elapsed, iteration {i}/{iters}, {int(i/iters * 100)}% complete")
                
        return self.history
    
'''
def earth_stable_orbit(r):
        """
        Inputs:
        - r: Altitude of orbit above earth's surface in meters
        Returns:
        - v: Velocity in m/s to sustain a stable circular orbit
        """
        G = 6.67408e-11
        massE = 5.974e24 # Mass in kg
        rE = 6.3781e6 # Radius in m
        return np.sqrt(G * massE / (r + rE))
'''

# Setting up initial state
names = ["Sun1", "Trisolaris"]
#X = np.array([
    #[863703709188.958, 0, 0, -1.12474E+1, 7.54876E+0, 0],
    #[-431851854594.4791, 7.48e+11, 0, 1.16497E+4, -4.14793E+4, 0],
    #[-431851854594.4791, -7.48e+11, 0, -3.22930E+04, 1.36960E+04, 0],
    #[-1.43778E+11, -4.00067E+10, 0, 7.65151E+03, -2.87514E+04, 0],
#]) 
masses = np.array([1.98854E+30, 5.97219E+24])

# You can subset out the planets
#n = 0 # Number of planets to remove from the end
# names = names[:-n]
#X = X[:-n]
#masses = masses[0:-n]

years = 1
timestep_days = 1/24

for p in range (1):
    SolarSystem = NBody(Xi=None, masses=masses)

    T, dt = years * 365 * 24 * 60**2, timestep_days * 24 * 60 ** 2 
   
    history = SolarSystem.run_simulation(T, dt)
    
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111,projection='3d')
    colors = ["yellow", "blue"]
    for i in range(history.shape[1]):
        c = colors[i]
        x = history[:, i, 0]
        y = history[:, i, 1]
        z = history[:, i, 2]
        if i < 3:
            ms = 20
        else:
            ms = 10
        ax.plot3D(x, y, z, color='gray', label=names[i], linewidth=0.2,
                 markevery=[0], marker='o', ms=ms, mfc=c, mec="black", mew=0.5)
        
    ax.set(xlabel='X')
    ax.set(ylabel='Y')
    ax.set(zlabel='Z')

    ax.grid(False) # Turn off grid
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0)) # No color on face x axis face
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0)) # No color on face y axis face
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0)) # No color on face z axis face
    
    plt.savefig("/home/willow/Python/Plots/simulation_result.png", figsize=(8, 10), dpi=600)
    plt.show()
    
    directory_path = "/home/willow/Python/three_star_csv_files"
    os.makedirs(directory_path, exist_ok=True)

    file_name_1 = f"system_states_{datetime.datetime.now().strftime('%Y%m%d_T%H_%M_%S')}.csv"
    file_path = os.path.join(directory_path, file_name_1)

    # first index is iteration
    # second index is planet/sun
    # third index is states

    states = np.hstack((SolarSystem.times, np.reshape(history, (history.shape[0], history.shape[1] * history.shape[2])), SolarSystem.temperature))

    state_names = ["time_s"]
    for name in names:
        for val in ["x", "y", "z", "vx", "vy", "vz"]:
            state_names.append(name + f"_{val}")
    state_names.append("T_K")

    #np.savetxt(file_path, states, delimiter=",", fmt="%.2f", header=",".join(state_names), comments='')

    print(f"Simulation {p} saved to {file_name_1}")