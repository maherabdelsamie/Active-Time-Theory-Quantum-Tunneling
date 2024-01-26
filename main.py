import numpy as np
import matplotlib.pyplot as plt

# Constants and System Parameters
hbar = 1.0  # Reduced Planck's constant
m = 1.0  # Particle mass
N = 2000  # Number of spatial points
x = np.linspace(-10, 10, N)  # Spatial coordinates
dx = x[1] - x[0]  # Spatial step size
k = np.fft.fftfreq(N, d=dx) * 2 * np.pi  # Fourier space

# Potential Barrier Parameters
V0 = 1.0  # Height of the potential barrier
L = 2  # Width of the potential barrier

# Time Evolution Parameters
dt = 0.01  # Time step
time_steps = 1000  # Total number of time steps

# Initial Wave Packet Parameters
x0 = -5.0  # Initial position of the center
k0 = 5.0  # Initial wave number
sigma = 1.0  # Width of the wave packet

# Initial wave function as a Gaussian wave packet
psi = np.exp(-0.5 * ((x - x0) / sigma)**2) * np.exp(1j * k0 * x)
psi /= np.sqrt(np.sum(abs(psi)**2) * dx)  # Normalize the wave function

# Potential barrier
V = np.zeros_like(x)
barrier_center = N // 2
barrier_start = barrier_center - int(L / (2 * dx))
barrier_end = barrier_center + int(L / (2 * dx))
V[barrier_start:barrier_end] = V0

# Active time functions
def ath_generative(t):
    return 0.1 * np.random.randn()

def ath_directive(psi, t):
    return np.sin(t) * abs(psi)

def ath_adaptive(t):
    return dt * (1 + 0.1 * np.cos(t))

# Function to compute the probability of finding the particle on the right side of the barrier
def compute_transmission_probability(psi, barrier_index):
    return np.sum(abs(psi[barrier_index:]) ** 2) * dx

# Store transmission probabilities
trans_probs_active = []
trans_probs_control = []

# Define kinetic energy operator in Fourier space for pseudo-spectral method
kinetic_energy_operator = -0.5 * hbar**2 * k**2 / m

# Simulation with active time effects using Pseudo-Spectral method
psi_active = psi.copy()
for t in range(time_steps):
    current_time = t * dt
    # Apply active time effects to the potential
    V_active = V * (1 + ath_generative(current_time) + ath_directive(psi_active, current_time))
    # Compute the Hamiltonian in Fourier space
    H_active = kinetic_energy_operator + np.fft.fft(V_active)
    # Evolve the wave function in Fourier space
    psi_active = np.fft.ifft(np.fft.fft(psi_active) * np.exp(-1j * ath_adaptive(current_time) * H_active * dt / hbar))
    # Normalize the wave function
    psi_active /= np.sqrt(np.sum(abs(psi_active)**2) * dx)
    # Compute transmission probability and append to list
    trans_prob_active = compute_transmission_probability(psi_active, barrier_end)
    trans_probs_active.append(trans_prob_active)

# Control simulation without active time effects using Pseudo-Spectral method
psi_control = psi.copy()
for t in range(time_steps):
    # Compute the Hamiltonian in Fourier space
    H_control = kinetic_energy_operator + np.fft.fft(V)
    # Evolve the wave function in Fourier space
    psi_control = np.fft.ifft(np.fft.fft(psi_control) * np.exp(-1j * H_control * dt / hbar))
    # Normalize the wave function
    psi_control /= np.sqrt(np.sum(abs(psi_control)**2) * dx)
    # Compute transmission probability and append to list
    trans_prob_control = compute_transmission_probability(psi_control, barrier_end)
    trans_probs_control.append(trans_prob_control)

# Plot the results
# Plotting the tunneling rates from the previous calculation
plt.figure(figsize=(14, 6))
plt.plot(range(time_steps), trans_probs_active, label='With Active Time Effects')
plt.plot(range(time_steps), trans_probs_control, label='Without Active Time Effects', alpha=0.7)
plt.xlabel('Time Step')
plt.ylabel('Tunneling Rate')
plt.title('Tunneling Rate Over Time Using Pseudo-Spectral Method')
plt.legend()
plt.show()
