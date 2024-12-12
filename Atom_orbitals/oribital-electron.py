import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm, genlaguerre, factorial

def hydrogenic_wavefunction(n, l, m, r, theta, phi, Z):
    # Radial part
    rho = 2 * Z * r / n
    norm_radial = np.sqrt((2 * Z / n)**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    radial = norm_radial * np.exp(-rho / 2) * rho**l * genlaguerre(n - l - 1, 2 * l + 1)(rho)
    
    # Angular part
    angular = sph_harm(m, l, phi, theta)
    
    # Total wavefunction
    psi = radial * angular
    return psi

def plot_orbital_3d(n, l, m, Z):
    # Create spherical grid with reduced size
    r = np.linspace(0, 20, 100)
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)
    r, theta, phi = np.meshgrid(r, theta, phi)
    
    # Compute probability density
    psi = hydrogenic_wavefunction(n, l, m, r, theta, phi, Z)
    prob_density = np.abs(psi)**2
    
    # Convert to Cartesian coordinates
    X = r * np.sin(theta) * np.cos(phi)
    Y = r * np.sin(theta) * np.sin(phi)
    Z_grid = r * np.cos(theta)
    
    # Flatten arrays
    X = X.flatten()
    Y = Y.flatten()
    Z_grid = Z_grid.flatten()
    prob_density = prob_density.flatten()
    
    # Filter data to plot regions with higher probability density
    threshold = np.percentile(prob_density, 99)
    indices = prob_density > threshold
    X = X[indices]
    Y = Y[indices]
    Z_grid = Z_grid[indices]
    prob_density = prob_density[indices]
    
    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(X, Y, Z_grid, c=prob_density, cmap='viridis', s=10)
    ax.set_xlabel('x (a₀)')
    ax.set_ylabel('y (a₀)')
    ax.set_zlabel('z (a₀)')
    ax.set_title(f'Atom Z={Z}, Orbital (n={n}, l={l}, m={m})')
    fig.colorbar(sc, label='Probability Density')
    plt.show()

if __name__ == "__main__":
    Z = 9
    n, l, m = 5, 2, -1
    
    psi = hydrogenic_wavefunction(n, l, m, 1.0, np.pi / 2, 0.0, Z)
    print(f"Orbital value at (r=1.0, θ=π/2, φ=0): {psi}")
    plot_orbital_3d(n, l, m, Z)