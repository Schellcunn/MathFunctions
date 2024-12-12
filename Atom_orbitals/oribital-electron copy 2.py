import numpy as np
from scipy.special import genlaguerre, sph_harm, factorial
import math
import matplotlib.pyplot as plt
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Bohr radius in atomic units
a0 = 1.0

def radial_wavefunction(n, l, r):
    """
    Compute the radial part of hydrogenic wavefunctions R_{n,l}(r).
    Parameters:
        n (int): Principal quantum number
        l (int): Azimuthal quantum number
        r (float or ndarray): Radial coordinate (units of a0)
    Returns:
        R_{n,l}(r)
    """
    
    if n <= l:
        raise ValueError("Must have n > l for hydrogen-like orbitals.")

    p = n - l - 1
    alpha = 2*l + 1
    laguerre_poly = genlaguerre(p, alpha)
    
    num = (2.0/(n*a0))**3 * factorial(n-l-1)
    den = 2*n*factorial(n+l)
    N = math.sqrt(num/den)
    
    x = 2.0*r/(n*a0)
    R = N * (r/(n*a0))**l * np.exp(-r/(n*a0)) * laguerre_poly(x)
    return R

def hydrogenic_wavefunction(n, l, m, r, theta, phi):
    """
    Compute the full hydrogenic wavefunction psi_{n,l,m}(r, theta, phi).
    """
    if abs(m) > l:
        raise ValueError("Must have |m| ≤ l")
        
    Y_lm = sph_harm(m, l, phi, theta)
    R_nl = radial_wavefunction(n, l, r)
    psi = R_nl * Y_lm
    return psi

def get_slater_screening(n, l, Z):
    """
    Calculate screening constant using Slater's rules
    """
    # Simplified Slater's rules
    if n <= 3:
        s = 0.35  # Inner shell screening
        return Z - (Z-1)*s
    else:
        s = 0.35  # Outer shell screening
        return Z - (Z-1)*s - 0.85*(n-3)

def radial_wavefunction_with_screening(n, l, r, Z):
    """
    Modified radial wavefunction with electron screening
    """
    if n <= l:
        raise ValueError("Must have n > l for atomic orbitals.")

    # Calculate effective nuclear charge
    Z_eff = get_slater_screening(n, l, Z)
    
    p = n - l - 1
    alpha = 2*l + 1
    laguerre_poly = genlaguerre(p, alpha)
    
    # Modified to use effective nuclear charge
    num = (2*Z_eff/(n*a0))**3 * factorial(n-l-1)
    den = 2*n*factorial(n+l)
    N = math.sqrt(num/den)
    
    x = 2*Z_eff*r/(n*a0)
    R = N * (2*Z_eff*r/(n*a0))**l * np.exp(-Z_eff*r/(n*a0)) * laguerre_poly(x)
    
    # Add correlation correction (simplified)
    correlation_factor = 1.0 - 0.1*(Z-1)/Z
    R *= correlation_factor
    
    return R

def hydrogenic_wavefunction_with_interactions(n, l, m, r, theta, phi, Z):
    """
    Modified wavefunction including electron interactions
    """
    if abs(m) > l:
        raise ValueError("Must have |m| ≤ l")
        
    Y_lm = sph_harm(m, l, phi, theta)
    R_nl = radial_wavefunction_with_screening(n, l, r, Z)
    psi = R_nl * Y_lm
    return psi

def plot_orbital_3d_with_interactions(n, l, m, Z, points=50, num_surfaces=4):
    """
    Plot atomic orbital with electron interactions using marching cubes algorithm
    """
    # Dynamic size scaling
    size = (n**2 + 2*n) / np.sqrt(Z) * 2
    points = max(40, min(100, points + 5*(n + l)))
    
    # Create 3D grid
    x = np.linspace(-size, size, points)
    y = np.linspace(-size, size, points)
    z = np.linspace(-size, size, points)
    X, Y, Z_grid = np.meshgrid(x, y, z)
    
    # Calculate wavefunction
    r = np.sqrt(X**2 + Y**2 + Z_grid**2)
    theta = np.arccos(Z_grid/(r + 1e-10))
    phi = np.arctan2(Y, X)
    psi = hydrogenic_wavefunction_with_interactions(n, l, m, r, theta, phi, Z)
    probability = np.abs(psi)**2
    
    # Create figure with proper spacing for colorbar
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Calculate isovalues
    max_prob = np.max(probability)
    min_prob = max_prob * 0.02
    iso_values = np.logspace(np.log10(min_prob), np.log10(max_prob), num_surfaces)
    
    # Plot isosurfaces using marching cubes seems to give better surface results
    colors = plt.cm.viridis(np.linspace(0, 1, len(iso_values)))
    
    for idx, iso_val in enumerate(iso_values):
        verts, faces, _, _ = measure.marching_cubes(probability, iso_val)
        verts = verts * (2*size/points) - size
        mesh = Poly3DCollection(verts[faces])
        mesh.set_facecolor(colors[idx])
        mesh.set_alpha(0.3)
        ax.add_collection3d(mesh)
    
    # Plot settings
    ax.set_xlabel('x (a₀)')
    ax.set_ylabel('y (a₀)')
    ax.set_zlabel('z (a₀)')
    ax.set_title(f'Atom Z={Z}, Orbital (n={n}, l={l}, m={m})')
    ax.set_box_aspect([1,1,1])
    
    # Set view limits
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_zlim(-size, size)
    
    # Add colorbar with proper spacing
    plt.subplots_adjust(right=0.85)
    cax = plt.axes([0.9, 0.15, 0.03, 0.7])
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    sm.set_array(iso_values)
    fig.colorbar(sm, cax=cax, label='|ψ|²')
    
    plt.show()

if __name__ == "__main__":
    Z = 9
    n, l, m = 5, 2, -1
    r = 1.0
    theta = np.pi/2
    phi = 0.0

    psi = hydrogenic_wavefunction_with_interactions(n, l, m, r, theta, phi, Z)
    print(f"Orbital value at (r={r}, θ={theta}, φ={phi}): {psi}")
    plot_orbital_3d_with_interactions(n, l, m, Z)