import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
import math

# Parameters for the Earth and orbit
earth_radius = 1  # Radius of the Earth (arbitrary units)
orbit_radius = 2  # Radius of the orbit (arbitrary units)
inclination_angle = 30  # Inclination angle in degrees

# Create the Earth
def create_sphere(radius, resolution):
    phi = np.linspace(0, np.pi, resolution)  # Polar angle
    theta = np.linspace(0, 2 * np.pi, resolution)  # Azimuthal angle
    phi, theta = np.meshgrid(phi, theta)

    # Cartesian coordinates
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    return x, y, z

x, y, z = create_sphere(6400000, 100)

def perifiocal_point( mu, a, e, true_anomaly ):

    #h**2/mu * 1/(1+e*cos(theta))( Y cos(theta) + X * sin( theta ))

    r = a * (1 - e**2) / (1 + e * np.cos( true_anomaly ))
    x =  r * np.cos( true_anomaly )
    y =  r * np.sin( true_anomaly )
    z = np.zeros( shape = len( true_anomaly ) )

    v_x = math.sqrt( mu / ( a * (1 - e ** 2 ) )) * - np.sin( true_anomaly )
    v_y = math.sqrt( mu / ( a * (1 - e ** 2 ) )) * ( e + np.cos( true_anomaly ))
    v_z = np.zeros( shape = len( true_anomaly ) )

    r_perifocal = np.array( [ x, y, z ] )
    v_perifocal = np.array( [ v_x, v_y, v_z ] )
    return r_perifocal,v_perifocal

def perifocal_to_eci( r, v, Omega_RAAN_degrees, Inclination_degrees, argument_periapsis_degrees ):

    RZ_omega = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[0,0,Omega_RAAN_degrees], degrees=True)
    RX_inclination = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[Inclination_degrees,0,0], degrees= True)
    RZ_argument_periapsis = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[0,0,argument_periapsis_degrees], degrees=True)

    R = RZ_omega.as_matrix() @ RX_inclination.as_matrix() @ RZ_argument_periapsis.as_matrix()

    r_ECI = R @ r
    v_ECI = R @ v

    return r_ECI, v_ECI

# View it.
from mayavi import mlab
s = mlab.mesh( x, y, z, color= (0.5,0.5,0.5) )

r,v = perifiocal_point( 3.986004418e14, 30000000, 0.002, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r_ECI, v_ECI = perifocal_to_eci(r,v, 100, 80, 10)
mlab.plot3d( r_ECI[0], r_ECI[1], r_ECI[2], tube_radius = None, color = (0,0.5,1) )

r2,v2 = perifiocal_point( 3.986004418e14, 20000000, 0.25, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r2_ECI, v2_ECI = perifocal_to_eci(r2,v2, 100, 80, 10)
mlab.plot3d( r2_ECI[0], r2_ECI[1], r2_ECI[2], tube_radius = None, color = (1,0.5,0) )

r2,v2 = perifiocal_point( 3.986004418e14, 30000000, (30000000 - 20000000)/(30000000 + 20000000), np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r2_ECI, v2_ECI = perifocal_to_eci(r2,v2, 100, 80, 10)
mlab.plot3d( r2_ECI[0], r2_ECI[1], r2_ECI[2], tube_radius = None, color = (0,1,0.5) )

s.scene.background = (1, 1, 1)  # for white
mlab.show()

# import matplotlib.pyplot
# import numpy as np

# fig = matplotlib.pyplot.figure()
# ax = fig.add_subplot(projection='3d')

# # Make data
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x = 6400000 * np.outer(np.cos(u), np.sin(v))
# y = 6400000 * np.outer(np.sin(u), np.sin(v))
# z = 6400000 * np.outer(np.ones(np.size(u)), np.cos(v))

# # Plot the surface
# ax.plot_surface(x, y, z )
# ax.plot(  *r_ECI )
# # Set an equal aspect ratio
# ax.set_aspect('equal')

# matplotlib.pyplot.show()

# import numpy
# from mayavi.mlab import *

# def test_plot3d():
#     """Generates a pretty set of lines."""
#     n_mer, n_long = 6, 11
#     dphi = np.pi / 1000.0
#     phi = np.arange(0.0, 2 * np.pi + 0.5 * dphi, dphi)
#     mu = phi * n_mer
#     x = np.cos(mu) * (1 + np.cos(n_long * mu / n_mer) * 0.5)
#     y = np.sin(mu) * (1 + np.cos(n_long * mu / n_mer) * 0.5)
#     z = np.sin(n_long * mu / n_mer) * 0.5

#     l = plot3d(x, y, z, np.sin(mu), tube_radius=0, colormap='Spectral')
#     return l
# test_plot3d()
# show()
# earth = plt.Circle((0, 0), earth_radius, color='blue', label='Earth')

# # Create the orbit
# theta = np.linspace(0, 2 * np.pi, 100)  # Angle for the orbit
# x_orbit = orbit_radius * np.cos(theta)
# y_orbit = orbit_radius * np.sin(theta)

# # Rotate the orbit by the inclination angle
# inclination_radians = np.radians(inclination_angle)
# x_orbit_inclined = x_orbit
# y_orbit_inclined = y_orbit * np.cos(inclination_radians) - x_orbit * np.sin(inclination_radians)



# # Plotting
# fig, ax = plt.subplots(figsize=(8, 8))
# ax.add_artist(earth)
# ax.plot(x_orbit_inclined, y_orbit_inclined, color='orange', label='Inclined Orbit')
# ax.set_xlim(-3, 3)
# ax.set_ylim(-3, 3)
# ax.set_aspect('equal', adjustable='box')
# ax.axhline(0, color='gray', lw=0.5, ls='--')
# ax.axvline(0, color='gray', lw=0.5, ls='--')
# ax.set_title('Inclined Orbit Around the Earth')
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.legend()
# plt.grid()
# plt.show()