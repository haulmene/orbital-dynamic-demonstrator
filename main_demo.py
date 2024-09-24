import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
import math

import time

# Parameters for the Earth and orbit
earth_radius = 6_400_000  # Radius of the Earth (arbitrary units)
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

def calculate_R_perifocal_to_ECI(Omega_RAAN_degrees, Inclination_degrees, argument_periapsis_degrees):

    RZ_omega = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[0,0,Omega_RAAN_degrees], degrees=True)
    RX_inclination = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[Inclination_degrees,0,0], degrees= True)
    RZ_argument_periapsis = scipy.spatial.transform.Rotation.from_rotvec( rotvec=[0,0,argument_periapsis_degrees], degrees=True)

    R_perifocal_to_ECI = RZ_omega.as_matrix() @ RX_inclination.as_matrix() @ RZ_argument_periapsis.as_matrix()

    return R_perifocal_to_ECI

# ra = 30000000
# rp = 20000000
mu = 3.986004418e14

# Omega_RAAN = 120
# inclination = 80
# omega_argument_perigee = 10
# true_anomaly = 20.0 / 180.0 * np.pi

# e = (ra - rp)/(ra + rp)
# a = (ra + rp) / 2

from mayavi import mlab
s = mlab.mesh( x, y, z, color= (0.5,0.5,0.5) )

# r,v = perifiocal_point( mu, a, e, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
# r_ECI, v_ECI = perifocal_to_eci( r, v, Omega_RAAN, inclination, omega_argument_perigee )

# r_ECI = np.array([-10370354.754472725, 16755833.74297841, 3420201.4332566867])
# v_ECI = np.array([-299.6592399981582, -1153.5907820138782, 4742.937819067031])
r_ECI = np.array([-10265508.28871291,  14266052.74092895,   9965375.59516136])
v_ECI = np.array([  440.14817962, -2266.16768951,  4264.26141223])
#119.9999999999735 79.99999999995487 9.999999999774763 0.1999999999987604 24999999.999956783 0.3500381786769686
#
# R_perifocal_to_ECI = calculate_R_perifocal_to_ECI( Omega_RAAN, inclination, omega_argument_perigee )
# r_periapsis = R_perifocal_to_ECI @ np.array([ a * (1 - e ), 0, 0 ])
# r_apoapsis = R_perifocal_to_ECI @ np.array([ -a * (1 + e ), 0, 0 ])

# node line
h = np.cross( r_ECI , v_ECI )
nodes = np.cross( h , [0, 0, 1] )
nodes = nodes / np.linalg.norm( nodes )

# draw whole orbit 
mlab.plot3d( r_ECI[0], r_ECI[1], r_ECI[2], tube_radius = None, color = (0,0.5,1) )
# mlab.plot3d( r[0], r[1], r[2], tube_radius = None, color = (0,0.7,1) )
# draw apse line
# mlab.plot3d( [r_periapsis[0], r_apoapsis[0]], [r_periapsis[1], r_apoapsis[1]], [r_periapsis[2], r_apoapsis[2]],  tube_radius = None, color = (0,0.5,1), )

# draw the equator
r_equator, v_equator = perifiocal_point( mu, earth_radius, 0, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
mlab.plot3d( r_equator[0], r_equator[1], r_equator[2], tube_radius = None, color = (0.7,0.5,0))

# draw the equator plane
# grid_1, grid_2 = np.mgrid[0:2, 0:2]
# mlab.mesh( 10 * earth_radius * ( 0.5 - grid_1 ), 10 * earth_radius * ( 0.5 - grid_2 ) , [[0,0],[0,0]], opacity = 0.2 )

# draw vernal equinox
# mlab.plot3d( [0, 5 * earth_radius ], [0, 0], [0, 0 ], tube_radius = None, color = (0.2,0.7,0.2), )

# draw frame
# mlab.plot3d( [0, 0 ], [-5 * earth_radius, 5 * earth_radius], [ 0, 0 ], tube_radius = None, color = (0.7,0.7,0.7), )
# mlab.plot3d( [0, 0 ], [0, 0], [ -5 * earth_radius, 5 * earth_radius ], tube_radius = None, color = (0.7,0.7,0.7), )

import ElementsPropagation
import ElementsConvertion

# true_anomaly_start = true_anomaly
# r_start,v_start = perifiocal_point( mu, a, e, [true_anomaly_start]  )
r_ECI_start = r_ECI
v_ECI_start = v_ECI

# radius vector
# mlab.plot3d( [0, r_ECI_start[0] ], [0, r_ECI_start[1] ], [ 0, r_ECI_start[2] ], tube_radius = None, color = (0.2,0.2,0.2), )
# mlab.quiver3d( r_ECI_start[0][0], r_ECI_start[1][0], r_ECI_start[2][0], 10*v_ECI_start[0][0], 10*v_ECI_start[1][0], 10*v_ECI_start[2][0] )
# mlab.quiver3d( 0, -5_000_000, 0, 0, 0, 10000, extent = [-10_000_000, 10_000_000, -10_000_000, 10_000_000, -10_000_000, 10_000_000] )

# r_current,v_current = perifiocal_point( mu, a, e, [true_anomaly_start]  )
# r_ECI_current, v_ECI_current = perifocal_to_eci(r_current,v_current, Omega_RAAN, inclination, omega_argument_perigee)

# satellite start
mlab.points3d( r_ECI_start[0], r_ECI_start[1], r_ECI_start[2], 1, color = (0,0,0), scale_factor= 300000 )
mlab.plot3d( [r_ECI_start[0], r_ECI_start[0] + 1000 * v_ECI_start[0]], [r_ECI_start[1], r_ECI_start[1] + 1000 * v_ECI_start[1]], [r_ECI_start[2], r_ECI_start[2] + 1000 * v_ECI_start[2]], color = (0,0,0), tube_radius = None )

# v and r -> elements
Omega_RAAN, inclination, omega_argument_perigee, eccentricity, semi_major_axis, true_anomaly = ElementsConvertion.rv2Elements( [r_ECI_start], [v_ECI_start], mu )
print("rv -> elements")
print( 180 / np.pi * Omega_RAAN, 180 / np.pi * inclination, 180 / np.pi * omega_argument_perigee, eccentricity, semi_major_axis, true_anomaly )

# orbit
r_periapsis = semi_major_axis * ( 1 - eccentricity )
r_apoapsis = semi_major_axis * ( 1 + eccentricity )

R_perifocal_to_ECI = calculate_R_perifocal_to_ECI( Omega_RAAN * 180 / np.pi, inclination * 180 / np.pi, omega_argument_perigee * 180 / np.pi )
r_periapsis = R_perifocal_to_ECI @ np.array([ r_periapsis, 0, 0 ])
r_apoapsis = R_perifocal_to_ECI @ np.array([ -r_apoapsis, 0, 0 ])

r,v = perifiocal_point( mu, semi_major_axis, eccentricity, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r_ECI_orbit, v_ECI_orbit = perifocal_to_eci( r, v, Omega_RAAN * 180 /np.pi, inclination * 180 /np.pi, omega_argument_perigee * 180 /np.pi)

# draw orbit trajectory
# mlab.plot3d( r_ECI_orbit[0], r_ECI_orbit[1], r_ECI_orbit[2], tube_radius = None, color = (0,0.5,1) )

# draw node line
# mlab.plot3d( [-2 * semi_major_axis * nodes[0], 2 * semi_major_axis * nodes[0]], [-2 * semi_major_axis * nodes[1], 2 * semi_major_axis * nodes[1]], [-2 * semi_major_axis * nodes[2], 2 * semi_major_axis * nodes[2] ], tube_radius = None, color = (0.7,0.2,0.2), )

# draw radius vector
# mlab.plot3d( [0, r_ECI_start[0] ], [0, r_ECI_start[1] ], [ 0, r_ECI_start[2] ], tube_radius = None, color = (0.2,0.2,0.2), )

# draw apse line
# mlab.plot3d( [r_periapsis[0], r_apoapsis[0]], [r_periapsis[1], r_apoapsis[1]], [r_periapsis[2], r_apoapsis[2]],  tube_radius = None, color = (0,0.5,1), )

# propagate

time_orbital = 3000

true_anomaly_current = ElementsPropagation.PropagateTrueAnomaly(   true_anomaly = true_anomaly,
                                                            e = eccentricity, 
                                                            T_period = ElementsConvertion.period( mu, semi_major_axis ),
                                                            delta_time= time_orbital )

r_current, v_current = perifiocal_point( mu, semi_major_axis, eccentricity, [true_anomaly_current]  )
r_ECI_current, v_ECI_current = perifocal_to_eci(r_current,v_current, Omega_RAAN * 180 / np.pi, inclination * 180 / np.pi, omega_argument_perigee * 180 / np.pi)
# mlab.points3d( r_ECI_current[0], r_ECI_current[1], r_ECI_current[2], 1, color = (0,0,0), scale_factor= 300000 )
# mlab.plot3d( [r_ECI_current[0], r_ECI_current[0] + 1000 * v_ECI_current[0]], [r_ECI_current[1], r_ECI_current[1] + 1000 * v_ECI_current[1]], [r_ECI_current[2], r_ECI_current[2] + 1000 * v_ECI_current[2]], color = (0,0,0), tube_radius = None )

# satellite flight
# pts = mlab.points3d( r_ECI_current[0], r_ECI_current[1], r_ECI_current[2], 1, color = (1,0,0), scale_factor= 300000 )

s.scene.background = (1, 1, 1)  # white


# print( "period: ", ElementsConvertion.period( mu, a ) )

# v and r -> elements
Omega_RAAN, inclination, omega_argument_perigee, eccentricity, semi_major_axis, true_anomaly = ElementsConvertion.rv2Elements( [r_ECI_start], [v_ECI_start], mu )
print("rv -> elements")
print( 180 / np.pi * Omega_RAAN, 180 / np.pi * inclination, 180 / np.pi * omega_argument_perigee, eccentricity, semi_major_axis, true_anomaly )
# print("original")
# print( Omega_RAAN, inclination, omega_argument_perigee, e, a, true_anomaly )

# global time_orbital

# time_orbital = 0

# animate orbital motion

# @mlab.animate(delay=100)
# def anim():
#     global time_orbital
#     while True:
#         f = mlab.gcf()

#         time_orbital += 100*0.1

#         true_anomaly_current = ElementsPropagation.PropagateTrueAnomaly(   true_anomaly = true_anomaly,
#                                                                     e = eccentricity, 
#                                                                     T_period = ElementsConvertion.period( mu, semi_major_axis ),
#                                                                     delta_time= time_orbital )

#         r_current, v_current = perifiocal_point( mu, semi_major_axis, eccentricity, [true_anomaly_current]  )
#         r_ECI_current, v_ECI_current = perifocal_to_eci(r_current,v_current, Omega_RAAN, inclination, omega_argument_perigee)
#         pts.mlab_source.set( x = r_ECI_current[0], y = r_ECI_current[1], z = r_ECI_current[2] )

#         print( true_anomaly_current )
#         yield

# anim()
mlab.show()