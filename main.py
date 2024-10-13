import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
import math

import time

from J2AnalyticalPertrubations import Omega_RAAN_derivative, omega_argument_perigee_derivative

# Parameters for the Earth and orbit
earth_radius = 6_378_150  # Radius of the Earth (arbitrary units)
orbit_radius = 2  # Radius of the orbit (arbitrary units)
inclination_angle = 30  # Inclination angle in degrees
J2 = 1082.64e-6
mu = 3.986004418e14

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

# calculating an orbit

ra = 30000000
rp = 20000000


global Omega_RAAN
Omega_RAAN = 120
global inclination
inclination = 80
global omega_argument_perigee
omega_argument_perigee = 10
global true_anomaly
true_anomaly = 20.0 / 180.0 * np.pi

e = (ra - rp)/(ra + rp)
a = (ra + rp) / 2

from mayavi import mlab
s = mlab.mesh( x, y, z, color= (0.5,0.5,0.5) )

r,v = perifiocal_point( mu, a, e, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r_ECI, v_ECI = perifocal_to_eci( r, v, Omega_RAAN, inclination, omega_argument_perigee )

R_perifocal_to_ECI = calculate_R_perifocal_to_ECI( Omega_RAAN, inclination, omega_argument_perigee )
r_periapsis = R_perifocal_to_ECI @ np.array([ a * (1 - e ), 0, 0 ])
r_apoapsis = R_perifocal_to_ECI @ np.array([ -a * (1 + e ), 0, 0 ])

# node line
h = np.cross( r_ECI[:, 0] , v_ECI[:, 0] )
nodes = np.cross( h , [0, 0, 1] )
nodes = nodes / np.linalg.norm( nodes )

# draw whole orbit 
orbi_dynamic = mlab.plot3d( r_ECI[0], r_ECI[1], r_ECI[2], tube_radius = None, color = (0.5,0.0,1) )

# draw whole orbit 
mlab.plot3d( r_ECI[0], r_ECI[1], r_ECI[2], tube_radius = None, color = (0,0.5,1) )

# mlab.plot3d( r[0], r[1], r[2], tube_radius = None, color = (0,0.7,1) )
# draw apse line
mlab.plot3d( [r_periapsis[0], r_apoapsis[0]], [r_periapsis[1], r_apoapsis[1]], [r_periapsis[2], r_apoapsis[2]],  tube_radius = None, color = (0,0.5,1), )

# draw the equator
r_equator, v_equator = perifiocal_point( mu, earth_radius, 0, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
mlab.plot3d( r_equator[0], r_equator[1], r_equator[2], tube_radius = None, color = (0.7,0.5,0))

# draw the equator plane
grid_1, grid_2 = np.mgrid[0:2, 0:2]
mlab.mesh( 10 * earth_radius * ( 0.5 - grid_1 ), 10 * earth_radius * ( 0.5 - grid_2 ) , [[0,0],[0,0]], opacity = 0.2 )

# draw node line
mlab.plot3d( [-ra * nodes[0], ra * nodes[0]], [-ra * nodes[1], ra * nodes[1]], [-ra * nodes[2], ra * nodes[2] ], tube_radius = None, color = (0.7,0.2,0.2), )

# draw vernal equinox
mlab.plot3d( [0, 5 * earth_radius ], [0, 0], [0, 0 ], tube_radius = None, color = (0.2,0.7,0.2), )

# draw frame
mlab.plot3d( [0, 0 ], [-5 * earth_radius, 5 * earth_radius], [ 0, 0 ], tube_radius = None, color = (0.7,0.7,0.7), )
mlab.plot3d( [0, 0 ], [0, 0], [ -5 * earth_radius, 5 * earth_radius ], tube_radius = None, color = (0.7,0.7,0.7), )

import ElementsPropagation
import ElementsConvertion

true_anomaly_start = true_anomaly
r_start,v_start = perifiocal_point( mu, a, e, [true_anomaly_start]  )
r_ECI_start, v_ECI_start = perifocal_to_eci(r_start, v_start, Omega_RAAN, inclination, omega_argument_perigee)

# radius vector
mlab.plot3d( [0, r_ECI_start[0][0] ], [0, r_ECI_start[1][0] ], [ 0, r_ECI_start[2][0] ], tube_radius = None, color = (0.2,0.2,0.2), )

r_current,v_current = perifiocal_point( mu, a, e, [true_anomaly_start]  )
r_ECI_current, v_ECI_current = perifocal_to_eci(r_current,v_current, Omega_RAAN, inclination, omega_argument_perigee)

# satellite start
mlab.points3d( r_ECI_start[0], r_ECI_start[1], r_ECI_start[2], 1, color = (0,0,1), scale_factor= 800000 )

# satellite flight
pts = mlab.points3d( r_ECI_current[0], r_ECI_current[1], r_ECI_current[2], 1, color = (1,0,0), scale_factor= 800000 )

s.scene.background = (1, 1, 1)  # white


print( "period: ", ElementsConvertion.period( mu, a ) )

global time_orbital

time_orbital = 0

# animate orbital motion

@mlab.animate(delay=100)
def anim():
    global time_orbital, Omega_RAAN, omega_argument_perigee
    while True:
        f = mlab.gcf()
        
        delta_t = 10000
        time_orbital += delta_t

        true_anomaly_current = ElementsPropagation.PropagateTrueAnomaly(   true_anomaly = true_anomaly,
                                                                    e = e, 
                                                                    T_period = ElementsConvertion.period( mu, a ),
                                                                    delta_time= time_orbital )

        dOmega_RAAN_dt = Omega_RAAN_derivative( J2=J2, R_Earth= earth_radius, mu_earth= mu, a=a, e = e, i = inclination / 180 * math.pi )
        domega_argument_periapsis_dt = omega_argument_perigee_derivative(  J2=J2, R_Earth= earth_radius, mu_earth= mu, a=a, e = e, i = inclination / 180 * math.pi )

        Omega_RAAN_new = Omega_RAAN / 180 * math.pi + dOmega_RAAN_dt * delta_t
        omega_argument_perigee_new  = omega_argument_perigee / 180 * math.pi + domega_argument_periapsis_dt * delta_t

        #account for argument perigee change
        true_anomaly_new = true_anomaly_current - domega_argument_periapsis_dt * delta_t


        r_current, v_current = perifiocal_point( mu, a, e, [true_anomaly_new]  )
        r_ECI_current, v_ECI_current = perifocal_to_eci(r_current,v_current, Omega_RAAN_new / math.pi * 180, inclination, omega_argument_perigee_new / math.pi * 180 )
        pts.mlab_source.set( x = r_ECI_current[0], y = r_ECI_current[1], z = r_ECI_current[2] )

        Omega_RAAN              = Omega_RAAN_new / math.pi * 180
        omega_argument_perigee  = omega_argument_perigee_new / math.pi * 180

        #orbit 
        r,v = perifiocal_point( mu, a, e, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
        r_ECI, v_ECI = perifocal_to_eci( r, v, Omega_RAAN_new / math.pi * 180, inclination, omega_argument_perigee_new / math.pi * 180 )
        orbi_dynamic.mlab_source.reset( x = r_ECI[0], y = r_ECI[1], z = r_ECI[2] )

        print( true_anomaly_current )
        print( "dOmega", dOmega_RAAN_dt * delta_t )
        print( "domega", domega_argument_periapsis_dt * delta_t )
        print( "Omega", Omega_RAAN )
        print( "omega", omega_argument_perigee )

        yield

anim()
mlab.show()