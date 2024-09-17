import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
import math

import time

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




mu = 3.986004418e14
e = (30000000 - 20000000)/(30000000 + 20000000)
a = (30000000 + 20000000) / 2
Omega_RAAN = 100
inclination = 80
omega_argument_perigee = 10
true_anomaly = 20.0 / 180.0 * np.pi

from mayavi import mlab
s = mlab.mesh( x, y, z, color= (0.5,0.5,0.5) )

r,v = perifiocal_point( mu, a, e, np.linspace( 0, 360, 360 ) / 180.0 * np.pi  )
r_ECI, v_ECI = perifocal_to_eci( r, v, Omega_RAAN, inclination, omega_argument_perigee )
mlab.plot3d( r_ECI[0], r_ECI[1], r_ECI[2], tube_radius = None, color = (0,0.5,1) )

import ElementsPropagation
import ElementsConvertion

r1,v1 = perifiocal_point( mu, a, e, [true_anomaly]  )
r_ECI1, v_ECI1 = perifocal_to_eci(r1,v1, Omega_RAAN, inclination, omega_argument_perigee)

r2,v2 = perifiocal_point( mu, a, e, [true_anomaly]  )
r_ECI2, v_ECI2 = perifocal_to_eci(r2,v2, Omega_RAAN, inclination, omega_argument_perigee)

mlab.points3d( r_ECI1[0], r_ECI1[1], r_ECI1[2], 1, color = (0,0,1), scale_factor= 300000 )


pts = mlab.points3d( r_ECI2[0], r_ECI2[1], r_ECI2[2], 1, color = (1,0,0), scale_factor= 300000 )

s.scene.background = (1, 1, 1)  # white


print( "period: ", ElementsConvertion.period( mu, a ) )
global time_orbital

time_orbital = 0

@mlab.animate(delay=100)
def anim():
    global time_orbital
    while True:
        f = mlab.gcf()

        time_orbital += 100*0.1

        true_anomaly2 = ElementsPropagation.PropagateTrueAnomaly(   true_anomaly = true_anomaly,
                                                                    e = e, 
                                                                    T_period = ElementsConvertion.period( mu, a ),
                                                                    delta_time= time_orbital )

        r2,v2 = perifiocal_point( mu, a, e, [true_anomaly2]  )
        r_ECI2, v_ECI2 = perifocal_to_eci(r2,v2, Omega_RAAN, inclination, omega_argument_perigee)
        pts.mlab_source.set( x = r_ECI2[0], y = r_ECI2[1], z = r_ECI2[2] )

        print( true_anomaly2 )
        yield

anim()
mlab.show()