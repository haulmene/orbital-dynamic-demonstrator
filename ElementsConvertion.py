# based on https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html

import numpy as np
import math

def rv2Elements( r, v, mu ):

    h_vec = np.cross(r, v)
    h = np.linalg.norm(h_vec)

    inclination = np.arccos(h_vec[0][2] / h)

    K = np.array((0, 0, 1))
    N_vec = np.cross(K, h_vec)
    N = np.linalg.norm(N_vec)

    if N_vec[0][1] >= 0:
        Omega_RAAN = np.arccos( N_vec[0][0] / N )
    else:
        Omega_RAAN = 2 * np.pi - np.arccos( N_vec[0][0] / N )

    # [section-5]
    r_norm = np.linalg.norm( r )
    e_vec = np.cross( v, h_vec ) / mu - r / r_norm
    eccentricity = np.linalg.norm( e_vec )

    if e_vec[0][2] >= 0:
        omega_argument_perigee = np.arccos( np.dot( N_vec[0], e_vec[0] ) / ( N * eccentricity ) )
    else:
        omega_argument_perigee = 2 * np.pi - np.arccos( np.dot( N_vec[0], e_vec[0] ) / ( N * eccentricity ) )
        
    # true anomaly
    v_r = np.dot( v[0], r[0] ) / np.linalg.norm( r )
    if v_r >= 0:
        true_anomaly = np.arccos( np.dot( r[0] / r_norm, e_vec[0] / eccentricity) )
    else:
        true_anomaly = 2 * np.pi - np.arccos( np.dot( r[0] / r_norm, e_vec[0] / eccentricity ) )

    # a = - mu / (2 * energy)
    # v^2/mu - 1/r
    v_norm = np.linalg.norm( v )
    energy = v_norm ** 2 / 2 - mu / r_norm
    semi_major_axis = -mu / ( 2 * energy )

    return Omega_RAAN, inclination, omega_argument_perigee, eccentricity, semi_major_axis, true_anomaly


def period( mu, semi_major_axis ):

    return 2 * math.pi * math.sqrt( semi_major_axis ** 3 / mu )