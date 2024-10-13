from math import cos, sin, sqrt

def Omega_RAAN_derivative(
    J2,
    R_Earth,
    mu_earth,
    a,
    e,
    i
):
    return -3/2 * J2 * ( R_Earth / (a * (1 - e ** 2 ) ) ) ** 2 * sqrt( mu_earth / a ** 3 ) * cos( i )

def omega_argument_perigee_derivative(
    J2,
    R_Earth,
    mu_earth,
    a,
    e,
    i
):
    return -3 / 2 * J2 * ( R_Earth / ( a * ( 1 - e ** 2 ) ) ) ** 2 * sqrt( mu_earth / a ** 3 ) * ( 5 * cos( i ) ** 2 - 1 )