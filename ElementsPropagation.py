
# based on https://orbital-mechanics.space/time-since-periapsis-and-keplers-equation/elliptical-orbit-example.html
import numpy as np
from scipy.optimize import newton
import math

def kepler(E, M_e, e):
    """Kepler's equation, to be used in a Newton solver."""
    return E - e * np.sin(E) - M_e

 
def d_kepler_d_E(E, M_e, e):
    """The derivative of Kepler's equation, to be used in a Newton solver.

    Note that the argument M_e is unused, but must be present so the function
    arguments are consistent with the kepler function.
    """
    return 1 - e * np.cos(E)


def KeplerM2E( M, e ):

    E = newton(func=kepler, fprime=d_kepler_d_E, x0=np.pi, args=(M, e), tol=1e-10)
    
    return E


def KeplerE2M( E, e ):

    return E - e * np.sin(E)


def PropagateTrueAnomaly( true_anomaly, e, T_period, delta_time ):

    E = math.sqrt( 1 - e ** 2 ) * math.sin( true_anomaly ) / ( e + math.cos( true_anomaly ) ) 

    M = KeplerE2M( E = E, e = e )
    
    M_new = M + 2 * np.pi / T_period * delta_time

    E_new = KeplerM2E( M= M_new, e = e )

    true_anomaly_new = 2 * math.atan( math.sqrt( 1 + e / ( 1 - e ) ) * math.tan( E_new / 2 ) )

    return true_anomaly_new