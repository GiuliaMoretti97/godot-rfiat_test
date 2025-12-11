# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
import numpy as np
from godot import cosmos, core
from godot.core import tempo

class EpochXYZ:
    def __init__(self, epoch: core.tempo.Epoch, x: float, y: float, z: float):
        self.epoch = epoch
        self.x = x
        self.y = y
        self.z = z
        
def to_rad(degrees : float):
    return degrees * np.pi / 180.0
def to_deg(radiants : float):
    return radiants * 180.0 / np.pi

def KeplerianPropagatorJ2(
    uni: cosmos.Universe,
    center: str,
    axes: str,
    ephoch0: tempo.Epoch,
    duration: float,
    dt: float,
    semimajoraxis: float,
    eccentricity: float,
    inclination: float,
    mean_anomaly_0: float = 0.0,
    raan_0: float = 0.0,
    arg_perigee_0: float = 0.0,
    tol: float = 1.0e-9,
) -> list[EpochXYZ]:
    
    # Convert angles to rad
    mean_anomaly_0 = to_rad(mean_anomaly_0)
    inclination = to_rad(inclination)
    raan_0 = to_rad(raan_0)
    arg_perigee_0 = to_rad(arg_perigee_0)
    
    # constants
    GM = uni.constants.get(center + "GM")
    Radius = uni.constants.get(center + "Radius")
    J2 = uni.constants.get(center + "J2")
    # mean motion
    n = np.sqrt(GM / (semimajoraxis ** 3))

    # auxiliary variables
    aux1 = (1 - eccentricity * eccentricity)
    A = J2 * n * ((Radius / semimajoraxis) ** 2.0) / aux1 ** 2.0
    B = J2 * n * (Radius / semimajoraxis) ** 2.0 / aux1 ** (3.0 / 2.0)
    
    # Mean Values
    true_anomaly_0 = core.astro.trueFromMean(mean_anomaly_0, eccentricity, tol)
    #! For this test the mean_n is not used so it is equal to n
    # aux = semimajoraxis **2 - 6.0 * J2 * Radius **2 * np.sin(inclination) **2 * np.cos( 2.0 * arg_perigee_0 + 2 * true_anomaly_0)
    #mean_semimajoraxis = (semimajoraxis + np.sqrt(aux)) / 2.0
    #mean_n = np.sqrt(GM / (mean_semimajoraxis**3))
    mean_n = n
    
    # [X, Y, Z]
    def computeCoordinates(dt: float):
        # RAAN
        raan = raan_0 - 3.0 / 2.0 * A * np.cos(inclination) * dt
        # omega
        arg_perigee = (
            arg_perigee_0
            + 3.0 / 4.0 * A * (5.0 * np.cos(inclination) * np.cos(inclination) - 1.0) * dt
        )
        # mean anomaly
        mean_anomaly = (
            mean_anomaly_0
            + mean_n * dt
            + 3.0 / 4.0 * B * (3.0 * np.cos(inclination) * np.cos(inclination) - 1) * dt
        )
        #print("raan ", raan)
        #print("arg_perigee ", arg_perigee)
        #print("mean_anomaly ", mean_anomaly)
        #mean_anomaly = mean_anomaly-int(mean_anomaly/(2.0*np.pi))*(2.0*np.pi)
        # True Anomaly
        true_anomaly = core.astro.trueFromMean(mean_anomaly, eccentricity, tol)  
        #print("ecc_anomaly ", ecc_anomaly)
        #print("true_anomaly ", true_anomaly)
        # compute Coordinates XYZ    
        # coe = [a(Km), e, I(rad), RAAN(rad), argP(rad), ta(rad)]
        coe = [semimajoraxis, eccentricity, inclination, raan, arg_perigee, true_anomaly]
        return core.astro.cartFromKep(coe, GM)
    
    # Propagate in ECI
    results = []
    for t in np.arange(0.0, duration + dt, dt):
        xyz = computeCoordinates(t)
        epoch = ephoch0 + t
        results.append(EpochXYZ(epoch, xyz[0], xyz[1], xyz[2]))
    return results
    

    
    
if __name__ == "__main__":
    # optionally avoid verbose logging messages
    cosmos.util.suppressLogger()
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    
    semimajoraxis = 7078 
    eccentricity = 0.0
    inclination = 98.0
    mean_anomaly_0=0.0
    raan_0 = 0.0
    arg_perigee_0 = 0.0

    propagated_orbit = KeplerianPropagatorJ2(
        uni=uni,
        center="Moon",
        axes="EME2000",
        ephoch0=core.tempo.parseEpoch("2025-12-31T23:43:53.273 UTC"),
        duration=10.0 * 24 * 3600.0,
        dt=60.0,
        semimajoraxis=semimajoraxis,
        eccentricity=eccentricity,
        inclination=inclination,
        mean_anomaly_0=mean_anomaly_0,
        raan_0=raan_0,
        arg_perigee_0=arg_perigee_0,
    )
    


