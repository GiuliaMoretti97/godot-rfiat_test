# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
from godot import core, cosmos
from Lander import addStationtoUniverse
from Propagator import KeplerianPropagatorJ2, to_rad, to_deg
import numpy as np
from Lander import check_occultation as check_moon_occultation
# Input Set to be used for the RFIAT-TEST

tj2000 = core.tempo.Epoch("2000-01-01T12:00:00.000 UTC") 
EME2000 = "EME2000"
ME = "Moon"
J2000_EPOCH = 2451545.0

if __name__ == "__main__":
    
    # INPUT ####################################################################
    # Initial Time Tag for the simulation
    t0 = core.tempo.Epoch("2025-12-31T23:43:53.273 UTC")
    t0 = core.tempo.convert(core.tempo.TimeScale.TDB, t0)
    Duration = 1 * 24 * 36000.0  # s
    dT = 30.0  # s
    # Initial coe = [a(Km), e, I(rad), RAAN(rad), argP(rad), ta(rad)]
    sma = 2000.0 # km Semi-major axes
    ecc = 0.1  # Eccentricity
    inc = 30.0  # Inclination
    ran = 0.0 # Right ascension of the ascending node
    aop = 0.0  # Argument of pericentre
    tan = 0.0  # True anomaly at reference epoch
    # Laner Time Series : 30 min/1 s
    lander_Duration = Duration  # s
    lander_dt = dT
    
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    
    # Radius
    MoonRadius = uni.constants.getRadius("Moon")
    EarthRadius = uni.constants.getRadius("Earth")
    
    
    # ###########################################################################
    ma0 = to_deg(core.astro.meanFromTrue(to_rad(tan), ecc))
    # Propagate Main Body (Moon) Centered in EME2000
    moon_orbit_eme2000 = KeplerianPropagatorJ2(
        uni=uni,
        center="Moon",
        axes=EME2000,
        ephoch0=t0,
        duration=Duration,
        dt=dT,
        semimajoraxis=sma,
        eccentricity=ecc,
        inclination=inc,
        mean_anomaly_0=ma0,
        raan_0=ran,
        arg_perigee_0=aop,
    )
    print("Orbiter Propagated, Moon Centered in EME2000")
    
    # ###########################################################################

    ### add Moon Lander to universe
    addStationtoUniverse(
        uni = uni,
        name = "MoonLander",
        center = "Moon",
        axis = ME,
        latitude = 0.0,
        longitude = 0.0,
        altitude = 0.0
    )
    ### add Earth Lander to universe
    addStationtoUniverse(
        uni = uni,
        name = "EarthLander",
        center = "Earth", 
        axis = "ITRF",
        latitude = 0,
        longitude = 0,
        altitude = 0.010
    )
    
    tL = t0 + lander_Duration
    epochs_lander = core.tempo.EpochRange(t0, tL).createGrid(lander_dt)
    tot = len(epochs_lander)
    
    
    moon_godot_1 = 0
    moon_godot_2 = 0
    moon_rfiat = 0
    earth_godot_1 = 0
    earth_godot_2 = 0
    earth_rfiat = 0
    with open("occultations_godot_rfiat.txt", "w") as f:
        f.write("MoonLanderGodot MoonLanderRFIAT EarthLanderGodot_Moon EarthLanderRFIAT_Moon EartLanderGodot_Earth EarthLanderRFIAT_Earth EarthLanderGodot EarthLanderRFIAT\n")
        for i, timetag in enumerate(epochs_lander):
            # MoonLander position wrt Earth
            moonlander = uni.frames.vector3("Earth", "MoonLander", EME2000, timetag)
            # EarthLander
            earthlander = uni.frames.vector3("Earth", "EarthLander", EME2000, timetag)
            # Moon Ephemeridis in EME2000 Earth centered 
            moon = uni.frames.vector3("Earth", "Moon", EME2000, timetag)
            # SC position wrt Earth
            spacecraft = moon + [
                moon_orbit_eme2000[i].x,
                moon_orbit_eme2000[i].y,
                moon_orbit_eme2000[i].z,
            ]
            # Position of the MoonLander wrt Moon (observer_to_occulter)
            moonlander_to_moon = uni.frames.vector3("MoonLander", "Moon", "EME2000", timetag)
            # SC to Moon 
            spacecraft_to_moon = moon - spacecraft
            # MoonLander position wrt SC (observer_to_spacecraft)
            moonlander_to_spacecraft = np.array(spacecraft) - np.array(moonlander)
            # Earth Lander occulted by Moon
            earthlander_to_moon = uni.frames.vector3("EarthLander", "Moon", "EME2000", timetag)
            # Earth Lander occulted wrt Earth
            earthlander_to_earth = uni.frames.vector3("EarthLander", "Earth", EME2000, timetag)
            earthlander_to_spacecraft = np.array(spacecraft) - np.array(earthlander)
            spacecraft_to_earthlander = - earthlander_to_spacecraft
            spacecraft_to_earth = - np.array(spacecraft)
            
            # Occultation with respect Moon (1)
            occ_coef_moon_1 = core.astro.computeOccultationCoefficient(
                earthlander_to_moon,
                earthlander_to_spacecraft,
                MoonRadius,
                0.0,
            )
            occ_moon_1 = core.astro.computeOccultationMargin(occ_coef_moon_1).partial >= 0.0
            # Occultation with respect Moon (1)
            occ_coef_moon_2 = core.astro.computeOccultationCoefficient(
                spacecraft_to_moon,
                spacecraft_to_earthlander,
                MoonRadius,
                0.0,
            )
            occ_moon_2 = (
                core.astro.computeOccultationMargin(occ_coef_moon_2).total >= 0.0
            )
            # Occultation computed with RFIAT function / Moon
            occ_moon_rfiat = check_moon_occultation(
                earthlander,
                spacecraft,
                moon,
                MoonRadius,
            )
            
            # Occultation with respect Earth (2)
            occ_coeff_earth_1 = core.astro.computeOccultationCoefficient(
                earthlander_to_earth,
                earthlander_to_spacecraft,
                EarthRadius,
                0.0,
            )
            occ_earth_1 = (
                core.astro.computeOccultationMargin(occ_coeff_earth_1).total >= 0.0
            )
            # Occultation with respect Earth (2)
            occ_coeff_earth_2 = core.astro.computeOccultationCoefficient(
                spacecraft_to_earth,
                spacecraft_to_earthlander,
                EarthRadius,
                0.0,
            )
            occ_earth_2 = (
                core.astro.computeOccultationMargin(occ_coeff_earth_2).total >= 0.0
            )
            # RFIAT / Earth
            occ_earth_rfiat = check_moon_occultation(
                earthlander,
                spacecraft,
                [0,0,0],
                EarthRadius,
            )
            
            if (occ_earth_1):
                earth_godot_1 += 1
            if (occ_earth_2):
                earth_godot_2 += 1
            if (occ_earth_rfiat):
                earth_rfiat +=1
            if (occ_moon_1):
                moon_godot_1 +=1
            if (occ_moon_2):
                moon_godot_2 +=1
            if (occ_moon_rfiat):
                moon_rfiat +=1            
            
    # print results

    print("EarthLander/Moon Godot (Lander OBS): ", 100 * moon_godot_1 / tot)
    print("EarthLander/Moon Godot (SC OBS): ", 100 * moon_godot_2 / tot)
    print("EarthLander/Moon RFIAT: ", 100 * moon_rfiat / tot)
    print("EarthLander/Earth Godot (Lander OBS): ", 100 * earth_godot_1 / tot)
    print("EarthLander/Earth Godot (SC OBS): ", 100 * earth_godot_2 / tot)
    print("EarthLander/Earth RFIAT: ", 100 * earth_rfiat / tot)
    



        
    