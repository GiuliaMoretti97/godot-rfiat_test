# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
from godot import core, cosmos
from Lander import addStationtoUniverse
from Propagator import EpochXYZ, KeplerianPropagatorJ2, to_rad, to_deg
from MoonRotations import get_rotation_matrix_ME_to_EME
import numpy as np
# Input Set to be used for the RFIAT-TEST

tj2000 = core.tempo.Epoch("2000-01-01T12:00:00.000 UTC") 
EME2000 = "EME2000"
ME = "Moon"
J2000_EPOCH = 2451545.0

def check_moon_occultation(rx_pos, tx_pos, moon_pos, moon_radius):
    """
    Verifica l'occultazione della Luna.
    Ritorna True se NON c'è occultazione (visibile), False se è occultato.
    Args:
        rx_pos (array-like): Posizione Ricevitore [x, y, z]
        tx_pos (array-like): Posizione Trasmettitore [x, y, z]
        moon_pos (array-like): Posizione Luna [x, y, z]
        moon_radius (float): Raggio della Luna
    """
    # Conversione in array numpy per calcoli vettoriali
    rx_pos = np.array(rx_pos, dtype=float)
    tx_pos = np.array(tx_pos, dtype=float)
    moon_pos = np.array(moon_pos, dtype=float)
    # 1. Calcolo vettori relativi al centro della Luna (Moon Centered Frame)
    # VB: MathModule.VecDiff(Receiver.Position, getMoonPosition_ECEF)
    rx_moon = rx_pos - moon_pos
    tx_moon = tx_pos - moon_pos
    # 2. Calcolo delle distanze (Norme)
    # VB: MathModule.GetNorm(...)
    range_mc_rx = np.linalg.norm(rx_moon)
    range_mc_tx = np.linalg.norm(tx_moon)
    # Check di sicurezza come nel codice originale (evita errori di dominio nella radice quadrata)
    if range_mc_rx < moon_radius:
        range_mc_rx = moon_radius
    if range_mc_tx < moon_radius:
        range_mc_tx = moon_radius
    # 3. Distanza all'orizzonte (Grazing distance)
    # Teorema di Pitagora: sqrt(ipotenusa^2 - cateto^2)
    m_graze_dist_rx = np.sqrt(range_mc_rx**2 - moon_radius**2)
    m_graze_dist_tx = np.sqrt(range_mc_tx**2 - moon_radius**2)
    # 4. Calcolo vettore e distanza tra TX e RX
    # Nota: Nel codice VB 'RangeTXRX' è usato prima di essere definito esplicitamente,
    # ma logicamente è la distanza tra i due punti.
    tx_rx_vec = rx_pos - tx_pos  # Vettore da TX a RX
    range_tx_rx = np.linalg.norm(tx_rx_vec)  # Distanza scalare
    # 5. Controlli preliminari (Early Exits)
    # Se la distanza tra i punti è inferiore alla distanza dell'orizzonte,
    # la curvatura della luna non può bloccarli.
    if range_tx_rx <= m_graze_dist_rx:
        return True
    if range_tx_rx <= m_graze_dist_tx:
        return True
    # 6. Calcolo Angoli di "Grazing" (Angolo limite)
    # asin(Opposto / Ipotenusa)
    m_graze_angle_rx = np.arcsin(moon_radius / range_mc_rx)
    m_graze_angle_tx = np.arcsin(moon_radius / range_mc_tx)

    # 7. Calcolo Angoli tra vettori
    # Funzione helper interna per calcolare l'angolo tra due vettori
    def get_angle(v1, v2):
        unit_v1 = v1 / np.linalg.norm(v1)
        unit_v2 = v2 / np.linalg.norm(v2)
        dot_product = np.dot(unit_v1, unit_v2)
        # Clip per evitare errori numerici fuori da [-1, 1]
        dot_product = np.clip(dot_product, -1.0, 1.0)
        return np.arccos(dot_product)

    # AngleMCRXTX: Angolo tra (RX-Moon) e (RX-TX)
    angle_mc_rx_tx = get_angle(rx_moon, tx_rx_vec)
    # AngleECTXRX: PiGreco - Angolo tra (RX-TX) e (TX-Moon)
    # Nota: nel VB usa VecDiff(Receiver, Transmitter) quindi il vettore punta verso RX.
    angle_tx_rx_v_tx_moon = get_angle(tx_rx_vec, tx_moon)
    angle_ec_tx_rx = np.pi - angle_tx_rx_v_tx_moon
    # 8. Verifiche Finali
    # Se l'angolo reale è maggiore dell'angolo di grazing, siamo "sopra" l'ostacolo
    if angle_mc_rx_tx >= m_graze_angle_rx:
        return True
    if angle_ec_tx_rx >= m_graze_angle_tx:
        return True
    # Se nessuno dei check passa, c'è occultazione
    return False

        
if __name__ == "__main__":
    
    # INPUT ####################################################################
    # Initial Time Tag for the simulation
    t0 = core.tempo.Epoch("2025-12-31T23:43:53.273 UTC")
    t0 = core.tempo.convert(core.tempo.TimeScale.TDB, t0)
    Duration = 1 * 24 * 36000.0  # s
    dT = 30.0  # s
    # Spacecraft Name
    sc_name = "SC"
    # Main Body 
    center = "Moon"
    gravity = "MoonGravityOnly"
    # Initial coe = [a(Km), e, I(rad), RAAN(rad), argP(rad), ta(rad)]
    sma = 2000.0 # km Semi-major axes
    ecc = 0.1  # Eccentricity
    inc = 30.0  # Inclination
    ran = 0.0 # Right ascension of the ascending node
    aop = 0.0  # Argument of pericentre
    tan = 0.0  # True anomaly at reference epoch
    # Lander Position
    lander_latitude = 0.0
    lander_longitude = 0.0
    lander_altitude = 0.0
    # Laner Time Series : 30 min/1 s
    lander_Duration = Duration  # s
    lander_dt = dT
    
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    
    # ###########################################################################
    ma0 = to_deg(core.astro.meanFromTrue(to_rad(tan), ecc))
    # Propagate Main Body (Moon) Centered in EME2000
    moon_orbit_eme2000 = KeplerianPropagatorJ2(
        uni=uni,
        center=center,
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
    # Compute Earth Centered EME2000
    earth_orbit_eme2000 : list[EpochXYZ] = []
    for i, tt in enumerate(moon_orbit_eme2000):
        earth_moon_eme2000_xyz = uni.frames.vector3("Earth", center, EME2000, tt.epoch)
        earth_orbit_eme2000_xyz = earth_moon_eme2000_xyz + [tt.x, tt.y, tt.z]
        earth_orbit_eme2000.append(EpochXYZ(tt.epoch, earth_orbit_eme2000_xyz[0], earth_orbit_eme2000_xyz[1], earth_orbit_eme2000_xyz[2]))
    
    # Write Coordinates in Moon EME2000
    with open("xyz_moon_eme2000.txt", "w") as f:
        for i, tt in enumerate(moon_orbit_eme2000):
            f.write(
                f"{tt.x} {tt.y} {tt.z} \n"
            )
    # Wirte Coordinates in Earth EME2000
    with open("xyz_earth_eme2000.txt", "w") as f:
        for i, tt in enumerate(earth_orbit_eme2000):
            f.write(f"{tt.x} {tt.y} {tt.z} \n")
    
    print("Propagation Files written")
    # ###########################################################################

    ### add Moon Lander to universe
    addStationtoUniverse(
        uni = uni,
        name = "Lander",
        center = center,
        axis = ME,
        latitude = lander_latitude,
        longitude = lander_longitude,
        altitude = lander_altitude
    )
    ### add Earth Lander to universe
    addStationtoUniverse(
        uni = uni,
        name = "EarthLander",
        center = center,
        axis = "ITRF",
        latitude = 0,
        longitude = 0,
        altitude = 0
    )
    
    tL = t0 + lander_Duration
    epochs_lander = core.tempo.EpochRange(t0, tL).createGrid(lander_dt)
    
    moon_lander_eme2000 : list[EpochXYZ] = []
    earth_lander_eme2000 : list[EpochXYZ] = []
    earth_earthlander_eme2000 : list[EpochXYZ] = []
    with open("occultations_godot_rfiat.txt", "w") as f:
        f.write("TimeTag MoonLanderGodot MoonLanderRFIAT EarthLanderGodot_Moon EarthLanderRFIAT_Moon EartLanderGodot_Earth EarthLanderRFIAT_Earth EarthLanderGodot EarthLanderRFIAT\n")
        for i, timetag in enumerate(epochs_lander):
            # Epochs vatables
            epoch_jd = timetag.jd(
                core.tempo.TimeScale.TDB, jdType=core.tempo.JulianDay.Type.JD
            )
            d = epoch_jd - J2000_EPOCH
            T = d / 36525.0
            # get Moon/Lander position in EME2000
            moon_lander_xyz = uni.frames.vector3(center, "Lander", EME2000, timetag)
            # get Earth/Lander position in EME2000
            earth_lander_xyz = uni.frames.vector3("Earth", "Lander", EME2000, timetag)
            # get Lander position in ME Frame (moon Centered)
            me_xyz = uni.frames.vector3(
                center,
                "Lander",
                ME,
                timetag
            )
            # Moon Ephemeridis in EME2000 Earth centered 
            earth_moon_eme2000_xyz = uni.frames.vector3("Earth", center, EME2000, timetag)
            # RFIAT rotation matrix computation VS Godot
            godot_me_2_eme = uni.frames.rotationMatrix("Moon", "EME2000", 0, epoch=timetag)
            rfiat_me_2_eme = get_rotation_matrix_ME_to_EME(T, d)
            rfiat_eme_xyz = rfiat_me_2_eme @ me_xyz
            moon_lander_eme2000.append(
                EpochXYZ(
                    timetag, moon_lander_xyz[0], moon_lander_xyz[1], moon_lander_xyz[2]
                )
            )
            earth_lander_eme2000.append(
                EpochXYZ(
                    timetag, earth_lander_xyz[0], earth_lander_xyz[1], earth_lander_xyz[2]
                )
            )
            # Position of the MoonLander wrt Moon (observer_to_occulter)
            moonlander_to_moon_occ = uni.frames.vector3(
                "Lander", "Moon", "EME2000", timetag
            )
            # SC position wrt Earth
            spacecraft = [
                earth_orbit_eme2000[i].x,
                earth_orbit_eme2000[i].y,
                earth_orbit_eme2000[i].z,
            ]
            # MoonLander position wrt Earth
            moonlander = [earth_lander_eme2000[i].x, earth_lander_eme2000[i].y, earth_lander_eme2000[i].z]
            # MoonLander position wrt SC (observer_to_spacecraft)
            moonlander_to_spacecraft = np.array(spacecraft) - np.array(moonlander)
            # Compute occultation coefficient for SC
            occ = core.astro.computeOccultationCoefficient(
                moonlander_to_moon_occ,
                moonlander_to_spacecraft,
                uni.constants.getRadius("Moon"),
                0.0,
            )
            occ_margin = core.astro.computeOccultationMargin(occ).total
            occ_godot = True
            if occ_margin < 0.0:
                occ_godot = False
            occ_rfiat = check_moon_occultation(
                moonlander,
                spacecraft,
                earth_moon_eme2000_xyz,
                uni.constants.getRadius("Moon"),
            )
            if (occ_godot is not occ_rfiat):
                print("Moon Lander - Discrepancy found at epoch ", timetag)

            # Earth Lander occulted by Moon
            earthlander_to_moon_occ = uni.frames.vector3(
                "EarthLander", "Moon", "EME2000", timetag
            )
            # Earth Lander occulted wrt Earth
            earthlander_to_earth_occ = uni.frames.vector3(
                "EarthLander", "Earth", "EME2000", timetag
            )
            earthlander = - earthlander_to_earth_occ
            earth_earthlander_eme2000.append(
                EpochXYZ(
                    timetag, earthlander[0], earthlander[1], earthlander[2]
                )
            )
            earthlander_to_spacecraft = np.array(spacecraft) - np.array(earthlander)
            # Occultation with respect Moon (1)
            occ_1 = core.astro.computeOccultationCoefficient(
                earthlander_to_moon_occ,
                earthlander_to_spacecraft,
                uni.constants.getRadius("Moon"),
                0.0,
            )
            occ_tf_1 = core.astro.computeOccultationMargin(occ_1).total >= 0.0
            # Occultation with respect Earth (2)
            occ_2 = core.astro.computeOccultationCoefficient(
                earthlander_to_earth_occ,
                earthlander_to_spacecraft,
                uni.constants.getRadius("Earth"),
                0.0,
            )
            occ_tf_2 = core.astro.computeOccultationMargin(occ_2).total >= 0.0
            occ_earthlander = occ_tf_1 and occ_tf_2 
            # Occultation computed with RFIAT function
            occ_rfiat_1 = check_moon_occultation(
                earthlander,
                spacecraft,
                earth_moon_eme2000_xyz,
                uni.constants.getRadius("Moon"),
            )
            occ_rfiat_2 = check_moon_occultation(
                earthlander,
                spacecraft,
                [0, 0, 0],
                uni.constants.getRadius("Earth"),
            )
            occ_earthlander_rfiat = occ_rfiat_1 and occ_rfiat_2
            if (occ_earthlander is not occ_earthlander_rfiat):
                print("Discrepancy found at epoch ", timetag)
            f.write(
                f"{timetag} {occ_godot} {occ_rfiat} {occ_tf_1} {occ_rfiat_1} {occ_tf_2} {occ_rfiat_2} {occ_earthlander} {occ_earthlander_rfiat} \n"
            )
    print("Lander Propagated, Moon Centered in EME2000")
        
    # Files
    with open("xyz_moon_lander_eme2000.txt", "w") as f:
        for i, tt in enumerate(moon_lander_eme2000):
            f.write(f"{tt.x} {tt.y} {tt.z} \n")
    with open("xyz_earth_lander_eme2000.txt", "w") as f:
        for i, tt in enumerate(earth_lander_eme2000):
            f.write(f"{tt.x} {tt.y} {tt.z} \n")
    with open("xyz_earth_earthlander_eme2000.txt", "w") as f:
        for i, tt in enumerate(earth_earthlander_eme2000):
            f.write(f"{tt.x} {tt.y} {tt.z} \n")
    print("Lander Files written")
    
    tj2000_tdb = core.tempo.convert(core.tempo.TimeScale.TDB, tj2000)
    tj2000_jd = tj2000_tdb.jd(core.tempo.TimeScale.TDB, jdType=core.tempo.JulianDay.Type.JD)
    tj2000_mjd = tj2000_tdb.jd(
        core.tempo.TimeScale.TDB, jdType=core.tempo.JulianDay.Type.MJD
    )

        
    