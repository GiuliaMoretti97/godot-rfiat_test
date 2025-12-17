# pyright: reportRedeclaration=false
# pyright: reportAttributeAccessIssue=false
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from godot.core import astro
import numpy as np
from godot import cosmos

from models.station_models import MountingType, StationPlugin
import json
import tempfile
import os
import math

MoonGM = 398600.4418
MoonRadius = 1737.4
MoonJ2 = 0.0010826359
arc = 180.0 / np.pi
# Constants
J2000_EPOCH = 2451545.0
W_DOT = 13.17635815  # deg/day

longitude = 0.0 / arc
latitude = 0.0 / arc
altitude = 0.0 / arc

x = MoonRadius * (np.cos(latitude) * np.cos(longitude))
y = MoonRadius * (np.cos(latitude) * np.sin(longitude))
z = MoonRadius * (np.sin(latitude))

def check_occultation(rx_pos, tx_pos, occ_pos, occ_radius):
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
    moon_pos = np.array(occ_pos, dtype=float)
    # 1. Calcolo vettori relativi al centro della Luna (Moon Centered Frame)
    # VB: MathModule.VecDiff(Receiver.Position, getMoonPosition_ECEF)
    rx_moon = rx_pos - moon_pos
    tx_moon = tx_pos - moon_pos
    # 2. Calcolo delle distanze (Norme)
    # VB: MathModule.GetNorm(...)
    range_mc_rx = np.linalg.norm(rx_moon)
    range_mc_tx = np.linalg.norm(tx_moon)
    # Check di sicurezza come nel codice originale (evita errori di dominio nella radice quadrata)
    if range_mc_rx < occ_radius:
        range_mc_rx = occ_radius
    if range_mc_tx < occ_radius:
        range_mc_tx = occ_radius
    # 3. Distanza all'orizzonte (Grazing distance)
    # Teorema di Pitagora: sqrt(ipotenusa^2 - cateto^2)
    m_graze_dist_rx = np.sqrt(range_mc_rx**2 - occ_radius**2)
    m_graze_dist_tx = np.sqrt(range_mc_tx**2 - occ_radius**2)
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
    m_graze_angle_rx = np.arcsin(occ_radius / range_mc_rx)
    m_graze_angle_tx = np.arcsin(occ_radius / range_mc_tx)

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


def calculate_W(julian_day):
    """
    Calculates the value of W (in degrees) for a given Julian Day based on the
    provided orbital formulation.
    """
    # Calculate D (Days past epoch of J2000)
    D = julian_day - J2000_EPOCH
    print("D = ", D)
    # Calculate T (Julian centuries past epoch J2000)
    # Note: T is defined in the image but not explicitly used in the W equation
    # (unless the D^2 term is a substitution), but we calculate it for completeness.
    T = D / 36525.0

    # Base constant terms
    term_constant = 38.3213
    term_linear = W_DOT * D
    term_quadratic = -1.4e-12 * (D**2)

    # Dictionary of Sine Terms
    # Format: [Amplitude, E_base, E_factor]
    # where E_i = E_base + E_factor * D
    sine_terms_data = [
        # i=1
        (3.5610, 125.045, -0.0529921),
        # i=2
        (0.1208, 250.089, -0.1059842),
        # i=3
        (-0.0642, 260.008, 13.0120009),
        # i=4
        (0.0158, 176.625, 13.3407154),
        # i=5
        (0.0252, 357.529, 0.9856003),
        # i=6
        (-0.0066, 311.589, 26.4057084),
        # i=7
        (-0.0047, 134.963, 13.0649930),
        # i=8
        (-0.0046, 276.617, 0.3287146),
        # i=9
        (0.0028, 34.226, 1.7484877),
        # i=10
        (0.0052, 15.134, -0.1589763),
        # i=11
        (0.0040, 119.743, 0.0036096),
        # i=12
        (0.0019, 239.961, 0.1643573),
        # i=13
        (-0.0044, 25.053, 12.9590088),
    ]

    # Summation of sine terms
    sum_sine_components = 0.0

    for coeff, base, factor in sine_terms_data:
        # Calculate angle E in degrees
        E_deg = base + (factor * D)

        # Convert E to radians for math.sin()
        E_rad = math.radians(E_deg)

        # Add to sum
        sum_sine_components += coeff * math.sin(E_rad)

    # Final Calculation
    W = term_constant + term_linear + term_quadratic + sum_sine_components

    # Normalize W to 0-360 range (optional, but standard for angles)
    W_normalized = W % 360.0

    return W_normalized


def get_rotation_matrix_LME_to_EME():
    """
    Returns the rotation matrix from Lunar Mean Equator (LME2000)
    to Earth Mean Equator (EME2000).

    Matrix R =
    [ 0.998...   0.049...  -0.022... ]
    [-0.054...   0.909...  -0.412... ]
    [ 0          0.412...   0.910... ]

    Note: The sign of R[1][2] (-0.412...) is inferred as negative to ensure
    the matrix is orthogonal (determinant = 1).
    """
    return [
        [0.9984965052050879, 0.0499357293985326, -0.0226086714041825],
        [-0.0548154092680678, 0.9096101252380440, -0.4124510189026893],
        [0.0, 0.4124510189026893, 0.9109797785934293],
    ]


def apply_rotation(matrix, vector):
    """
    Multiplies a 3x3 matrix by a 3D vector (list or tuple of length 3).
    Returns the resulting rotated vector as a list [x, y, z].
    """
    x, y, z = vector

    result = [0.0, 0.0, 0.0]

    # Row 0
    result[0] = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z)
    # Row 1
    result[1] = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z)
    # Row 2
    result[2] = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z)

    return result


def get_rotation_matrix_ME_to_LME(W_degrees : float):
    """
    Returns the rotation matrix from Mean Earth (ME) to
    Lunar Mean Equator (LME2000) based on angle W.

    Matrix R =
    [ cos(W)   -sin(W)   0 ]
    [ sin(W)    cos(W)   0 ]
    [ 0         0        1 ]
    """
    # Convert W to radians for trigonometric functions
    W_rad = math.radians(W_degrees)
    cos_w = math.cos(W_rad)
    sin_w = math.sin(W_rad)

    return [[cos_w, -sin_w, 0.0], [sin_w, cos_w, 0.0], [0.0, 0.0, 1.0]]





"""""
Latitude : xxx deg xx ' xx '' N/S
Longitude : xxx deg xx ' xx '' E/W
Altitude : float [km]
"""

def create_station_template(
        center : str,
        axis: str,        
) -> StationPlugin:
    station_object = StationPlugin(
        refEpoch="2000-01-01T12:00:00.000 UTC",
        refCenter=center,
        refAxis=axis,
        mountingType=MountingType.AZ_EL.value,
        heightAboveAntennaFoot="0.0 m",
        dryOffset="0 mm",
        stationComplex="xx",
        stationAlias=["XX"],
        coordinates=["1 km", "1 km", "1 km"],
        plateMotion=["0 cm/year", "0 cm/year", "0 cm/year"],
        semiMajorAxis=None,
        inverseOfFlattening=None,
        antennaDiameter=None,
        mounts=None,
        rangeCalLongLoopCorrection=None,
        maxAzimuthRate=None,
        rangeCalMediumLoopCorrection=None,
        groupDelay=None,
        phaseDelay=None,
        mappings=None,
        receivers=None,
        rangeUnitConversion=None,
        horizonMaskFile=None,
    )
    return station_object

def createStationConfig(
    uni: cosmos.Universe,
    name: str,
    center : str,
    axis: str,
    latitude: float,
    longitude: float,
    altitude: float,
) -> dict:
    """
    Args:
        uni: cosmos.Universe object with constants
        name: station name
        latitude: xxx deg xx ' xx '' N/S
        longitude: xxx deg xx ' xx '' E/W
        altitude: altitude in km
    """
    try:
        lat_rad = latitude * np.pi / 180.0
        lon_rad = longitude * np.pi / 180.0
        geo = [lat_rad, lon_rad, altitude]
        # Earth parameters
        radius = uni.constants.getRadius(center)
        flattening = uni.constants.getFlattening(center)
        xyz = astro.cartesianFromGeodetic(geo, radius, flattening)
        coordinates_string = [
            str(xyz[0]) + " km",
            str(xyz[1]) + " km",
            str(xyz[2]) + " km",
        ]
        station_object = create_station_template(center, axis)
        station_object.coordinates = coordinates_string
        station_object.stationAlias = [name]
        station_config = station_object.model_dump(exclude_none=True)
        return station_config
    except Exception as e:
        raise ValueError(f"GODOT Error : Error in createStationConfig: {e}")

def addStationtoUniverse(
    uni: cosmos.Universe,
    center: str,
    axis: str,
    name: str,
    latitude: float,
    longitude: float,
    altitude: float,
) -> None:
    """
    Defined a cosmos.Universe with a filled core.constants.ConstantBook
    It loads into the uni a Station with "name" and input coordinates
    It generates a tmp json file of the Station dictionary
    It creates Point and Axes associated to the Station into the
    uni.frames , having Station name
    The Axes is the Topocentric (ENU) at the station coordinates
    Args:
        uni: GODOT cosmos.Universe object with constants
        name: station name
        latitude: xxx deg xx ' xx '' N/S
        longitude: xxx deg xx ' xx '' E/W
        altitude: altitude in km

    """
    try:
        station_config = createStationConfig(
            uni=uni,
            name=name,
            center=center,
            axis=axis,
            latitude=latitude,
            longitude=longitude,
            altitude=altitude,
        )
        station_dict = {name: station_config}

        # create tmp input file
        fd, custom_station_path = tempfile.mkstemp(suffix=".json")
        with os.fdopen(fd, "w") as tmp_file:
            json.dump(station_dict, tmp_file, indent=4)
        print(f"Ground station tmp file created : {custom_station_path}")

        # Update Universe
        uni.frames.addStations(custom_station_path, uni.constants)

        # romove tmp file
        os.remove(custom_station_path)
        print(
            f"Ground station saved into Universe and tmp file removed : {custom_station_path}"
        )
    except Exception as e:
        raise RuntimeError(f"GODOT Error : Error in addStationtoUniverse: {e}")
