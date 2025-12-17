import numpy as np
from functions.conversions import Rx, Ry, Rz


def matrix_PA2ME():
        rot_x = Rx( - (0.2785 / 3600.0) * np.pi / 180.0)
        rot_y = Ry( - (78.6944 / 3600.0) * np.pi / 180.0)
        rot_z = Rz( - (67.8526 / 3600.0) * np.pi / 180.0)
        return rot_x @ (rot_y @ rot_z)

def calculate_astronomic_values(T, d):
    """
    Calculates the mean right ascension (alpha_0), mean declination (delta_0),
    and Greenwich Mean Sidereal Time (W) based on the provided trigonometric series.

    The equations are based on the input image and assume T and d are
    input in units compatible with the coefficients (likely T in Julian centuries
    since J2000.0 and d in days).

    Args:
        T (float): Time variable, e.g., in Julian centuries.
        d (float): Time variable, e.g., in days.

    Returns:
        dict: A dictionary containing the calculated values for
              'alpha_0', 'delta_0', and 'W'.
    """

    # --- 1. Calculate the arguments E1 through E13 (in degrees) ---
    # The coefficients for E1-E13 seem to relate to the motions of specific
    # celestial bodies or nodes.

    # Note: The original image has '0.0529921d', which is interpreted as 0.0529921 * d,
    # and the units are assumed to be in degrees.
    E1 = 125.045 - 0.0529921 * d
    E2 = 250.089 - 0.1059842 * d
    E3 = 260.008 + 13.0120009 * d
    E4 = 176.625 + 13.3407154 * d
    E5 = 357.529 + 0.9856003 * d
    E6 = 311.589 + 26.4057084 * d
    E7 = 134.963 + 13.0649930 * d
    E8 = 276.617 + 0.3287146 * d
    E9 = 34.226 + 1.7484877 * d
    E10 = 15.134 - 0.1589763 * d
    E11 = 119.743 + 0.0036096 * d
    E12 = 239.961 + 0.1643573 * d
    E13 = 25.053 + 12.9590088 * d

    # Store E values in a list for easier access and verification
    E_values_deg = [E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13]

    # Convert all E values from degrees to radians for NumPy's trigonometric functions
    E_values_rad = np.deg2rad(E_values_deg)

    # --- 2. Calculate Mean Right Ascension ($\alpha_0$) ---
    # The terms are: Constant + T term + sin(E) terms

    alpha_0 = (
        269.9949
        + 0.0031 * T
        - 3.8787 * np.sin(E_values_rad[0])  # E1
        - 0.1204 * np.sin(E_values_rad[1])  # E2
        + 0.0700 * np.sin(E_values_rad[2])  # E3
        - 0.0172 * np.sin(E_values_rad[3])  # E4
        + 0.0072 * np.sin(E_values_rad[5])  # E6 (E5 is skipped)
        - 0.0052 * np.sin(E_values_rad[9])  # E10
        + 0.0043 * np.sin(E_values_rad[12])  # E13
    )

    # --- 3. Calculate Mean Declination ($\delta_0$) ---
    # The terms are: Constant + T term + cos(E) terms

    delta_0 = (
        66.5392
        + 0.0130 * T
        + 1.5419 * np.cos(E_values_rad[0])  # E1
        + 0.0239 * np.cos(E_values_rad[1])  # E2
        - 0.0278 * np.cos(E_values_rad[2])  # E3
        + 0.0068 * np.cos(E_values_rad[3])  # E4
        - 0.0029 * np.cos(E_values_rad[5])  # E6 (E5 is skipped)
        + 0.0009 * np.cos(E_values_rad[6])  # E7
        + 0.0008 * np.cos(E_values_rad[9])  # E10
        - 0.0009 * np.cos(E_values_rad[12])  # E13
    )

    # --- 4. Calculate Greenwich Mean Sidereal Time (W) ---
    # The terms are: Constant + d term + d^2 term + sin(E) terms
    # Note: 1.4 \times 10^{-12} is represented as 1.4e-12 in Python.

    W = (
        38.3213
        + 13.17635815 * d
        - 1.4e-12 * (d**2)
        + 3.5610 * np.sin(E_values_rad[0])  # E1
        + 0.1208 * np.sin(E_values_rad[1])  # E2
        - 0.0642 * np.sin(E_values_rad[2])  # E3
        + 0.0158 * np.sin(E_values_rad[3])  # E4
        + 0.0252 * np.sin(E_values_rad[4])  # E5
        - 0.0066 * np.sin(E_values_rad[5])  # E6
        - 0.0047 * np.sin(E_values_rad[6])  # E7
        - 0.0046 * np.sin(E_values_rad[7])  # E8
        + 0.0028 * np.sin(E_values_rad[8])  # E9
        + 0.0052 * np.sin(E_values_rad[9])  # E10
        + 0.0040 * np.sin(E_values_rad[10])  # E11
        + 0.0019 * np.sin(E_values_rad[11])  # E12
        - 0.0044 * np.sin(E_values_rad[12])  # E13
    )
    W = W % 360.0
    phi = alpha_0 + 90.0
    theta = 90.0 - delta_0
    psi = W
    return [phi, theta, psi]
    

def get_rotation_matrix_EME_to_ME(T, d):
    """
    Returns the rotation matrix from Lunar Mean Equator (LME2000)
    to Mean Earth (ME) based on angle W.

    Matrix R =
    [ 0.998...   0.049...  -0.022... ]
    [-0.054...   0.909...  -0.412... ]
    [ 0          0.412...   0.910... ]

    Note: The sign of R[1][2] (-0.412...) is inferred as negative to ensure
    the matrix is orthogonal (determinant = 1).
    """
    [phi, theta, psi] = calculate_astronomic_values(T, d)
    rot_z_psi = Rz(psi * np.pi / 180.0)
    rot_x_theta = Rx(theta * np.pi / 180.0)
    rot_z_phi = Rz(phi * np.pi / 180.0)
    tmp = rot_x_theta @ rot_z_phi
    rot_EME2PA = rot_z_psi @ tmp
    #! It is in the paper but formulation does not fit
    # = matrix_PA2ME()
    #rot_EME2ME = rot_PA2ME @ rot_EME2PA
    rot_EME2ME =  rot_EME2PA
    return rot_EME2ME

def get_rotation_matrix_ME_to_EME(T, d):
    """
    Returns the rotation matrix from Lunar Mean Equator (LME2000)
    to Mean Earth (ME) based on angle W.

    Matrix R =
    [ 0.998...   0.049...  -0.022... ]
    [-0.054...   0.909...  -0.412... ]
    [ 0          0.412...   0.910... ]

    Note: The sign of R[1][2] (-0.412...) is inferred as negative to ensure
    the matrix is orthogonal (determinant = 1).
    """
    EME_to_ME = get_rotation_matrix_EME_to_ME(T, d)
    return EME_to_ME.T