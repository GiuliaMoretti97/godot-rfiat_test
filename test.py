import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions.conversions import Rx, Ry, Rz

if __name__== "__main__":
        # Define angles in radians 
        alpha = np.pi / 2  # Z rotation
        beta = 32 * np.pi / 180  # Y rotation 
        gamma = 11 * np.pi / 180  # X rotation 
        # Get individual matrices
        R3 = Rz(alpha)
        R2 = Ry(beta)
        R1 = Rx(gamma)
        # Order: Rz * Ry * Rx
        Rot_A_to_B = R3 @ R2 @ R1
        tmp1 = R2 @ R1
        Rot_tmp = R3 @ tmp1
        print("Matrix 1 ", np.round(Rot_A_to_B, 2))
        print("Matrix 2 ", np.round(Rot_tmp, 2))

