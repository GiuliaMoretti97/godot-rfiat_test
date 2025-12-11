# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from godot import core, cosmos
from MoonRotations import get_rotation_matrix_ME_to_EME
# Constants
J2000_EPOCH = 2451545.0
    
if __name__ == "__main__":
    
    # J2000 + 1GG
    #j2000 = 2451545.0
    timetag = core.tempo.parseEpoch("2025-12-31T23:43:53.000 UTC")
    #timetag = core.tempo.parseEpoch("2000-01-01T12:00:00.000 UTC")
    timetag = core.tempo.convert(core.tempo.TimeScale.TDB, timetag)
    timetag_jd = timetag.jd(core.tempo.TimeScale.TDB, jdType=core.tempo.JulianDay.Type.JD)
    print("TimeTag ", timetag)
    print("TimeTag JD ", timetag_jd)
        
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
        
    print("\nAxes : ", uni.frames.listAxesNames())
    print("\nPoints : ", uni.frames.listPointNames())
    ############################################################################

    t0 = timetag
    Duration = 1 * 24 * 36000.0  # s
    dT = 30.0  # s
    tf = t0 + Duration
    epochs = core.tempo.EpochRange(t0, tf).createGrid(dT)
    matrix_diff = []
    time_points = []
    with open("matrix_diff.txt", "w") as f:
        for time in epochs:
            time_jd = time.jd(
                core.tempo.TimeScale.TDB, jdType=core.tempo.JulianDay.Type.JD
            )
            T = (time_jd - J2000_EPOCH) / 36525.0
            godot_me_2_eme = uni.frames.rotationMatrix("Moon", "EME2000", 0, epoch=time)
            T = (time_jd - J2000_EPOCH) / 36525.0
            d = time_jd - J2000_EPOCH
            matrix_me2eme_test = get_rotation_matrix_ME_to_EME(T, d)
            matrix_diff.append(matrix_me2eme_test - godot_me_2_eme)
            time_points.append(time_jd)
            f.write(f"{d} ")
            for row in matrix_diff[-1]:
                for elem in row:
                    f.write(f"{elem} ")
            f.write("\n")
        
    # Convert to numpy array for easier indexing
    matrix_diff = np.array(matrix_diff)
    
    # Create 3x3 subplot for 9 matrix elements
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    fig.suptitle('Error in ME2EME2000 Rotation Matrix (Godot VS computed)')
    
    for i in range(3):
        for j in range(3):
            axes[i, j].plot(np.array(time_points) - J2000_EPOCH, matrix_diff[:, i, j])
            axes[i, j].set_title(f'Element [{i},{j}]')
            axes[i, j].set_xlabel('Time (MJD)')
            axes[i, j].grid(True)
            
    
    plt.tight_layout()
    plt.savefig('matrix_evolution.png', dpi=300, bbox_inches='tight')
    plt.show(block=True)
    
    print(f"Plotted evolution of {len(matrix_diff)} matrices over time")
    