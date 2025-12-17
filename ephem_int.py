# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
from godot import core, cosmos
import numpy as np
import matplotlib.pyplot as plt

def test_moon_interpolation(
    t_start: core.tempo.Epoch,
    t_end: core.tempo.Epoch,
    t_interval: float = 3600.0,
    dt: float = 60.0,  # s
):
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    
    # Interpolating Points
    time_points = core.tempo.EpochRange(t_start, t_end).createGrid(t_interval)
    mjd_points = [t.mjd() for t in time_points]
    xyz_points = [uni.frames.vector3("Earth", "Moon", "EME2000", t) for t in time_points]
    x_points = [x[0] for x in xyz_points]
    y_points = [x[1] for x in xyz_points]
    z_points = [x[2] for x in xyz_points]
    # Ephemeris Time Series 
    timetag = core.tempo.EpochRange(t_start, t_end).createGrid(dt)
    mjd_timetag = [t.mjd() for t in timetag]
    xyz_godot = [uni.frames.vector3("Earth", "Moon", "EME2000", t) for t in timetag]
    
    # Error 
    XYZ_error = []
    t_interval_min = t_interval / 60.0
    filename = "interpolation_" + str(t_interval_min) +".txt" 
    # Linear Interpolation of Ephemeris Time Series
    with open(filename, "w") as f:
        for i, tt in enumerate(mjd_timetag):
            x_interpolated = np.interp(tt, mjd_points, x_points)
            y_interpolated = np.interp(tt, mjd_points, y_points)
            z_interpolated = np.interp(tt, mjd_points, z_points)
            XYZ_error.append([x_interpolated - xyz_godot[i][0], y_interpolated - xyz_godot[i][1], z_interpolated - xyz_godot[i][2]])
           
            f.write(
                f"{tt} {x_interpolated} {y_interpolated} {z_interpolated} {xyz_godot[i][0]} {xyz_godot[i][1]} {xyz_godot[i][2]} {XYZ_error[i][0]} {XYZ_error[i][1]} {XYZ_error[i][2]}\n"
            )

    with open("moon_ephem.txt", "w") as f:
        for i, tt in enumerate(time_points):
            timestring = tt.calStr("UTC")
            f.write(
                f"{timestring}  {xyz_godot[i][0]} {xyz_godot[i][1]} {xyz_godot[i][2]} {0.0} {0.0} {0.0}\n"
            )
            
    # Plot Error
    plt.scatter(mjd_timetag, [XYZ_error[i][0] for i, tt in enumerate(mjd_timetag)], label="X Error (km)", marker = ".")
    plt.scatter(mjd_timetag, [XYZ_error[i][1] for i, tt in enumerate(mjd_timetag)], label="Y Error (km)", marker = ".")
    plt.scatter(mjd_timetag, [XYZ_error[i][2] for i, tt in enumerate(mjd_timetag)], label="Z Error (km)", marker = ".")
    plt.title("Linear Interpolation Error " + str(t_interval_min) + " min interval")
    plt.legend()
    plt.xlabel("MJD")
    plt.ylabel("Position Error (km)")
    plt.show(block = True)

if __name__ == "__main__":
    start = core.tempo.parseEpoch("2025-11-12T14:43:53.0 UTC")
    t1 = core.tempo.convert(core.tempo.TimeScale.TDB, start)
    t2 = t1 + (90.0 * 24.0* 3600.0)  # 28 days later
    interval_interpolation = 30 * 60.0  # 30 minutes
    dt = 30 * 60.0  # 30 minutes
    t1 = core.tempo.convert(core.tempo.TimeScale.TDB, t1)
    t2 = core.tempo.convert(core.tempo.TimeScale.TDB, t2)
    test_moon_interpolation(t1, t2, interval_interpolation, dt)
