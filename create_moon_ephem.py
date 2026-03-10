# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
from godot import core, cosmos
import datetime

def date_formatting(time_tag: str) -> str:
    """
    output : time tag in format %Y-%m-%dT%H:%M:%S.
    Args:
        time_tag : format input  %Y-%m-%dT%H:%M:%S.%f UTC
    """
    try:
        formatin = "%Y-%m-%dT%H:%M:%S.%f"
        formatout = "%Y-%m-%dT%H:%M:%S."
        date_part = time_tag.replace(" UTC", "")
        date = datetime.datetime.strptime(date_part, formatin)
        milliseconds = str(date.microsecond // 1000).zfill(2)
        return date.strftime(formatout) + milliseconds
    except Exception as e:
        raise ValueError(f"Error formatting time '{time_tag}': {e}")

def date_title(time_tag: str) -> str:
    """
    output : time tag in format %Y-%m-%dT%H:%M:%S.
    Args:
        time_tag : format input  %Y-%m-%dT%H:%M:%S.%f UTC
    """
    try:
        formatin = "%Y-%m-%dT%H:%M:%S.%f"
        formatout = "%Y%m%d"
        date_part = time_tag.replace(" UTC", "")
        date = datetime.datetime.strptime(date_part, formatin)
        return date.strftime(formatout)
    except Exception as e:
        raise ValueError(f"Error formatting time '{time_tag}': {e}")
    
def write_header(f, start: str, end: str):
    start_time_str = start.replace("UTC", "")
    stop_time_str = end.replace("UTC", "")
    f.write(f"CCSDS_OEM_VERS = 2.0\n")
    f.write(f"META_START\n")
    f.write(f"OBJECT_NAME         = Moon\n")
    f.write(f"OBJECT_ID           = Moon\n")
    f.write(f"CENTER_NAME         = EARTH\n")
    f.write(f"REF_FRAME           = EME2000\n")
    f.write(f"TIME_SYSTEM         = UTC\n")
    f.write(f"START_TIME          = {start_time_str}\n")
    f.write(f"STOP_TIME           = {stop_time_str}\n")
    f.write(f"ORIGIN              = ARPSOFT\n")
    f.write(f"META_STOP\n\n")

def moon_ephem_file(
    start: str,
    end: str,
    dt: float, 
    filename : str = "moon_ephemeris.txt"
):
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    t1 = core.tempo.parseEpoch(start)
    t2 = core.tempo.parseEpoch(end)
    
    # Moon Position in EME2000
    time_grid = core.tempo.EpochRange(t1, t2).createGrid(dt)
    xyz_timeseries = [uni.frames.vector3("Earth", "Moon", "EME2000", t) for t in time_grid]
        
    with open(filename, "w") as f:
        write_header(f, start, end)
        
        # Write data
        for i, tt in enumerate(time_grid):
            timestring = tt.calStr("UTC")
            timestring = date_formatting(timestring)
            x = xyz_timeseries[i][0]
            y = xyz_timeseries[i][1]
            z = xyz_timeseries[i][2]
            f.write(
                f"{timestring} {x:.6f} {y:.6f} {z:.6f} {0} {0} {0}\n"
                #f"{timestring} {x:.6f} {y:.6f} {z:.6f}\n"
            )
    print(filename + " created successfully!")
    return
            
if __name__ == "__main__":
    
    start = "2010-01-01T00:00:00.000 UTC"
    end = "2050-01-01T00:00:00.000 UTC"
    
    # 10 min / 30 min / 60 min
    #dt_min = [10.0, 30.0, 60.0]
    dt_min = [30, 60] 
    for dt in dt_min:
        filename = "Moon_" + date_title(start) + "_" + str(dt) + "_min.txt"
        moon_ephem_file(
            start=start,
            end=end,
            dt=dt * 60,
            filename=filename
        )
    
