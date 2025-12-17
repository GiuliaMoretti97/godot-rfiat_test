# pyright: reportRedeclaration=false
# pyright: reportUnknownArgumentType=false
# pyright: reportAttributeAccessIssue=false
import numpy as np
from godot import cosmos
from godot import core
from Lander import check_occultation


if __name__ == "__main__":
    # INPUT ####################################################################
    # optionally avoid verbose logging messages
    core.util.suppressLogger()
    # Create Universe
    uni_cfg = cosmos.util.load_yaml("./universe.yml")
    uni = cosmos.Universe(uni_cfg)
    
    # Read file
    with open("moon_test_orbiter.txt", "r") as f:
        lines = f.readlines()
    # Parse data
    data = []
    for line in lines:
        # Split by tab and take first 9 elements (3 positions x 3 coordinates each)
        parts = line.strip().split('\t')
        if len(parts) >= 9:
            # Convert comma decimal separators to dots and parse as floats
            row = [float(part) for part in parts[0:9]]
            row.append(False)
            row.append(False)
            if parts[9] == "True":
                row[9] = True
            if parts[10] == "True":
                row[10] = True
            data.append(row)

    # Convert to numpy array
    
    # GODOT VS FILE
    errors_moon = 0
    errors_earth = 0
    # RFIAT VS FILE
    errors_moon_rfiat = 0
    errors_earth_rfiat = 0
    # GODOT VS RFIAT
    errors_moon_x = 0
    errors_earth_x = 0
    tested = len(data)
    for i, dataline in enumerate(data):
        # Positions Earth Centered
        lander = dataline[0:3]
        orbiter = dataline[3:6]
        moon = dataline[6:9]
        
        # Occultation Results
        occ_earth_ref = dataline[9]
        occ_moon_ref = dataline[10]
        
        # Lander wrt Moon
        lander_to_moon = np.array(moon) - np.array(lander)
        lander_to_orbiter = np.array(orbiter) - np.array(lander)
        # Lander wrt Earth
        lander_to_earth = - np.array(lander)      
        # Occultation with respect Moon (1)
        moon_coeff = core.astro.computeOccultationCoefficient(
            lander_to_moon,
            lander_to_orbiter,
            uni.constants.getRadius("Moon"),
            0.0,
        )
        moon_margin = core.astro.computeOccultationMargin(moon_coeff).total
        occ_moon = (moon_margin >= 0.0)
        # Occultation with respect Earth (2)
        earth_coeff = core.astro.computeOccultationCoefficient(
            lander_to_earth,
            lander_to_orbiter,
            uni.constants.getRadius("Earth"),
            0.0,
        )
        earth_margin = core.astro.computeOccultationMargin(earth_coeff).total
        occ_earth = (earth_margin >= 0.0)
        if (occ_earth_ref != occ_earth):
            errors_earth += 1
        if (occ_moon_ref != occ_moon):
            errors_moon += 1
        
        # Try to use RFIAT functions 
        occ_rfiat_moon = check_occultation(
            lander,
            orbiter,
            moon,
            uni.constants.getRadius("Moon"),
        )
        occ_rfiat_earth = check_occultation(
            lander,
            orbiter,
            [0, 0, 0],
            uni.constants.getRadius("Earth"),
        )
        if (occ_moon_ref != occ_rfiat_moon):
            errors_moon_rfiat += 1
        if (occ_earth_ref != occ_rfiat_earth):
            errors_earth_rfiat += 1
        if (occ_moon != occ_rfiat_moon):
            errors_moon_x += 1
        if (occ_earth != occ_rfiat_earth):
            errors_earth_x += 1
        if ((occ_moon != occ_moon_ref) and (abs(moon_margin) > 1.0e-2)):
            print(
                i,
                # occ_earth_ref,
                occ_moon_ref,
                # occ_earth,
                occ_moon,
                moon_margin,
                # earth_coeff
            )
    print(f"Errors GODOT VS FILE Moon: {errors_moon/tested * 100}")
    print(f"Errors GODOT VS FILE Earth: {errors_earth/tested * 100}")
    print(f"Errors RFIAT VS FILE Moon RFIAT: {errors_moon_rfiat/tested * 100}")
    print(f"Errors RFIAT VS FILE Earth RFIAT: {errors_earth_rfiat/tested * 100}")
    print(f"Errors GODOT VS RFIAT Moon: {errors_moon_x/tested * 100}")
    print(f"Errors GODOT VS RFIAT Earth: {errors_earth_x/tested * 100}")
    