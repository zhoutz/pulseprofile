import subprocess

subprocess.run("make sd1", shell=True)
subprocess.run("rm sd1*.txt", shell=True)

for s in ["a", "b", "c", "d", "e", "f"]:
    if s == "a":
        nu = 1
        spot_theta = 90
        obs_theta = 90
        angular_radius = 0.01
    elif s == "b":
        nu = 1
        spot_theta = 90
        obs_theta = 90
        angular_radius = 1
    elif s == "c":
        nu = 200
        spot_theta = 90
        obs_theta = 90
        angular_radius = 0.01
    elif s == "d":
        nu = 200
        spot_theta = 90
        obs_theta = 90
        angular_radius = 1
    elif s == "e":
        nu = 400
        spot_theta = 60
        obs_theta = 30
        angular_radius = 1
    elif s == "f":
        nu = 400
        spot_theta = 20
        obs_theta = 80
        angular_radius = 1
    subprocess.run(
        f"./build/sd1 {nu} {spot_theta} {obs_theta} {angular_radius} {s}", shell=True
    )
subprocess.run("python3 p_sd1_all.py s", shell=True)
