# robobee3d

Robobee 3D sim (and maybe files for planar ones as well)

## Installing

Common: https://github.com/avikde/controlutils

Pybullet 3D sim:
- `pip install pybullet`: need version 2.8.2 or above (see https://github.com/bulletphysics/bullet3/issues/2152 and https://pybullet.org/Bullet/phpBB3/viewtopic.php?f=9&t=12998). It may need some VS stuff on windows and takes a while to build, but works fine.

## Running

Pybullet 3D sim
- Generate the urdf by running `python xacro.py sdab.xacro > sdab.xacro.urdf` in the `urdf/` subdirectory (or run the task in VS code)
- For the pybullet example, run `python pybullet_hello.py`

## Contributing

- Create a branch and push to it
- Create a pull request to merge
- master should always have working code for everyone to use

# Julia

- Download from official or use choco
- install Ipopt.jl -- have to do `Pkg.build("Ipopt")` after following the instructions.
- try the hs071.jl test in flapopt
