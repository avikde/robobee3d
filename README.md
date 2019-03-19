# robobee3d

Robobee 3D sim (and maybe files for planar ones as well)

## Installing

Common: https://github.com/avikde/controlutils

Pybullet 3D sim:
- `pip install pybullet`: need version 2.4.8 or above (see https://github.com/bulletphysics/bullet3/issues/2152). It may need some VS stuff on windows and takes a while to build, but works fine.

## Running

Pybullet 3D sim
- Generate the urdf by running `python xacro.py sdab.xacro > sdab.xacro.urdf` in the `urdf/` subdirectory (or run the task in VS code)
- For the pybullet example, run `python pybullet_hello.py`

## Contributing

- Create a branch and push to it
- Create a pull request to merge
- master should always have working code for everyone to use
