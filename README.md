# robobee3d

Robobee 3D sim (and maybe files for planar ones as well)

## Installing

Common:
- I am using python 3.7
- Numpy (use wheels from https://www.lfd.uci.edu/~gohlke/pythonlibs/ on windows)

Pybullet 3D sim:
- `pip install pybullet` (it may need some VS stuff on windows and takes a while to build, but works fine)

Planar MPC:
- scipy
- `pip install osqp`

## Running

Pybullet 3D sim
- Generate the urdf by running `python xacro.py sdab.xacro > sdab.xacro.urdf` in the `urdf/` subdirectory (or run the task in VS code)
- For the pybullet example, run `python pybullet_hello.py`

## Contributing

- Create a branch and push to it
- Create a pull request to merge
- master should always have working code for everyone to use
