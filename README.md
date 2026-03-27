# Heat1D
A code solving the 1D heat equation with Dirichlet boundary conditions by implicit time integration (backward Euler)

## Minimal requirements
- A C++ compiler supporting the `C++17` standard
- `cmake`
- The HDF5 library
- `python3` with `numpy`, `matplotlib`, and `h5py` to generate evolution snapshots
- `ffmpeg` to generate movies

## Usage
1. Tune the parameters in `Parameters.hh`
2. Compile:
   ```
   ./build.sh
   ```
3. Run:
   ```
   ./install/bin/Heat1D_exe
   ```
4. To plot evolution snapshots and make movies, see `Utils/README.md`
