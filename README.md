# HeatSolvers
C++ codes to solve the heat equation in one and two spatial dimensions on uniform Cartesian grids


## Description
The heat equation in one spatial dimension with Dirichlet boundary conditions is

$$
\begin{align*}
    \frac{\partial u}{\partial t}\left(t, x\right) &= \alpha\nabla^2 u\left(t, x\right) \\
    u\left(t, 0\right) &= u_0 \\
    u\left(t, L\right) &= u_L
\end{align*}
$$

where $\alpha > 0$. The codes in this repository solve this equation and its trivial 2D extension by implicit time integration (backward Euler) on a uniform Cartesian grid.


## Minimal requirements
- A C++ compiler supporting the `c++17` standard
- `cmake`
- The HDF5 library
- `python3` with `numpy`, `matplotlib`, and `h5py` to generate evolution snapshots
- `ffmpeg` to generate movies
- VTK if you want to run the interactive visualization tool for the 2D case (see `2D/Utils/VTK_visualize`)


## Usage
1. Tune the parameters in `Parameters.hh`
2. Compile:
   ```
   ./build.sh
   ```
3. Run:
   ```
   ./install/bin/<executable>
   ```
4. To plot evolution snapshots and make movies, see `Utils/README.md` within each subdirectory
5. Check out `2D/Utils/VTK_visualize` for a simple interactive visualization tool based on VTK for the 2D case
