# PowerLLEL_Channel

PowerLLEL_Channel is a **Power**ful para**LLEL** solver for incompressible turbulent **Channel** flow. It is developed and optimized for massively-parallel Direct Numerical Simulation (DNS) of turbulent channel flow at high Reynolds number, e.g., $Re_{\tau} \approx 10^4$. The Navier-Stokes (NS) equations are solved by an explicit second-order Runge-Kutta based projection method. The second-order central difference scheme is used for the spatial discretization. The Pressure Poisson Equation (PPE) is efficiently solved by the FFT-based direct method, combined with the Parallel Diagonal Dominant (PDD) algorithm for the tridiagonal system.

PowerLLEL_Channel is mainly written in Fortran 90/2003 with reasonable modularity so that it can be extended in a sustainable manner. For performance reason, a C version of the PPE solver is also provided. Using a 2D pencil-like domain decomposition and the MPI+OpenMP hybrid parallel programming model, PowerLLEL_Channel shows excellent parallel performance with up to $10^4$ CPU cores, on a number of HPC systems, e.g., Tianhe-2A and TACC Frontera.

**Reference**

Jiabin Xie, Jianchao He, Yun Bao & Xi Chen (2021) A Low-Communication-Overhead Parallel DNS Method for the 3D Incompressible Wall Turbulence, International Journal of Computational Fluid Dynamics, 35:6, 413-432, DOI: [10.1080/10618562.2021.1971202](https://doi.org/10.1080/10618562.2021.1971202) [[arXiv preprint](https://arxiv.org/abs/2104.08863v2)]

## Building

The prerequisites are as follows:
- [CMake](https://cmake.org/) (3.16.0 or newer)
- MPI ([MPICH](https://www.mpich.org), [Open MPI](https://www.open-mpi.org) or [Intel MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html))
- [HDF5](https://www.hdfgroup.org/downloads/hdf5) (1.10.4 or newer)
- [FFTW](http://www.fftw.org) (3.3.4 or newer) / [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) (2018 or newer)
- [GPTL](https://github.com/jmrosinski/GPTL) (optional, 8.1.1 or newer)

Assume that the PowerLLEL_Channel source package is unzipped to the directory `$POWERLLEL_DIR`.

First, copy the template build script `$POWERLLEL_DIR/examples/template/build_th2a.sh` to `$POWERLLEL_DIR`.

Then, modify the build script according to the prompts in it.

Finally, execute the script.
The script will firstly build the executable `PowerLLEL` in the newly created directory `$POWERLLEL_DIR/build/install/bin`,
and then launch a simple test. Without any errors, the test will end successfully.

## Running

In general, the main `PowerLLEL` executable reads the parameter file `param.in` and outputs results `*.out/*.h5` at working directory `$WORK_DIR`. The best practice is:

- After a successful build, copy the `PowerLLEL` executable to the working directory `$WORK_DIR`.
- Copy the template parameter file `$POWERLLEL_DIR/examples/template/param.in` to `$WORK_DIR`.
- Modify simulation parameters in `param.in` as required.
- Launch a parallel simulation with the following command:
```bash
mpirun -n $nprocs ./PowerLLEL
```
Note that the total number of MPI processes `nprocs` should be equal to the product of parameters `p_row` and `p_col` in `param.in`.

## Visualizing

PowerLLEL_Channel stores 3D field data in HDF5 format (*.h5), and the file name format is `inst_xxxxxxxx_y.h5` (the string of eight 'x's indicates the timestep, character 'y' indicates a flow variable such as u, v, w or p). In order to facilitate users to visualize 3D field data directly by visualization software such as ParaView and VisIt, we provide a shell script `$POWERLLEL_DIR/utils/gen_xdmf.sh` to generate a XDMF descriptive file (\*.xdmf) for the HDF5 file. The generated XDMF file can be read directly by the visualization software mentioned above. Please refer to `gen_xdmf.sh` for detailed instructions.

## Contributors

PowerLLEL_Channel was originally developed by [Jiabin Xie](https://xiejb6.github.io/) and Jianchao He, under the guidance of Prof. Yun Bao. Later, other contributors participated in the subsequent development and optimization. Their names are listed below in alphabetical order.

- Yun Bao
- Guangnan Feng
- Junxuan Feng
- Jianchao He
- [Kan Wu](https://wu-kan.cn/)
- [Jiabin Xie](https://xiejb6.github.io/)

## Contributing

We appreciate any contributions and feedback that can improve PowerLLEL_Channel. If you wish to contribute to the code, please get in touch with the maintainers or open an Issue in the repository. Pull Requests are also welcome.

## License

PowerLLEL_Channel is distributed under the MIT license.