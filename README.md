## Approximating Q to accelerate modeling

The repo focuses on approximating Q effect to accelerate finite difference in both 2D/3D, acoustic/elastic modeling.

### Highlights:
- varying grid spacing in different time blocks
- 2D and 3D implementations
- acoustic and elastic implementations
- fortran and C/CUDA implementations

## Compile dependency & instrutions
This project is first written by Bob Clapp in fortran using SEPLib and Intel compilers. Conghui He continues its work by moving the 2D acoustic modeling on CPU to other difference implementations, including 3D acoustic in fortran and C, 3D elastic in C, 3D elastic in CUDA (GPU). Conghui is more familiar with Madagascar than SEPLib. So the following dependencies should be meet in order to build the project:

- Intel compilers
- GCC
- CUDA
- Madagascar
- scons

You need edit `config/config.sh` and source it via `. config/config.sh`. `.envrc` is a symbolic link to `config/config.sh`. It works if you have [direnv](http://direnv.net/) installed that will update the environment variables once you enter this project.

To compile the project, run `scons` in the root of the project after `. config.sh`

## Run
I create a two test cases in `test/` directory. Enter one of them and run `scons`.

## Data dependency
In the `src/acoustic-2d-fortran` and `src/acoustic-3d-fortran`, additional data are needed to run the program. The data is not public to others. You may ask [Bob](bob@sep.stanford.edu)'s permission if necessary.

## Acknowledgment
1. Many thanks Bob for handing over the 2D acoustic code and answering all questions.
2. Thanks Yi Shen for deriving 3D elastic Q approximation formulations.
3. Thanks Gustavo for discussing 3D elastic modeling.
4. Thanks SEPLib and Madagascar repo for the code. Part of the code in this project is borrowed from Madagascar and SEPLib. Then I add my work on it.

