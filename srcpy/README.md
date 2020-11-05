# pyADCIRCDGSWEM - The Python interface of ADCIRC DG-SWEM coupler

pyADCIRCDGSWEM is ADCIRC DG-SWEM coupler itself wrapped into Python. It allows
users to access the coupler's variables and subroutines directly in Python.
The ADCIRC DG-SWEM coupler is compiled into a shared library, `pyadcircdgswem`,
that can be imported into Python as a module. `pyADCIRCDGSWEM` is a package that
imports `pyadcircdgswem` and decorates it Pythonically for easily using the
coupler in Python. [f2py](https://numpy.org/doc/stable/f2py/),
which is a part of [NumPy](https://numpy.org/), is the tool used to
implement the ADCIRC DG-SWEM coupler's Python interface.

**Note: It is quite possible that this Python interface might be useless to you
if your only goal is to run the ADCIRC DG-SWEM Coupler without hacking into it
through Python. If that is the case, stick using the compiled binary. If that
isn't the case, then we hope that `pyADCIRCDGSWEM` is of some use to you.**

## Compiling and setting up things

### CMake way
Currently, users should be able to use either the GNU Compiler Collection or
Intel compilers for compiling the code. The general CMake compilation workflow
is:
```bash
cd <adcircdgswem_root_dir>
mkdir build
cd build/
cmake  <CMake_build_options>  ..
make
```
The coupler code may be built with or without the Python interface. This README
is only concerned with the Python interface, `pyadcircdgswem`, and the Python
package built on top of it, `pyADCIRCDGSWEM`. These are built using the CMake
option given next.

#### Python Interface build
For building the Python shared library, `pyadcircdgswem`, which is used by the
Python interface/package, `pyADCIRCDGSWEM`, pass the `-DPYTHON_INTERFACE=ON`
flag to CMake. On successful compilation, the directory, 
`<adcircdgswem_root_dir>/pyADCIRCDGSWEM`,
should contain two shared libraries named `libadcircdgswem_shared.so` and
`pyadcircdgswem*.so` in case of builds on the Linux operating system.
`pyadcircdgswem*.so` is the Python interface that is accessible in Python by
running `import pyadcircdgswem`. `pyadcircdgswem` depends on
`libadcircdgswem_shared.so`, which is the shared library created out of the
original Fortran source code of ADCIRC DG-SWEM coupler. Because of that, the
`LD_LIBRARY_PATH` environment variable must be set so that `pyadcircdgswem*.so`
can find and use `libadcircdgswem_shared.so`.

The steps involved are given below:
```bash
cd <adcircdgswem_root_dir>
mkdir build
cd build/
cmake -DPYTHON_INTERFACE=ON ..
make
export LD_LIBRARY_PATH="<adcircdgswem_root_dir>/pyADCIRCDGSWEM:$LD_LIBRARY_PATH"
```

##### Optional: `pip install`
After compiling the Python interface, you can pip install `pyADCIRCDGSWEM` by
running `python3 -m pip install .` from `<adcircdgswem_root_dir>`. To uninstall
the package, you can run `python3 -m pip uninstall pyADCIRCDGSWEM`.

#### Parallel/MPI build
For building the `pyADCIRCDGSWEM` in parallel, additionally pass the
`-DUSE_MPI=ON` flag to CMake.

#### Debug build
For building the `adcircdgswem` binary or `pyADCIRCDGSWEM` in debug mode,
additionally pass the `-DBUILD_DEBUG=ON` flag to CMake.

## Using the adcircdgswem Python Interface
There two ways of using the Python interface: `pyadcircdgswem` and
`pyADCIRCDGSWEM`. `pyadcircdgswem` is the CMake-compiled Python module in the
`pyadcircdgswem*.so` file, and `pyADCIRCDGSWEM` is a Python package in the
[pyADCIRCDGSWEM](../pyADCIRCDGSWEM) folder that imports `pyadcircdgswem` and
Pythonically organizes the imported library (under development). 

The first way of using the Python interface is by directly importing the
compiled shared library into Python, `import pyadcircdgswem`. In this case, the
main ADCIRC DG-SWEM coupler program can be run as follows, for example:
```python3
import pyadcircdgswem as pyad
pyad.pyadmain.pyad_main() # Main ADCIRC DG-SWEM coupler program.
# Or you could do the following for better readability:
# from pyadcircdgswem import pyad_main as pmain
# pmain.pyad_main()
```

The second way is using by importing the `pyADCIRCDGSWEM` package,
i.e., `import pyADCIRCDGSWEM`:
```python3
import pyADCIRCDGSWEM
#import pyADCIRCDGSWEM.pyadcircdgswem # This is the same as the first way
                                      # of accessing the coupler.
pyADCIRCDGSWEM.main()
```

The recommended way of importing the coupler in Python is the second way, i.e.,
`import pyADCIRCDGSWEM`. This is intended to be the dominant way that the
library will be used in the future as the `f2py`-compiled `pyadcircdgswem`
shared library is decorated Pythonically inside the `pyADCIRCDGSWEM` package.
Use the first way, i.e., `import pyadcircdgswem`, only if you know what you are
doing.

Currently, though, `pyADCIRCDGSWEM` mainly only contains minor, incomplete unit
tests and a `main()` function that calls the ADCIRC DG-SWEM coupler's
main/driver program defined in
[the coupler's adcircdgswem\_main.F95 file](../src/adcircdgswem_main.F95), that
allows you to invoke the main function on the Linux command line using
`python3 pyADCIRCDGSWEM` or `python3 -m pyADCIRCDGSWEM`. If MPI build is also
enabled, then you can run `mpirun -np <numprocs> python3 -m pyADCIRCDGSWEM` from
the command line as well. These commands are equivalent to running the original
coupler binary `adcircdgswem_serial` and `adcircdgswem_parallel` from the
command line.

### Notes
* For now, do not move `libadcircdgswem_shared.so` and `pyadcircdgswem*.so` from
  the `<adcircdgswem_root_dir>/pyADCIRCDGSWEM` directories to other locations.
  It is better to use symlinks/shortcuts instead of moving the files around.
* Run the `pyADCIRCDGSWEM` unit tests located at
  [<adcircdgswem_root_dir>/pyADCIRCDGSWEM/tests/unit/](../pyADCIRCDGSWEM/tests/unit)
  using the appropriate version of Python. Make sure they are all successful.
  Run the command, `python3 -m unittest discover`, from the `pyADCIRCDGSWEM`
  folder or <adcircdgswem_root_dir> for automatically running serial build
  tests. Manually run all tests in case of parallel builds for now. Do not use
  the library if the tests fail. To run the tests, you will need to set the
  `LD_LIBRARY_PATH` environment variable, as mentioned in previous sections.

## License
`pyADCIRCDGSWEM` is licensed under the
[BSD 3-clause "new" or "revised" license](../LICENSE), same as that of the
original ADCIRC DG-SWEM coupler code.

