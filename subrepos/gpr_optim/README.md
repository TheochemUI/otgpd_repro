# gpr_dimer

Port of the GPR-dimer code from Matlab to C++

# Setup

## Conda / Micromamba

``` bash
micromamba create -p $(pwd)/.tmp -c conda-forge compilers meson cmake gfortran eigen==3.4.0 pkg-config gtest
micromamba activate $(pwd)/.tmp
```

Or perhaps if this is used very often, install the environment globally.

``` bash
micromamba env create -f environment.yml
micromamba activate gprd
```

Alternatively, `pixi shell` can be used.

# Meson Usage

If for some reason no envionment manager is used, then one can manually setup
`gtests` as below.

1. Setup `gtests` (once)

``` bash
mkdir subprojects
meson wrap install gtest
```

Note that optionally, if not using `nix` or something else, `meson wrap install eigen` will also help

2. Compile (`input` is copied over at configuration)

``` bash
meson setup builddir
meson compile -C builddir
```

Recall that a release build with `meson` is via:

``` bash
meson setup builddir --buildtype=release --optimization=3
```

3. Profit

``` bash
cd builddir
./gprd
# Sample output
Final convergence obtained after 2 relaxation phases (total number of image evaluations: 7).
Elapsed time: 18.092s

8.98514 9.948 7.88447 7.64819 9.94644 7.88399 


-0.000319603 


0.00123983 0.0022487 0.00325383 -1.93536e-07 -0.00136041 0.000147053 
```

4. Test

To run the tests, the EAM potential binary must be in the source root.

``` bash
meson setup bbdir -Dwith_tests=True --reconfigure --buildtype="debug"
ln -sf bbdir/mEAMCUH2 .
meson test -C bbdir
```


# Documentation

To build and view the documentation, we have to obtain the theme and tag files:

```bash
# To build
pixi r mkdoxydoc
# To serve directly
pixi r srv_docs
```

Note that:

- We expect `doxygen` to be at **1.9.1**
- Use Javadoc

Once built, they can be perused with any HTTP/S server:

``` bash
# Either
darkhttpd html
# OR
python -m http.server 8888 --directory html
```
