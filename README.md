# dyablo-analysis

Analysis tools for dyablo

## Requirements

pydy requires the following to be able to compile :
  * A C++14/17 compiler.
  * HDF5 libraries compiled for serial with the C backend

## Getting the repo

Clone using submodules to get all the sources :

```
git clone --recursive https://github.com/mdelorme/pydy
```

If you directly cloned the repo without using the `--recursive` option, you can get the submodules independently :

```
git submodule init
git submodule update
```

## Compilation

Using Cmake, compilation is straightforward. In th root of the project :

```
make build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

## Implementation notes :

This library is used as a temporary measure before we have a real analysis tool in dyablo. Implementation choices can be questionable at times but feel free to add an issue if something is not working "as it should".

Please note that the AMR is "not really dealt with" for the moment, so there is no "projection" method allowing us to reproject everything on a regular cube as would be done for other analysis tools. Instead the probing method here is to probe cells by positions or id to get the relevant quantities as well as their position and bounding boxes in the domain. This can ben (as you can easily imagine) highly inefficient.

## Usage

There are two ways of using pydy : Either you make a C++ code, in which case you'll need to include the `include` directory, and link with the library created at compilation time (`libdyablo-analysis.so`), or use the python module.

To use the python module, just add the build directory to your `PYTHONPATH` environment variable (eg on Linux : `export PYTHONPATH=/my/installation/folder/build/:$PYTHONPATH`) and then all the functionalities can be used by importing the `pydy` module :

## Example

This example can be found in the `example` folder. In the following, we open a 2D set and plot the data :

Basic imports : 
```python
import pydy
import numpy as np
import matplotlib.pyplot as plt
```

Creating an XDMF/HDF5 reader, and loading a snapshot :
```python
reader = pydy.XdmfReader()
d = reader.readSnapshot('example.xmf')
```

Extracting scalar values on the domain :
```python
print('Simulation time is :',      d.getTime())
print('Total energy is :',         d.getTotalEnergy())
print('Total kinetic energy is :', d.getTotalKineticEnergy())
```

Calculating position to probe in the domain :
```python
print('Plotting slice of the data :')

### Getting all cell indices
# Probing positions, domain is 256x64 :
Nx = 256
Ny = 64
x = np.linspace(0.0, 4.0, Nx)
y = np.linspace(0.0, 1.0, Ny)
xx, yy = np.meshgrid(x, y)
xx = xx.ravel()
yy = yy.ravel()
zz = np.zeros_like(xx)
pos = np.array((xx, yy, zz)).T
```

Getting corresponding cell ids, and extracting density
```python
### Extracting density
cids = d.getCellsFromPositions(pos)
rho  = np.array(d.getDensity(cids))
rho = rho.reshape((Ny, Nx))
```

And plotting:
```python
plt.imshow(rho, extent=[0.0, 4.0, 0.0, 1.0], cmap='bwr')
plt.colorbar(label=r'$\rho$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()
```
