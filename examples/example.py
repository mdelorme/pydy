import pydy
import numpy as np
import matplotlib.pyplot as plt

reader = pydy.XdmfReader()
d = reader.readSnapshot('example.xmf')

print('Simulation time is :',      d.getTime())
print('Number of dimensions is :', d.getNDim())
print('Total energy is :',         d.getTotalEnergy())
print('Total kinetic energy is :', d.getTotalKineticEnergy())

print('Plotting slice of the data :')

### Getting all cell indices
# Probing positions, domain is 256x64 :
x = np.linspace(0.0, 4.0, 256)
y = np.linspace(0.0, 1.0, 64)
xx, yy = np.meshgrid(x, y)
xx = xx.ravel()
yy = yy.ravel()
zz = np.zeros_like(xx)
pos = np.array((xx, yy, zz)).T

### Extracting density
cids = d.getCellsFromPositions(pos)
rho  = np.array(d.getDensity(cids))

rho = rho.reshape((64, 256))

plt.imshow(rho, extent=[0.0, 4.0, 0.0, 1.0], cmap='bwr')
plt.colorbar(label=r'$\rho$')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()
