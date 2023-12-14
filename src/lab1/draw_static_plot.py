from matplotlib import pyplot as plt
import numpy as np

matrix = np.loadtxt('dataT.bin', dtype='f')

plt.imshow(matrix, cmap=plt.cm.get_cmap('rainbow', 5096), vmin=matrix.min(), vmax=matrix.max())
plt.colorbar()
plt.show()

array = np.loadtxt('dataMDD.bin', dtype='f')

x = np.arange(len(array))
y = array

plt.plot(x, y)
plt.show()
