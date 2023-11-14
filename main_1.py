import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()
def animate(i):
   ax.clear()
   matrix = np.loadtxt(f'./data_1/T{i+1}.bin', dtype='f')
   plt.imshow(matrix, cmap=plt.cm.get_cmap('rainbow', 5096), vmin=matrix.min(), vmax=matrix.max())

ani = animation.FuncAnimation(fig, animate, frames=300, interval=100)
plt.show()
