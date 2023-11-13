import numpy as np
import matplotlib.pyplot as plt

NX, NY, NT = 15 * 3 + 1, 12 * 3 + 1, 2
Size = [NX, NY]

# Чтение файла T1.dat
try:
    with open('c:/work/T1.dat', 'rb') as f:
        U = np.fromfile(f, dtype=np.float64, count=NX * NY).reshape(Size)
except FileNotFoundError:
    print('File "T1.dat" not found')

x = np.arange(1, NX + 1)
y = np.arange(1, NY + 1)
xx, yy = np.meshgrid(y, x)

plt.figure()
plt.surf(xx, yy, U)
plt.axis([1, NX, 1, NY, 0, 15])
plt.xlabel('X')
plt.ylabel('Y')
plt.zlabel('U')
plt.show()

# Чтение и отображение остальных файлов
basename = './data_T/T'
for i in range(2, NT + 2):
    filename = f'{basename}{i}.bin'
    try:
        with open(filename, 'rb') as f:
            U = np.fromfile(f, dtype=np.float64, count=NX * NY).reshape(Size)
    except FileNotFoundError:
        print(f'File "{filename}" not found')

    plt.figure()
    plt.surf(xx, yy, U)
    plt.axis([1, NX, 1, NY, 0, 15])
    plt.title(f'n={i} N={NT + 1}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.zlabel('U')
    plt.show()
    plt.pause(0.05)
