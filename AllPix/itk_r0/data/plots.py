import matplotlib.pyplot as plt
import numpy as np

data_front = np.genfromtxt('r0_uirrad_data_front.csv', delimiter=',', skip_header=1, skip_footer=1)
deposited, collected = data_front[:, 0], data_front[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.hist(deposited, bins=35, histtype='bar', facecolor='none', edgecolor='black')
ax2.hist(collected, bins=35, histtype='bar', facecolor='none', edgecolor='red')

ax1.set_title('Deposited')
ax2.set_title('Collected')

fig.suptitle('Front R0')

plt.show()

data_back = np.genfromtxt('r0_uirrad_data_back.csv', delimiter=',', skip_header=1, skip_footer=1)  
deposited, collected = data_back[:, 0], data_back[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.hist(deposited, bins=35, histtype='bar', facecolor='none', edgecolor='black')
ax2.hist(collected, bins=35, histtype='bar', facecolor='none', edgecolor='red')

ax1.set_title('Deposited')
ax2.set_title('Collected')

fig.suptitle('Back R0')

plt.show()
