#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

std = np.loadtxt(fname = './std-9.txt')
q20 = np.loadtxt(fname = './q20-9.txt')

f = plt.figure(1)
plt.plot(std, '-', label='standard')
plt.plot(q20, 'r-', label='$Q_s=Q_p=20$')
plt.xlabel('Wavenumber (1/)')
plt.ylabel('Amplitude')
plt.title('Wavefield spectrum at t=0.9s')
plt.legend(loc='lower right')

plt.show()
f.savefig('spec-9.pdf')
