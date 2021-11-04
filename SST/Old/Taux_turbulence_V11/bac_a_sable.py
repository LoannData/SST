# -*- coding: utf-8 -*-
"""
Created on Fri May  5 09:04:19 2017

@author: root
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-10, 10, 100)
y = np.arctan(x)
plt.plot(x, y)
plt.show()
