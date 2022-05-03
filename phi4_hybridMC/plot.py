import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

data = pd.read_csv("data.csv")

plt.plot(range(len(data['Enew'])), data['Enew'])
plt.show()
plt.plot(range(len(data['Eold'])), data['Eold'])
plt.show()
plt.plot(range(len(data['Enew'])), data['Enew']-data['Eold'])
plt.show()
print(np.mean(data['Enew'] - data['Eold']))