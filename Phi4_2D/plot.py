import pandas as pd 
from matplotlib import pyplot as plt 

data = pd.read_csv("results.csv")

plt.plot(data['m2'], data['M'])
plt.xlabel("m2")
plt.ylabel("M")
plt.show()