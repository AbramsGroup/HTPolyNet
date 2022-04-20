from matplotlib.pyplot import get
from HTPolyNet.coordinates import get_atom_attribute
import pandas as pd
import numpy as np

np.random.seed(0)
df1 = pd.DataFrame(np.random.choice(10, (5, 4)), columns=list('ABCD'))

print(df1.to_string())

print(2,4,get_atom_attribute(df1,'D',{'A':2,'B':4}))

df1['A'].iloc[1]=8
print(df1.to_string())
