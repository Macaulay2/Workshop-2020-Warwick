import pandas as pd
import gtda 

df = pd.read_csv(r'/Users/k87125/Dropbox/Math conferences/Workshop-2020-Warwick/AlgebraicStatistics/ColoredGraphicalModels/PSDL5.csv')
#print(df.head(1))
#df.dtypes
#print(type(df))

from gtda.homology import VietorisRipsPersistence
VR = VietorisRipsPersistence(homology_dimensions=[0, 1, 2])  # Parameter explained in the text
diagrams = VR.fit_transform(df)
diagrams.shape