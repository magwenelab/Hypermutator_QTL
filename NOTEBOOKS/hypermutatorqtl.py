## Load in needed mods
import pandas as pd, numpy as np

## Set centromere dataframe and mat locus
Centromeres = pd.DataFrame([ ## This was taken from Vikas Yadav's paper in PNAS 2018
    [970169,835384,1370568,708804,1559983,780649,525714,## Chromosome 3 was edited here to include a 1 at beginning.
    451162,801830,199434,868824,139633,579772,441845],
    [1006931,889427,1409632,752337,1587231,821756,584338,
    512653,839446,243741,933658,171048,632362,477986]]).T

## MAT is on chromosome 5 in H99 reference genome
MAT = np.array([185700,275730])