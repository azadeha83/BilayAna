import numpy as np

def calculate_geometriccenter(self,coordinateinput):
    geocenter=[0,0,0]
    for atomcoords in coordinateinput:
        for dimension in atomcoords:
            geocenter[atomcoords.index(dimension)]+=dimension
    geocenter=[x/len(coordinateinput) for x in geocenter]
    return np.array(geocenter)

