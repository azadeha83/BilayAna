import numpy as np


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def calculate_geometriccenter(self,coordinateinput):
    geocenter=[0,0,0]
    for atomcoords in coordinateinput:
        for dimension in atomcoords:
            geocenter[atomcoords.index(dimension)]+=dimension
    geocenter=[x/len(coordinateinput) for x in geocenter]
    return np.array(geocenter)

