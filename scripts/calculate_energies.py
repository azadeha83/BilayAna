#!/bin/env python3

from bilana import gromacstoolautomator as gmx

myenergy1 = gmx.Energy('resindex_all', 'energy_recalculation.mdp', parts='head-tail')
myenergy2 = gmx.Energy('resindex_all', 'energy_recalculation.mdp', parts='head-tailhalfs')
myenergy1.write_energyfile()
myenergy2.write_energyfile()
