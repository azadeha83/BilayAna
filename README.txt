Program to automate certain analysis task on lipid bilayer trajectories.

Modules: 
	-systeminfo: Class systeminfo.SysInfo() initializes analysis and gathers all important information from inputfile.
	-gromacstoolautomator: Automates certain calculations using gromacs tools
		class Energy: All about calculating energies
		class Neighbors: All about determine arrangement of lipids; Also creates indexfile containing all lipids, their parts and whole system
	-energyfilecreator: Creates energyfiles in certain format from gromacstoolautomator.Energy(...).run_energycalculation(...)  --> all_energies.dat
	-mainanalysis:
	



File output: 1neighbor_info, 2resindex_all, 3scd_distribution, 4all_energies

Install with
pip install  --user -e /path/to/mypackage
