import numpy as np

#     def get_xy_angle(self,res,moltype,time,getcoords=None):
#         res=str(res)
#         tot_xyangle=0
#         if getcoords==None:
#             getcoords=self.trajectory_to_gro(time)[0]
#         for i in range(len(self.scd_tail_atoms_of[moltype])):
#             tailcoords+=([getcoords[(moltype,str(time),str(res),x)] for x in self.scd_tail_atoms_of[moltype][i]])
#         #tailcoords=[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][0]]+[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][1]]          #Saves coordinates of residue [tail1]+[tail2]
#         headcoords=np.array(getcoords[(moltype,str(time),str(res),self.central_atom_of[moltype])])
#         geocenter=self.calculate_geometriccenter(tailcoords)
#         tiltangle=np.arccos(np.dot((geocenter-headcoords),[1,0,0])/np.linalg.norm(geocenter-headcoords))
#         return (tiltangle*180/np.pi)
     
def calculate_geometriccenter(self,coordinateinput):
    geocenter=[0,0,0]
    for atomcoords in coordinateinput:
        for dimension in atomcoords:
            geocenter[atomcoords.index(dimension)]+=dimension
    geocenter=[x/len(coordinateinput) for x in geocenter]
    return np.array(geocenter)


def rotate_by(self,coordinatelist,axisvector,angle=180):
    pass

def mirror_molecules(self):
    pass
