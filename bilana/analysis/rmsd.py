'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
from .neighbors import get_neighbor_dict
import MDAnalysis as mda
import numpy as np
import scipy
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.align
##plotting
import pandas as pd
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pickle

LOGGER = log.LOGGER

class RMSD(SysInfo):

    def __init__(self,inputfilename="inputfile"):
        super().__init__(inputfilename)
    
        self.u = mda.Universe(self.gropath,self.trjpath)
        self.neiblist = get_neighbor_dict()

    def rmsd(self, sterol, start_frame, end_frame, step):
        
        ref = MDAnalysis.Universe('{}'.format(self.initgropath))         
        print(len(self.u.trajectory))   
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        resid_list = list(set(sterol_sel.resids)) 
        
        headatms = lipidmolecules.head_atoms_of(sterol)[1:18]
        tailatms = lipidmolecules.tail_atoms_of(sterol)[0][5:]
        
        refs_head = [[] for i in range(len(resid_list))]
        refs_tail = [[] for i in range(len(resid_list))]
        
        for i_res,res in enumerate(resid_list):

            ref_head_crds = ref.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms)))).positions
            ref_tail_crds = ref.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, tailatms)))).positions
            
            refs_head[i_res]= ref_head_crds
            refs_tail[i_res]= ref_tail_crds
        
        #rmsds = [[] for i in range(len(resid_list))]
        
        with open("rmsds.dat" , 'w') as fout:
            
            fout.write('Time\tResid\tHead_rmsd\tTail_rmsd\n')
            
            #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
            for ts in self.u.trajectory[start_frame:end_frame:step]:
                for i_res,res in enumerate(resid_list):
                    
                    ref_head_crds,ref_tail_crds = refs_head[i_res], refs_tail[i_res]
                    
                    head_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms)))).positions
                    tail_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, tailatms)))).positions
                    
                    # MDAnalysis.analysis.align.alignto(self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms)))).atoms,\
                    #      ref.atoms, select='resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms))))
                    # MDAnalysis.analysis.align.alignto(self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, tailatms)))).atoms,\
                    #      ref.atoms, select='resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, tailatms))))

                    #rmsds[i].append(((self.u.trajectory.time), (MDAnalysis.analysis.rms.rmsd(ref_tail_crds, tail_crds, center= True, superposition= True)), (MDAnalysis.analysis.rms.rmsd(ref_head_crds, head_crds, center= True, superposition= True))))
                    fout.write('{}\t{}\t{}\t{}\n'.format(self.u.trajectory.time, res,\
                        (MDAnalysis.analysis.rms.rmsd(ref_tail_crds, tail_crds, center= True, superposition= True)),\
                            (MDAnalysis.analysis.rms.rmsd(ref_head_crds, head_crds, center= True, superposition= True))))
    
    def rmsf(self, sterol, start_frame, end_frame, step):
        
        ref = MDAnalysis.Universe('{}'.format(self.initgropath))     
        self.u.trajectory[start_frame:end_frame:step]   
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        resid_list = list(set(sterol_sel.resids)) 
        
        headatms = lipidmolecules.head_atoms_of(sterol)[1:18]
        tailatms = lipidmolecules.tail_atoms_of(sterol)[0][5:]
        
        #rmsds = [[] for i in range(len(resid_list))]  
            
        #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
        
        rmsf_head = []
        rmsf_tail = []

        df_head = pd.DataFrame(columns = ['resid','atom_id','head_rmsf'])
        df_tail = pd.DataFrame(columns = ['resid','atom_id','tail_rmsf'])

        for i_res,res in enumerate(resid_list):
            
            head_sel = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms))))
            tail_sel = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, tailatms))))          
                           
            MDAnalysis.analysis.align.alignto(self.u.atoms, ref.atoms, select='resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, headatms))))
            
            R_head = MDAnalysis.analysis.rms.RMSF(head_sel, start=start_frame , stop=end_frame , step=step , verbose=True).run()
            R_tail = MDAnalysis.analysis.rms.RMSF(tail_sel, start=start_frame , stop=end_frame , step=step , verbose=True).run()
            
            # R_head = MDAnalysis.analysis.rms.RMSF(head_sel)
            # R_head.run()
            df1_head = pd.DataFrame(columns = ['resid','atom_id','head_rmsf'])
            df1_tail = pd.DataFrame(columns = ['resid','atom_id','tail_rmsf'])
 
            df1_head['resid'] = [res]*len(head_sel)
            df1_head['atom_id'] = head_sel.atoms.ids
            df1_head['head_rmsf'] = R_head.rmsf
            
            df1_tail['resid'] = [res]*len(tail_sel)
            df1_tail['atom_id'] = tail_sel.atoms.ids
            df1_tail['tail_rmsf'] = R_tail.rmsf
            # R_tail = MDAnalysis.analysis.rms.RMSF(tail_sel)
            # R_tail.run()

            df_head = df_head.append(df1_head)
            df_tail = df_tail.append(df1_tail)
            
            rmsf_head.append(R_head.rmsf)
            rmsf_tail.append(R_tail.rmsf)
            
            #fout.write('{}\t{}\t{}\t{}\n'.format(self.u.trajectory.time, res, )
        df_head.to_csv('rmsf_head.csv', sep='\t')
        df_tail.to_csv('rmsf_tail.csv', sep='\t')

        # np.savetxt('rmsf_head.dat', np.array(rmsf_head), delimiter=' ')
        # np.savetxt('rmsf_tail.dat', np.array(rmsf_tail), delimiter=' ')
        # np.savetxt('rmsf_head_resnums.dat', head_sel.resnums, delimiter=' ')
        # np.savetxt('rmsf_tail_resnums.dat', tail_sel.resnums, delimiter=' ')

    def plane_fitting(self, sterol, start_frame, end_frame, step):
        
        ringatms = lipidmolecules.head_atoms_of(sterol)[1:18]
        
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        resid_list = list(set(sterol_sel.resids)) 

        # The perpendicular distance (i.e shortest distance) from a given point to a Plane is the perpendicular 
        # distance from that point to the given plane. Let the co-ordinate of the given point be (x1, y1, z1)
        # and equation of the plane be given by the equation a * x + b * y + c * z + d = 0, where a, b and c are real constants.
        
        def rotate_via_numpy(x,y, radians):
            """Use numpy to build a rotation matrix and take the dot product."""
            c, s = np.cos(radians), np.sin(radians)
            j = np.matrix([[c, s], [-s, c]])
            m = np.dot(j, [x, y])
            return m.T[:,0],m.T[:,1]
        
        def shortest_distance(x1, y1, z1, a, b, c, d):  
      
            d = np.abs((a * x1 + b * y1 + c * z1 + d))  
            e = (np.sqrt(a * a + b * b + c * c)) 
            #print("Perpendicular distance is"), d/e 
            return d/e

        with open("plane_fitting.dat" , 'w') as fout:
            
            fout.write('Time\tResid\tSum_distances\t\n')
            
            #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
            for ts in self.u.trajectory[start_frame:end_frame:step]:
                for i_res,res in enumerate(resid_list):
                                        
                    ring_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, ringatms)))).positions
                    data = np.c_[ring_crds[:,0],ring_crds[:,1],ring_crds[:,2]]
                    
                    # print(data)
                    # print(np.shape(data))
                    
                    ring_crds_rotated_x, ring_crds_rotated_y = rotate_via_numpy(ring_crds[:,0],ring_crds[:,1],np.pi/3)
                    data_rotated = np.c_[ring_crds_rotated_x,ring_crds_rotated_y,ring_crds[:,2]]
                    #data = np.array(data_rotated) 
                    
                    # print(data_rotated)
                    # print(np.shape(data_rotated))
                    
                    mn = np.min(data, axis=0)
                    mx = np.max(data, axis=0)    

                    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
                    XX = X.flatten()
                    YY = Y.flatten()
                        
                    # best-fit linear plane (1st-order)
                    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
                    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
                        
                    # evaluate it on grid
                    Z = C[0]*X + C[1]*Y + C[2]

                    sum_least_squares = 0

                    for i in range(len(data[:,0])):
                        #print(shortest_distance(data[i,0],data[i,1],data[i,2], C[0], C[1], -1, C[2]))

                        sum_least_squares += shortest_distance(data[i,0],data[i,1],data[i,2], C[0], C[1], -1, C[2])
                        #print(sum_least_squares)
                    
                    #print(sum_least_squares)
                    fout.write('{}\t{}\t{}\n'.format(self.u.trajectory.time, res, sum_least_squares))
                    '''
                    # plot points and fitted surface using Matplotlib
                    #fig1 =  plt.figure(figsize=(10, 10))
                    fig1 =  plt.figure()
                    ax = fig1.gca(projection='3d')
                    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
                    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    ax.set_zlabel('Z')
                    ax.axis('equal')
                    ax.axis('tight')
                    plt.show()
                    '''

    def methyl_angle_plane(self, sterol, start, end, interval):
        
        planeatms = lipidmolecules.head_atoms_of(sterol)
        ringatms = lipidmolecules.head_atoms_of(sterol)[1:18]
        
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        atm1, atm2, atm3, atm4 = planeatms[10], planeatms[19], planeatms[13], planeatms[18] # C10, C19, C13, C18
        
        resid_list = list(set(sterol_sel.resids)) 
        # The perpendicular distance (i.e shortest distance) from a given point to a Plane is the perpendicular 
        # distance from that point to the given plane. Let the co-ordinate of the given point be (x1, y1, z1)
        # and equation of the plane be given by the equation a * x + b * y + c * z + d = 0, where a, b and c are real constants.

        with open("methyl_angle_plane.dat" , 'w') as fout:
            
            print('{: ^10}{: ^10}{: ^15}{: ^15}{: ^15}{: ^15}{: ^15}'.format('Time','Resid','n_sterol_neib','angle1','angle2','angle1_n','angle2_n'),file=fout)
            
            dt = self.u.trajectory.dt

            start_frame = int(start / dt)
            end_frame = int(end / dt)
            time_interval = int(interval / dt)
            
            #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
            for ts in self.u.trajectory[start_frame:end_frame:time_interval]:
                for i_res,res in enumerate(resid_list):
                                        
                    coords1 = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol,res,atm1)).positions
                    coords2 = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol,res,atm2)).positions
                    coords3 = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol,res,atm3)).positions
                    coords4 = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol,res,atm4)).positions
                    
                    v1 = np.subtract(coords2[0],coords1[0])
                    v2 = np.subtract(coords4[0],coords3[0])                    
                    
                    ring_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, ringatms)))).positions
                    data = np.c_[ring_crds[:,0],ring_crds[:,1],ring_crds[:,2]]
                        
                    # best-fit linear plane (1st-order)
                    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
                    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
                    
                    plane_v = [C[0],C[1],-1]

                    alpha1 = np.arccos(np.dot(v1, plane_v)/(np.linalg.norm(v1)*np.linalg.norm(plane_v))) * (180/np.pi)
                    alpha2 = np.arccos(np.dot(v2, plane_v)/(np.linalg.norm(v2)*np.linalg.norm(plane_v))) * (180/np.pi)

                    if alpha1 > 90:
                        alpha1 = 180 - alpha1
                    
                    if alpha2 > 90:
                        alpha2 = 180 - alpha2
                    
                    theta1 = 90 - alpha1
                    theta2 = 90 - alpha2
                    
                    neibs = self.neiblist[res][float(self.u.trajectory.time)]
                    
                    neib_list = []
                    
                    if len(neibs) != 0:
                        for r in neibs:
                            resn = self.resid_to_lipid[r]
                            if resn in lipidmolecules.STEROLS:
                                neib_list.append(r)
                    
                    n_sterol = len(neib_list)
                    #print(sum_least_squares)
                    print('{: <10}{: <10}{: <15}{: <15.2f}{: <15.2f}{: <15.2f}{: <15.2f}'.format(self.u.trajectory.time, res, n_sterol, theta1, theta2, alpha1, alpha2),file=fout)

    def plane_fitting_neib(self, sterol, start_frame, end_frame, step):
        
        ringatms = lipidmolecules.head_atoms_of(sterol)[1:18]
        
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        resid_list = list(set(sterol_sel.resids)) 

        # The perpendicular distance (i.e shortest distance) from a given point to a Plane is the perpendicular 
        # distance from that point to the given plane. Let the co-ordinate of the given point be (x1, y1, z1)
        # and equation of the plane be given by the equation a * x + b * y + c * z + d = 0, where a, b and c are real constants.
        
        def rotate_via_numpy(x,y, radians):
            """Use numpy to build a rotation matrix and take the dot product."""
            c, s = np.cos(radians), np.sin(radians)
            j = np.matrix([[c, s], [-s, c]])
            m = np.dot(j, [x, y])
            return m.T[:,0],m.T[:,1]
        
        def shortest_distance(x1, y1, z1, a, b, c, d):  
      
            d = np.abs((a * x1 + b * y1 + c * z1 + d))  
            e = (np.sqrt(a * a + b * b + c * c)) 
            #print("Perpendicular distance is"), d/e 
            return d/e

        with open("plane_fitting_neib.dat" , 'w') as fout:
            
            fout.write('Time\tResid\tn_sterol_neib\tSum_distances\t\n')
            
            #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
            for ts in self.u.trajectory[start_frame:end_frame:step]:
                for i_res,res in enumerate(resid_list):
                                        
                    ring_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, ringatms)))).positions
                    data = np.c_[ring_crds[:,0],ring_crds[:,1],ring_crds[:,2]]
                    
                    # print(data)
                    # print(np.shape(data))
                    
                    ring_crds_rotated_x, ring_crds_rotated_y = rotate_via_numpy(ring_crds[:,0],ring_crds[:,1],np.pi/3)
                    data_rotated = np.c_[ring_crds_rotated_x,ring_crds_rotated_y,ring_crds[:,2]]
                    #data = np.array(data_rotated) 
                    
                    # print(data_rotated)
                    # print(np.shape(data_rotated))
                    
                    mn = np.min(data, axis=0)
                    mx = np.max(data, axis=0)    

                    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
                    XX = X.flatten()
                    YY = Y.flatten()
                        
                    # best-fit linear plane (1st-order)
                    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
                    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
                        
                    # evaluate it on grid
                    Z = C[0]*X + C[1]*Y + C[2]

                    sum_least_squares = 0

                    for i in range(len(data[:,0])):
                        #print(shortest_distance(data[i,0],data[i,1],data[i,2], C[0], C[1], -1, C[2]))

                        sum_least_squares += shortest_distance(data[i,0],data[i,1],data[i,2], C[0], C[1], -1, C[2])
                        #print(sum_least_squares)
                    
                    neibs = self.neiblist[res][float(self.u.trajectory.time)]
                    
                    neib_list = []
                    
                    if len(neibs) != 0:
                        for r in neibs:
                            resn = self.resid_to_lipid[r]
                            if resn in lipidmolecules.STEROLS:
                                neib_list.append(r)
                    
                    n_sterol = len(neib_list)
                    #print(sum_least_squares)
                    fout.write('{}\t{}\t{}\t{}\n'.format(self.u.trajectory.time, res, n_sterol, sum_least_squares))
                    '''
                    # plot points and fitted surface using Matplotlib
                    #fig1 =  plt.figure(figsize=(10, 10))
                    fig1 =  plt.figure()
                    ax = fig1.gca(projection='3d')
                    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
                    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    ax.set_zlabel('Z')
                    ax.axis('equal')
                    ax.axis('tight')
                    plt.show()
                    '''    
    def angle_tworings(self, sterol, start_frame, end_frame, step):
        
        ringatms1 = lipidmolecules.head_atoms_of(sterol)[1:11]
        ringatms2 = lipidmolecules.head_atoms_of(sterol)[11:18] + lipidmolecules.head_atoms_of(sterol)[8:10]
        
        sterol_sel = self.u.select_atoms('resname {}'.format(sterol))
        resid_list = list(set(sterol_sel.resids)) 

        with open("angle_tworings.dat" , 'w') as fout:
            
            fout.write('Time\tResid\tAngle_tworings\t\n')
            
            #for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
            for ts in self.u.trajectory[start_frame:end_frame:step]:
                for i_res,res in enumerate(resid_list):
                                        
                    ring1_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, ringatms1)))).positions
                    ring2_crds = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, res, ' '.join(map(str, ringatms2)))).positions
                    
                    data1 = np.c_[ring1_crds[:,0],ring1_crds[:,1],ring1_crds[:,2]]
                    data2 = np.c_[ring2_crds[:,0],ring2_crds[:,1],ring2_crds[:,2]]
                    
                    mn1 = np.min(data1, axis=0)
                    mx1 = np.max(data1, axis=0)     
                    mn2 = np.min(data2, axis=0)
                    mx2 = np.max(data2, axis=0)    
                    # best-fit linear plane (1st-order)
                    A1 = np.c_[data1[:,0], data1[:,1], np.ones(data1.shape[0])]
                    A2 = np.c_[data2[:,0], data2[:,1], np.ones(data2.shape[0])]
                    
                    C1,_,_,_ = scipy.linalg.lstsq(A1, data1[:,2])    # coefficients
                    C2,_,_,_ = scipy.linalg.lstsq(A2, data2[:,2])    # coefficients
                    
                    X1,Y1 = np.meshgrid(np.linspace(mn1[0], mx1[0], 20), np.linspace(mn1[1], mx1[1], 20))
                    X2,Y2 = np.meshgrid(np.linspace(mn2[0], mx2[0], 20), np.linspace(mn2[1], mx2[1], 20))
                        
                    # evaluate it on grid
                    Z1 = C1[0]*X1 + C1[1]*Y1 + C1[2]
                    Z2 = C2[0]*X2 + C2[1]*Y2 + C2[2]
                    
                    v1,v2 = [C1[0],C1[1],-1], [C2[0],C2[1],-1]
                    #print(v1,v2)
                    angle = np.arccos(np.dot(v1, v2)/((np.linalg.norm(v1))*(np.linalg.norm(v2)))) * (180/np.pi)
                    
                    if angle < 90:
                         angle = 180-angle
                    print(angle)
                    fout.write('{}\t{}\t{}\n'.format(self.u.trajectory.time, res, angle))     
        
        def analysis(self, file_name):
        
            analysis_arglist = [GMXNAME, 'analysis', '-f', file_name, '-ac', file_name+'analysis']

            out, err = exec_gromacs(analysis_arglist)

def block_average(datastream, isplot=True, maxBlockSize=0):
    """This program computes the block average of a potentially correlated timeseries "x", and 
    provides error bounds for the estimated mean <x>. 
    As input provide a vector or timeseries "x", and the largest block size.
    
    Check out writeup in the following blog posts for more:
    http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty_14.html
    http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html
    """

    Nobs         = len(datastream)           # total number of observations in datastream
    minBlockSize = 1;                        # min: 1 observation/block

    if maxBlockSize == 0:
        maxBlockSize = int(Nobs/4);        # max: 4 blocs (otherwise can't calc variance)

    NumBlocks = maxBlockSize - minBlockSize   # total number of block sizes
    blockMean = np.zeros(NumBlocks)               # mean (expect to be "nearly" constant)
    blockVar  = np.zeros(NumBlocks)               # variance associated with each blockSize
    blockCtr  = 0
    
                #
                #  blockSize is # observations/block
                #  run them through all the possibilities
                #

    for blockSize in range(minBlockSize, maxBlockSize):
        Nblock    = int(Nobs/blockSize)               # total number of such blocks in datastream
        obsProp   = np.zeros(Nblock)                  # container for parcelling block 
        # Loop to chop datastream into blocks
        # and take average
        for i in range(1,Nblock+1):
            
            ibeg = (i-1) * blockSize
            iend =  ibeg + blockSize
            obsProp[i-1] = np.mean(datastream[ibeg:iend])
        blockMean[blockCtr] = np.mean(obsProp)
        blockVar[blockCtr]  = np.var(obsProp)/(Nblock - 1)
        blockCtr += 1

    v = np.arange(minBlockSize,maxBlockSize)

    if isplot:
        plt.subplot(2,1,1)
        plt.plot(v, np.sqrt(blockVar),'ro-',lw=2)
        plt.xlabel('block size')
        plt.ylabel('std')
        plt.subplot(2,1,2)
        plt.errorbar(v, blockMean, np.sqrt(blockVar))
        plt.ylabel('<x>')
        plt.xlabel('block size')
        print('<x> = {0:f} +/- {1:f}\n'.format(blockMean[-1], np.sqrt(blockVar[-1])))
        plt.tight_layout()
        plt.show()
        
    return v, blockVar, blockMean

