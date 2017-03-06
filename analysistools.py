import math
import numpy as np
import os


def get_inputlist():
    cwd=os.getcwd()
    ls=os.listdir(cwd)
    if 'energies' in cwd:
        dppc_selflist=[]
        dppc_dppc_headtaillist=[]
        chol_selflist=[]
        dppc_chollist=[]
        dupc_chollist=[]
        for datfile in ls:
            if 'headtailhalfs' in datfile and 'DPPC_DPPC' in datfile:
                dppc_dppc_headtaillist.append(datfile)
            if 'Eofscd_dppc_dppc' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dppc_selflist.append(datfile)
            if 'Eofscd_dppc_chol' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dppc_chollist.append(datfile)
            if 'Eofscd_dupc_chol' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dupc_chollist.append(datfile)
            if 'Eofscd_chol_chol' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                chol_selflist.append(datfile)
        return dppc_dppc_headtaillist,dppc_selflist,chol_selflist,dppc_chollist,dppc_chollist,dupc_chollist
        
    if 'neighbors' in cwd:
        dppc_list=[]
        chol_list=[]
        dupc_list=[]
        for datfile in ls:
            if 'Nofscd_dppc' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dppc_list.append(datfile)
            if 'Nofscd_chol' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                chol_list.append(datfile)
            if 'Nofscd_dupc' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dupc_list.append(datfile)
        return dppc_list,chol_list,dupc_list
        
    if 'scdplots' in cwd:
        dppc_list=[]
        chol_list=[]
        dupc_list=[]
        for datfile in ls:
            if 'scd_distribution_dppc' in datfile or 'scd_distributionLength_Based_dppc' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dppc_list.append(datfile)
            if 'scd_distribution_chol' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                chol_list.append(datfile)
            if 'scd_distribution_dupc' in datfile and '.dat' in datfile and not 'meanof' in datfile:
                dupc_list.append(datfile)
        return dppc_list,chol_list,dupc_list

def calc_avg_energy_deltascd(inputfilelist,totaloutputname):
    print("Working on",totaloutputname)
    total_data=[]
    for inputfile in inputfilelist:
        print("...on file",inputfile)
        data_of_file=[]
        with open(inputfile,"r") as inpf:
            inpf.readline()
            for line in inpf:
                cols=line.split()
                if 'Nofscd' in inputfile:
                    Countneib=cols[3]
                    Scd=cols[2]
                    data_of_file.append([Scd,Countneib])
                    total_data.append([Scd,Countneib])
                if 'Eofscd' in inputfile:
                    Etot=cols[7]
                    #VDW=cols[8]
                    #Coul=cols[9]
                    deltascd=cols[5]
                    Scd=cols[6]
                    host=cols[1]
                    neib=cols[3]
                    if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                        continue
                    #data_of_file.append([Scd,Etot])
                    total_data.append([Scd,Etot,deltascd])
    totaloutput=totaloutputname+"_totalmeandeltascd.dat"
    #with open(totaloutput+"Etot.dat","w") as Etotout, open(totaloutput+"VDW.dat","w") as VDWout, open(totaloutput+"Coul.dat","w") as Coulout:
    with open(totaloutput,"w") as totalout:
#        totalout.write("{}\t{}\t{}\t{}\t{}\n".format('Average Delta Scd','Average Scd','Average energy','Standard deviation','Sample size')
        for deltascd in range(-50,1001,100):
            floatdeltascd=deltascd/1000
            meandeltascd=floatdeltascd+0.05
            print("Working on, ",meandeltascd)

            for scd in range(-500,1000,50):
                floatscd=scd/1000
                meanscd=floatscd+0.025
                print("..on ",meanscd)
                Etotvals=[]
                #VDWvals=[]
                #Coulvals=[]
                for item in total_data:
                    #if (floatscd <= float(item[0]) <= (floatscd+0.05)) and (floatdeltascd<= float(item[2]) <= (floatdeltascd+0.1)):
                    if (scd <= int(float(item[0])*1000) <= (scd+50)) and (deltascd<= int(float(item[2])*1000) <= (deltascd+100)):
                        Etotvals+=[float(item[1])]
                if len(Etotvals)==0:
                    continue
                stderrortot=np.std(Etotvals)/(len(Etotvals))**0.5
                meantot=np.mean(Etotvals)
                totalout.write("{}\t{}\t{}\t{}\t{}\n".format(meandeltascd,meanscd,meantot,stderrortot,len(Etotvals)))
                #stderrorvdw=np.std(VDWvals)
                #meanvdw=np.mean(VDWvals)
                #stderrorcoul=np.std(Coulvals)
                #meancoul=np.mean(Coulvals)
                #totalfile.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meantot,meanoflist,variance,stderror,math.log(len(valpool),10),len(valpool)))
                #print("{}\t{}\t{}\t".format(meantot,stderrortot,len(meantot),file=Etotout))
                #print("{}\t{}\t{}\t".format(meanvdw,stderrorvdw,len(meanvdw),file=VDWout))
                #print("{}\t{}\t{}\t".format(meanvdw,stderrorcoul,len(meancoul),file=Coulout))





def calc_avg_energy_temperature(inputfilelist,totaloutputname,temprange):
    print("Working on",totaloutputname)
    data_of_low=[]
    data_of_high=[]
    #total_data=[]
    Trange=[]
    for item in inputfilelist:
        Trange.append(int(item[-7:-4]))
        Trange=sorted(Trange)
    midtemp=int((Trange[0]/10+Trange[-1]/10)/2)*10


    for inputfile in inputfilelist:
        if temprange=="low":
            if int(inputfile[-7:-4])<=midtemp:
                print("...on file (low)",inputfile)
                with open(inputfile,"r") as inpf:
                    inpf.readline()
                    for line in inpf:
                        cols=line.split()
                        #if 'Nofscd' in inputfile:
                        #    Countneib=cols[3]
                        #    Scd=cols[2]
                        #    data_of_file.append([Scd,Countneib])
                        #    total_data.append([Scd,Countneib])
                        if 'Eofscd' in inputfile:
                            Etot=cols[7]
                            Scd=cols[6]
                            host=cols[1]
                            neib=cols[3]
                            if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                                continue

                            data_of_low.append([Scd,Etot])
        elif temprange=="high":
            if int(inputfile[-7:-4])>midtemp:
                print("...on file (high)",inputfile)
                with open(inputfile,"r") as inpf:
                    inpf.readline()
                    for line in inpf:
                        cols=line.split()
                        #if 'Nofscd' in inputfile:
                        #    Countneib=cols[3]
                        #    Scd=cols[2]
                        #    data_of_file.append([Scd,Countneib])
                        #    #total_data.append([Scd,Countneib])
                        if 'Eofscd' in inputfile:
                            Etot=cols[7]
                            Scd=cols[6]
                            host=cols[1]
                            neib=cols[3]
                            if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                                continue
                            data_of_high.append([Scd,Etot])
                            #total_data.append([Scd,Etot])
    outputfile1=temprange+'_'+totaloutputname+'_totalmean.dat'
    if temprange=="low":
        data=data_of_low
    elif temprange=="high":
        data=data_of_high
    with open(outputfile1,"w") as out:
        for val in range(-500,1000,50):
            decval=val/1000
            meanscd=decval+0.025
            valpool=[]
            for item in data:
                if decval <= float(item[0]) <= (decval+0.05):
                    valpool+=[float(item[1])]
            if len(valpool)==0:
                continue
            stderror=np.std(valpool)/(len(valpool))**0.5
            variance=np.var(valpool)
            meanoflist=np.mean(valpool)
            out.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist,variance,stderror,math.log(len(valpool),10),len(valpool)))



def calc_avg_energy_nchol(inputfilelist,totaloutputname,parts):
    if parts=='complete':
        interactions=['']
    elif parts=='head-tail':
        interactions=['head-tail','head-head','tail-tail']
        interactionskey=['h_h','h_t','t_t','h_w','t_w','w_w']
    elif parts=='head-tailhalfs':
        interactions=['head-tail12','tail12-tail12','head-tail22','tail22-tail22']
        interactionskey=['h_h','h_t12','t12_t12','h_t22','t22_t22','h_w','t12_w','t22_w','w_w']
    print("Working on",totaloutputname)
    total_data=[]
    for inputfile in inputfilelist:
        print("...on file",inputfile)
        data_of_file=[]
        with open(inputfile,"r") as inpf:
            inpf.readline()
            for line in inpf:
                cols=line.split()
                if 'Nofscd' in inputfile:
                    Countneib=cols[3]
                    Scd=cols[2]
                    nchols=int(cols[4])
                    data_of_file.append([Scd,Countneib,nchols])
                    total_data.append([Scd,Countneib,nchols])
                if 'Eofscd' in inputfile:
                    if parts=='complete':
                        Etot=cols[7]
                        nchol=int(cols[10])
                        Scd=cols[6]
                        host=cols[1]
                        neib=cols[3]
                        if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                            continue
                        total_data.append([Scd,Etot,nchol])
                    else:
                        Etot=float(cols[8])
                        nchol=int(cols[11])
                        Scd=float(cols[7])
                        part=cols[5]
                        host=cols[1]
                        neib=cols[3]
                        if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                            continue
                        total_data.append([Scd,Etot,nchol,part])
                        print("On lipid {}".format(host),end='\r')
    '''
    totaloutput4=totaloutputname+parts+"_totalmean_4.dat"
    totaloutput1=totaloutputname+parts+"_totalmean_1.dat"
    totaloutput2=totaloutputname+parts+"_totalmean_2.dat"
    totaloutput0=totaloutputname+parts+"_totalmean_0.dat"
    totaloutput3=totaloutputname+parts+"_totalmean_3.dat"
    '''
    outputfiles=[]
    for i in range(5):
        outputfiles.append(''.join([totaloutputname,parts,"_totalmean_",str(i),".dat"]))
        f=open(outputfiles[i],"w")
        f.close()
        
    for val in range(-500,1000,50):
        decval=val/1000
        meanscd=decval+0.025
        valpools={}
        for i in range(5):
            valpools.update({i:{}})
            for inter in interactionskey:
                valpools[i].update({inter:[]})
                
        #valpool={0:[],1:[],2:[],3:[],4:[]}
        for item in total_data:
            #if item[3]!=inter:
            #    continue
            #print(item)
            if (decval <= float(item[0]) < (decval+0.05)) and item[2]==0:
                valpools[0][item[3]].append(float(item[1]))
            elif (decval <= float(item[0]) < (decval+0.05)) and item[2]==1:
                valpools[1][item[3]].append(float(item[1]))
            elif (decval <= float(item[0]) < (decval+0.05)) and item[2]==2:
                valpools[2][item[3]].append(float(item[1]))
            elif (decval <= float(item[0]) < (decval+0.05)) and item[2]==3:
                valpools[3][item[3]].append(float(item[1]))
            elif (decval <= float(item[0]) < (decval+0.05)) and item[2]==4:
                valpools[4][item[3]].append(float(item[1]))
        #print(valpool)

    #with open(totaloutput0,"w") as totalfile0, open(totaloutput1,"w") as totalfile1,open(totaloutput2,"w") as totalfile2,open(totaloutput3,"w") as totalfile3, open(totaloutput4,"w") as totalfile4:

        for i in range(5):
            with open(outputfiles[i],"a") as outf:
                for inter in valpools[i]:
                    if len(valpools[i][inter])!=0:
                        stderror0=np.std(valpools[i][inter])/(len(valpools[i][inter]))**0.5
                        variance0=np.var(valpools[i][inter])
                        meanoflist0=np.mean(valpools[i][inter])
                        print(meanscd,"mean",i,meanoflist0)
                        print("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}{: <10}".format(meanscd,meanoflist0,variance0,stderror0,math.log(len(valpools[0]),10),len(valpools[0]),inter),file=outf)
                            
            '''
            if len(valpool[0])!=0:
                stderror0=np.std(valpool[0])/(len(valpool[0]))**0.5
                variance0=np.var(valpool[0])
                meanoflist0=np.mean(valpool[0])
                print("mean0",meanoflist0)
                totalfile0.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist0,variance0,stderror0,math.log(len(valpool[0]),10),len(valpool[0])))
            if len(valpool[1])!=0:
                stderror1=np.std(valpool[1])/(len(valpool[1]))**0.5
                variance1=np.var(valpool[1])
                meanoflist1=np.mean(valpool[1])
                print("mean1",meanoflist1)
                totalfile1.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist1,variance1,stderror1,math.log(len(valpool[1]),10),len(valpool[1])))
            if len(valpool[2])!=0:
                stderror2=np.std(valpool[2])/(len(valpool[2]))**0.5
                variance2=np.var(valpool[2])
                meanoflist2=np.mean(valpool[2])
                print("mean2",meanoflist2)
                totalfile2.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist2,variance2,stderror2,math.log(len(valpool[2]),10),len(valpool[2])))
            if len(valpool[3])!=0:
                stderror3=np.std(valpool[3])/(len(valpool[3]))**0.5
                variance3=np.var(valpool[3])
                meanoflist3=np.mean(valpool[3])
                print("mean3",meanoflist3)
                totalfile3.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist3,variance3,stderror3,math.log(len(valpool[3]),10),len(valpool[3])))
            if len(valpool[4])!=0:
                stderror4=np.std(valpool[4])/(len(valpool[4]))**0.5
                variance4=np.var(valpool[4])
                meanoflist4=np.mean(valpool[4])
                print("mean4",meanoflist4)
                totalfile4.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist4,variance4,stderror4,math.log(len(valpool[4]),10),len(valpool[4])))
            '''



def calc_avg_energy(inputfilelist,totaloutputname):
    print("Working on",totaloutputname)
    total_data=[]
    for inputfile in inputfilelist:
        print("...on file",inputfile)
        data_of_file=[]
        with open(inputfile,"r") as inpf:
            inpf.readline()
            for line in inpf:
                cols=line.split()
                if 'Nofscd' in inputfile:
                    Countneib=cols[3]
                    Scd=cols[2]
                    data_of_file.append([Scd,Countneib])
                    total_data.append([Scd,Countneib])
                if 'Eofscd' in inputfile:
                    Etot=cols[7]
                    Scd=cols[6]
                    host=cols[1]
                    neib=cols[3]
                    if (host=='372' and neib=='242') or (host=='242' and neib=='372') and ('dppc_dupc_chol' in inputfile):
                        continue
                    data_of_file.append([Scd,Etot])
                    total_data.append([Scd,Etot])
        outputfile='meanof'+inputfile
        with open(outputfile,"w") as out:
            for val in range(-500,1000,50):
                decval=val/1000
                meanscd=decval+0.025
                valpool=[]
                for item in data_of_file:
                    if decval <= float(item[0]) <= (decval+0.05):
                        valpool+=[float(item[1])]
                if len(valpool)==0:
                    continue
                stderror=np.std(valpool)
                variance=np.var(valpool)
                meanoflist=np.mean(valpool)
                out.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist,variance,stderror,math.log(len(valpool),10),len(valpool)))
    totaloutput=totaloutputname+"_totalmean.dat"
    with open(totaloutput,"w") as totalfile:
        print("..on file",totaloutput)
        for val in range(-500,1000,50):
            decval=val/1000
            meanscd=decval+0.025
            valpool=[]
            for item in total_data:
                if decval <= float(item[0]) <= (decval+0.05):
                    valpool+=[float(item[1])]
            if len(valpool)==0:
                continue
            stderror=np.std(valpool)
            variance=np.var(valpool)
            meanoflist=np.mean(valpool)
            totalfile.write("{: <10.3f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20}{: <20}\n".format(meanscd,meanoflist,variance,stderror,math.log(len(valpool),10),len(valpool)))





def calc_scd_distribution(inputfilelist,totaloutputname):
    print("Working on",totaloutputname)

    for inputfile in inputfilelist:
        total_data=[]
        print("...on file",inputfile)
        data_of_file=[]
        with open(inputfile,"r") as inpf:
            inpf.readline()
            for line in inpf:
                cols=line.split()
                Scd=cols[2]
                data_of_file.append([Scd])
                #total_data.append([Scd,Etot])
        outputfile='meanof'+inputfile
        with open(outputfile,"w") as out:
            for val in range(-500,1000,50):
                decval=val/1000
                meanscd=decval+0.025
                valpool=[]
                for item in data_of_file:
                    if decval <= float(item[0]) <= (decval+0.05):
                        valpool+=[float(item[0])]
                if len(valpool)==0:
                    continue
                number_ocurrence=len(valpool)
                normalized_ocurrence=number_ocurrence/len(data_of_file)
                total_data.append(normalized_ocurrence)
                norm_to_bin_occurrence=normalized_ocurrence/0.05
                out.write("{: <10.3f}{: <20}{: <40.20f}{: <40.20f}{: <40.20f}\n".format(meanscd,number_ocurrence,normalized_ocurrence,norm_to_bin_occurrence,math.log(norm_to_bin_occurrence)))
        if round(sum(total_data))!=1.:
            print("Not normalized... :",sum(total_data))
        elif round(sum(total_data))==1:
            print("YES! Normalized!",sum(total_data))





def draw_func(inpfile,maxdeg,nfit):
    nfit=int(nfit)
    with open(inpfile,"r") as inpf:
        for fit in range(nfit):
            begin=0
            coeffs=[]
            for line in inpf:
                if begin==1:
                    cols1=line.split(')')
                    degree=cols1[1].split()[0]
                    coefficient=(cols1[1].split()[1])
                    if coefficient=='NA':
                        continue
                    else:
                        coefficient=float(coefficient)
                    coeffs.append(coefficient)
                    if degree==str(maxdeg):
                        break
                    #print(line[29:40])
                if '(Intercept)' in line:
                    begin=1
                    cols1=line.split(')')
                    coefficient=(cols1[1].split()[0])
                    if coefficient=='NA':
                        continue
                    else:
                        coefficient=float(coefficient)
                    #print(line[29:40])
                    coeffs.append(coefficient)
            for i in range(len(coeffs)):
                try:
                    float(coeffs[i])
                except:
                    coeffs.remove(coeffs[i])

            function=[]
            for i in coeffs:
                function.append(str(i)+'*x**'+str(coeffs.index(i)))
            for i in coeffs: print(i)
            function='+'.join(function)
            print(function)

inputlists=get_inputlist()
#calc_avg_energy_nchol(inputlists[0], "DPPC_DPPC", 'head-tail')
calc_avg_energy_nchol(inputlists[0], "DPPC_DPPC", 'head-tailhalfs')
#draw_func(sys.argv[1],sys.argv[2],sys.argv[3])
