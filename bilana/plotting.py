'''
Controls all plotting of EofScd:
    Functions: 
        --> bin_data: 
        --> merga_data:
        --> read_EofScd:
        --> temperature_avg:
        --> R plot specific functions
        
1. read_EofScd('name', 'file.dat')
2. bin_data('name', 'interactions_pair(e.g. w_w)')
'''
import os, sys
import re
import numpy as np
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
base = importr('base')
datatable = importr('data.table')
gr = importr('grDevices')

#possible_interaction = {'complete': 'w_W',\
#                        'headtail':''}


carbonlist = []
for i in range(7):
    for j in range(7):
        if j > i:
            continue
        else:
            carbonlist.append('C{0}_C{1}'.format(i,j))

carbonlist = ['C0_C0', 'C1_C1', 'C2_C2', 'C3_C3', 'C4_C4', 'C5_C5', 'C6_C6']

Tranges = {\
           'dppc':[290, 360],\
           'dppc_chol':[290, 350],\
#           'dppc_dupc':[290, 330],\
           'dppc_dupc_chol':[290, 330],\
           'dppc_dupc_chol05':[290, 330],\
           'dppc_dupc_chol40':[290, 330],\
#           'dupc_chol':[270, 330],\
           }
#===============================================================================
# interactions = {\
#                 'complete':['w_w'],\
#                 'head-tail': ['h_h', 'h_t', 't_t',],\
# #                'head-tailhalfs':['h_h', 'h_t12', 'h_t22', 'h_w',\
# #                                    't12_t12', 't12_t22', 't12_w',\
# #                                    't22_t22', 't22_w',\
# #                                    'w_w'],\
#                 'head-tailhalfs':[\
#                                     't12_t12', 't12_t22', 't12_w',\
#                                     't22_t22', 't22_w',\
#                                     'w_w'],\
#                 'carbons':carbonlist,\
#                 }
#===============================================================================

#def rdfplot(datafile, outputfile_name, ncholrange=None, temperature=None):
def rdfplot(system, temperature, interaction, lipid='DPPC', ncholrange=None):
    ''' creates an rdf-plot
        datafilename path: <system>_<temperature>/datafiles/rdf/rdf_<inter1>_<inter2>.xvg
        system:         systemname
        temperature:    can be either a list with [Tstart, Tend] or an integer
        interaction:    tuple ('HEAD', 'HEAD') or ('HEAD', 'TAIL')
        ncholrange:     range for datafiles with nchol specification can be given with [start, end, step]
        allsys:         list of systems
    '''
    print("inputparameters:", system, temperature, interaction, lipid, ncholrange)
    if isinstance(temperature, list):
        tempframes = []
        for temp in range(temperature[0], temperature[1]+1, temperature[2]):
            if ncholrange is not None:
                outputfile_name = 'rdf_{0}_tavg_{2}_{3}{4}_nchol.pdf'.format(system, temp, lipid, *interaction)
                #frame_regex = re.compile(r'(.*)(\..*)$')
                dataframes = {}
                framespecifiers = []
                for nchol in range(ncholrange[0], ncholrange[1]+1):
                    datafile_name = '{0}_{1}/datafiles/rdf/rdf_{3}-{2}-{4}-{2}{5}.xvg'\
                                    .format(system, temp, lipid, interaction[0], interaction[1], nchol)
                    #regmatch = re.match(frame_regex, datafile)
                    #filesuffix = regmatch.group(2)
                    #dtfile_name = regmatch.group(1)+str(nchol)+filesuffix
                    dataframes[(temp, nchol)] = ro.r(\
                                                     'dt{2}{0} <- fread("{1}")'
                                                     '\ndt{2}{0}<-dt{2}{0}[V1>0.1]'
                                                     '\ndt{2}{0}$V2 <- dt{2}{0}$V2/max(dt{2}{0}$V2)'
                                                     '\ndt{2}{0}$temp <- {2}'
                                                     '\ndt{2}{0}'.format(nchol, datafile_name, temp))
                    print("dataframe:", dataframes[(temp, nchol)])
                    framespecifiers.append(['dt'+str(temp)+str(nchol), nchol])
                dt = merge_data('dt'+str(temp), 'NChol', *framespecifiers)
                tempframes.append(['dt'+str(temp), temp])
                aes = ggplot2.aes_string(x='V1', y='V2', color='NChol')
            elif isinstance(system, list):
                outputfile_name = 'rdf_all_tavg_{1}_{2}{3}.pdf'.format(temp, lipid, *interaction)
                #frame_regex = re.compile(r'(.*)(\..*)$')
                dataframes = {}
                framespecifiers = []
                for syst in system:
                    datafile_name = '{0}_{1}/datafiles/rdf/rdf_{3}-{2}-{4}-{2}None.xvg'\
                                    .format(syst, temp, lipid, interaction[0], interaction[1])
                    #regmatch = re.match(frame_regex, datafile)
                    #filesuffix = regmatch.group(2)
                    #dtfile_name = regmatch.group(1)+str(nchol)+filesuffix
                    dataframes[(temp, syst)] = ro.r(\
                                                     'dt{2}{0} <- fread("{1}")'
                                                     '\ndt{2}{0}<-dt{2}{0}[V1>0.1]'
                                                     '\ndt{2}{0}$V2 <- dt{2}{0}$V2/max(dt{2}{0}$V2)'
                                                     '\ndt{2}{0}$temp <- {2}'
                                                     '\ndt{2}{0}'.format(syst, datafile_name, temp))
                    print("dataframe:", dataframes[(temp, syst)])
                    framespecifiers.append(['dt'+str(temp)+str(syst), syst])
                dt = merge_data('dt'+str(syst)+str(temp), 'System', *framespecifiers)
                print("THIS IS DT at", temp, "\n", dt)
                tempframes.append(['dt'+str(syst)+str(temp), temp])
                aes = ggplot2.aes_string(x='V1', y='V2', color='System')
            else:
                outputfile_name = 'rdf_{0}_tavg_{2}_{3}{4}.pdf'.format(system, temp, lipid, *interaction)
                datafile_name = '{0}_{1}/datafiles/rdf/rdf_{3}-{2}-{4}-{2}None.xvg'\
                                    .format(system, temp, lipid, interaction[0], interaction[1])
                dataframes = {}
                dataframes[temp] = ro.r('dt{0} <- fread("{1}")'
                      '\ndt{0}$V2 <- dt{0}$V2/max(dt{0}$V2)'
                      '\ndt{0}$temp <- {2}'.format(temp, datafile_name, temp))
                tempframes.append(['dt'+str(temp), temp])
                aes = ggplot2.aes_string(x='V1', y='V2')
        dt = merge_data('dt_tavg', 'Temperature', *tempframes)
        print(dt)
    else:
        temp = temperature
        if ncholrange is not None:
            outputfile_name = 'rdf_{0}_{1}_{2}_{3}{4}_nchol.pdf'.format(system, temp, lipid, *interaction)
            #frame_regex = re.compile(r'(.*)(\..*)$')
            dataframes = {}
            framespecifiers = []
            for nchol in range(ncholrange[0], ncholrange[1]+1):
                datafile_name = '{0}_{1}/datafiles/rdf/rdf_{3}-{2}-{4}-{2}{5}.xvg'\
                                    .format(system, temp, lipid, interaction[0], interaction[1], nchol)
                #regmatch = re.match(frame_regex, datafile)        
                #filesuffix = regmatch.group(2)
                #dtfile_name = regmatch.group(1)+str(nchol)+filesuffix
                dataframes[nchol] = ro.r('dt{0} <- fread("{1}")'
                                          '\ndt{0}<-dt{0}[V1>0.1]'
                                          '\ndt{0}$V2 <- dt{0}$V2/max(dt{0}$V2)'
                                          '\ndt{0}'.format(nchol, datafile_name))
                framespecifiers.append(['dt'+str(nchol), nchol])
            dt = merge_data('dt', 'NChol', *framespecifiers)
            aes = ggplot2.aes_string(x='V1', y='V2', color='NChol')
        else:
            datafile_name = '{0}_{1}/datafiles/rdf/rdf_{3}-{2}-{4}-{2}None.xvg'\
                                .format(system, temp, lipid, interaction[0], interaction[1])
            outputfile_name = 'rdf_{0}_{1}_{2}_{3}{4}.pdf'.format(system, temp, lipid, *interaction)
            dt = ro.r('{0} <- fread("{1}")'
                          '\n{0}<-{0}[V1>0.1]'
                          '\n{0}$V2 <- {0}$V2/max({0}$V2)'
                          '\n{0}'.format('dt_raw', datafile_name))
            aes = ggplot2.aes_string(x='V1', y='V2')
    plot = ggplot2.ggplot(dt)\
            + aes\
            + ggplot2.geom_point(alpha='0.6')\
            + add_axis('x', 'distance r / nm', breaks=[0.0, 1.5, 0.2])\
            + add_axis('y', 'g(r)', breaks=[0.0, 1.1, 0.2])\
            + ggplot2.theme_light() + my_theme()\
#            + ggplot2.geom_smooth(se='False', span='0.01', n=900)\
    if isinstance(temperature, list):
        plot += ggplot2.facet_grid('Temperature ~ .')
    gr.pdf(file=outputfile_name, width=10, height=6)
    print(plot)
    #save_plot(outputfile_name, plot, width=25, height=15)
    gr.dev_off()

class EofScdplots():
    ''' Controls plotting of EofScd plots
        Needs raw_input files with following structure and header: 
            "Time Host Host_Scd Neib Neib_Scd Interaction DeltaScd AvgScd Etot Evdw Ecoul NChol"
            
        --> Bug1: binAvgScd seems not to be numeric!?
    '''
    def __init__(self, interaction, systemnames=Tranges.keys()):
        pairs = [('DPPC_DPPC', (-90, -50, 5)), ('DPPC_CHL1', (-50, -20, 5)), ('DUPC_CHL1', (-50, -20, 5))]
        self.enthalpyranges = {}
        self.systemnames = systemnames
        self.Tlow = {}
        self.Thigh = {}
        for syst in systemnames:
            self.Tlow[syst], self.Thigh[syst] = Tranges[syst]
        self.interaction = interaction
        for pair in pairs:
            self.enthalpyranges[pair[0]] = pair[1]
    
    #===========================================================================
    # def _specialmethod_parts(self, systemname, Tlow, Thigh, interaction, pair):
    #     ''' plot script specialized for carbons '''
    #     framenames = []
    #     for temp in range(Tlow, Thigh+1, 10):
    #         system = '{}_{}'.format(systemname, temp)
    #         print(system)
    #         if interaction == 'complete':
    #             file_to_read = '{}/Eofscd{}.dat'.format(system, pair)
    #         else:
    #             file_to_read = '{}/Eofscd{}{}.dat'.format(system, pair, interaction)
    #         if os.path.exists(file_to_read):
    #             self.read_EofScd(system, file_to_read)
    #         else:
    #             print("File {} does not exist".format(file_to_read))
    #             return None
    #         #print(dt)
    #     #Averaging before everything's collected=============
    #         self.bin_raw_data(system, system, interaction_pair, nchol=nchol)
    #         framenames.append(system)
    #     framenames = [(i, i[-3:]) for i in framenames]
    #     print("Framenames:", framenames)
    #     dt = merge_data(outputtable, 'temperature', *framenames)
    #     print("MERGED DATA:", dt)
    #     if nchol is not None:
    #         df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
    #                      '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binAvgScd, NChol)]'
    #                      #'\n{1}$NChol <- factor({1}$NChol)'
    #                      .format(outputtable)) 
    #     else:
    #         df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
    #                         '{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binAvgScd]'.format(outputtable))
    #         print("FINALFRAME:", df_final)
    #     return df_final
    #===========================================================================
    #===========================================================================
    # def parts_together(self, pair, sysname, filename, nchol=None):
    #     '''
    #     Plot all interactionpairs together
    #     '''
    #     breaks_inter = {'head-tail':[-60, 10, 10],\
    #                     'head-tailhalfs':[-25, 10, 5],\
    #                     'carbons':[-2.5, 0, 0.1],\
    #                     }
    #     tables = {}
    #     tablelist = []
    #     for interaction_pair in interactions[self.interaction]:
    #         print("AT", interaction_pair)
    #         tables[sysname+interaction_pair] = self.temperature_avg(interaction_pair, sysname, pair, interaction_pair, nchol=nchol)
    #         print("TAVG=", tables[sysname+interaction_pair])
    #         tablelist.append((interaction_pair, interaction_pair))
    #         print("TABLELIST", tablelist)
    #     final_table = merge_data('final_table', 'Interaction', *tablelist)
    #     print("FINAL TABLE:", final_table)
    #     if nchol is None:
    #         plot = ggplot2.ggplot(final_table)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', color='Interaction', weight='weight')\
    #             + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
    #             + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
    #             + ggplot2.theme_light() + my_theme()
    #     else:
    #         plot = ggplot2.ggplot(final_table)\
    #                 + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
    #                 + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.6)\
    #                 + ggplot2.geom_smooth(ggplot2.aes_string(), size=1., se='False')\
    #                 + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
    #                 + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
    #                 + ggplot2.theme_light() + my_theme()\
    #                 + ggplot2.facet_grid('Interaction ~ .')
    #     gr.pdf(file=filename)
    #     print(plot)
    #     gr.dev_off()
    #===========================================================================
    def delta_scd(self, pair, systemnames, interaction_pair, fileoutputname, nchol=None):
        ''' '''
        tables = {}
        tablelist = []
        for sysname in systemnames:
            tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair, nchol=nchol)
            tablelist.append((sysname, sysname))
#'\n{0}fin$weight <- {0}fin$weight/sum({0}fin$weight)'
        final_table = merge_data('final_table', 'System', *tablelist)
        print(final_table)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='DeltaScd',  color='System', weight='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(shape='factor(NChol)'))\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'deltascd', breaks=[0, 0.3, 0.1], pair=pair)\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' pair')\
                + ggplot2.guides(color=ggplot2.guide_legend("System"), shape=ggplot2.guide_legend("NChol"))\
                + ggplot2.theme_light() + my_theme()\
            #save_plot(filename, plot, width=25, height=20)
        gr.pdf(file=fileoutputname, width=10, height=7 )
        print(plot)
        return gr.dev_off()
    def parts_together(self, pair, sysname, filename, nchol=None):
        '''
        Plot all interactionpairs together
        '''
        if self.interaction == 'complete':
            return
        breaks_inter = {'head-tail':[-60, 10, 20],\
                        'head-tailhalfs':[-20, 0, 5],\
                        'carbons':[-2.5, 0, 0.1],\
                        }
        dt = self.temperature_avg('dataframe', sysname, pair, nchol=nchol, interaction_pair='all')
        if self.interaction == 'carbons' and nchol is None:
            dt = ro.r(\
                      'dataframe[,paste0("int",1:2) := tstrsplit(Interaction,"_")]'
                      '\ndataframe'
                      )
            print(dt)
            plot = ggplot2.ggplot(dt)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', weight='weight')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[0.2, 1.0, 0.4], limits=[-0.2, 1.0], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-2.5, 0, 1], pair=pair)\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.theme_light() + my_theme()\
                + ggplot2.facet_grid('int1 ~ int2')
            #save_plot(filename, plot, width=25, height=20)
            gr.pdf(file=filename,  width=10, height=7)
            print(plot)
            return gr.dev_off()
        elif self.interaction == 'carbons' and nchol is not None:
            dt = ro.r(\
                      'dataframe[,paste0("int",1:2) := tstrsplit(Interaction,"_")]'
                      '\ndataframe'
                      )
            print(dt)
            plot = ggplot2.ggplot(dt)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
                + ggplot2.geom_point(alpha=0.3)\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[0.2, 1.0, 0.4], limits=[-0.2, 1.0], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-2.5, 0, 1], pair=pair)\
                + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.facet_grid('int1 ~ int2')\
                + ggplot2.theme_light() + my_theme()\
            #save_plot(filename, plot, width=25, height=20)
            gr.pdf(file=filename, width=10, height=7)
            print(plot)
            return gr.dev_off()
        elif self.interaction == 'head-tailhalfs' and nchol is not None:
            dt = ro.r(\
                      'dataframe <- dataframe[Interaction %in% c("t12_t12", "t22_t22")]'
                      '\ndataframe'
                      )
            print(dt)
            plot = ggplot2.ggplot(dt)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
                + ggplot2.geom_point(alpha=0.3)\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.4], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], limits=[-25, 0], pair=pair)\
                + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.facet_grid('Interaction ~ .')\
                + ggplot2.theme_light() + my_theme()\
            #save_plot(filename, plot, width=25, height=20)
            gr.pdf(file=filename, width=10, height=7)
            print(plot)
            return gr.dev_off()
        elif nchol is None:
            plot = ggplot2.ggplot(dt)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', color='Interaction', weight='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.theme_light() + my_theme()
        else:
            plot = ggplot2.ggplot(dt)\
                    + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
                    + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.6)\
                    + ggplot2.geom_smooth(ggplot2.aes_string(), size=1., se='False')\
                    + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                    + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
                    + ggplot2.theme_light() + my_theme()\
                    + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                    + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
                    + ggplot2.facet_grid('Interaction ~ .')
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()
    def all_systems(self, pair, interaction_pair, filename, nchol=None):
        '''
        Plot all systems together 
        '''
        tables = {}
        tablelist = []
        for sysname in self.systemnames:
            tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair, nchol=nchol)
            tablelist.append((sysname, sysname))
#'\n{0}fin$weight <- {0}fin$weight/sum({0}fin$weight)'
        final_table = merge_data('final_table', 'System', *tablelist)
        print(final_table)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', color='System', weight='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                + ggplot2.scale_size("weight")\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.theme_light() + my_theme()
#                + ggplot2.geom_smooth(se='False')\
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
    def energy_components(self, pair, systemname, interaction_pair, filename):
        '''
        Plot Vdw and Coulomb separately
        '''
        #=======================================================================
        # tables = {}
        # tablelist = []
        # for sysname in self.systemnames:
        #     tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair)
        #     tablelist.append((sysname, sysname))
        # final_table = merge_data('final_table', 'System', *tablelist)
        #=======================================================================
        final_table = self.temperature_avg(systemname, systemname, pair, interaction_pair)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(size='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Evdw'), color='red')\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Ecoul'), color='blue')\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Evdw'), color='red')\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Ecoul'),color='blue')\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.theme_light() + my_theme()
                #+ add_axis('y', 'enthalpy', breaks=[-50, 0, 5], pair=pair)\
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()

    def standard_plot(self, systemname, pair, interaction_pair, filename, nchol=None):
        ''' Creates a standard EofScd plot of one system '''
        tavg = self.temperature_avg(systemname, systemname, pair, interaction_pair, nchol=nchol)
        if tavg is None:
            return
        if nchol is None:
            plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x='AvgScd', y='Etot')\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                + ggplot2.theme_light() + my_theme()
        else:
            plot = ggplot2.ggplot(tavg)\
                    + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
                    + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                    + ggplot2.geom_smooth(se='False')\
                    + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                    + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                    + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
                    + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
                    + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()

    def bin_raw_data(self, outputtablename, datatablename, interaction_pair, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05, nchol=None):
        ''' Discretize data by aggregating specified column '''
        print("Binning... {} : {}".format(datatablename, interaction_pair))
        labelbreakh = breakh-breakdx
        print("interaction pair is:", interaction_pair)
        if interaction_pair == 'all' and nchol is not None:
            mydf = ro.r(\
                '\n{0}$bin{1} <- cut({0}${1},'
                                'breaks=seq({2}, {3}, {4}),'
                                'labels=seq({2}, {5}, {4}))'
                '\nbckp_data <- {0}[NChol<5]'
                '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(bin{1}, Interaction, NChol)][NChol<5]'
                '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Interaction, NChol)][,N]'
                '\n{6}'\
                .format(datatablename, column,\
                        breakl, breakh, breakdx, labelbreakh,\
                        outputtablename)\
                )
        elif interaction_pair == 'all' and nchol is None:
            mydf = ro.r(\
                '\n{0}$bin{1} <- cut({0}${1},'
                                'breaks=seq({2}, {3}, {4}),'
                                'labels=seq({2}, {5}, {4}))'
                '\nbckp_data <- {0}'
                '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(bin{1}, Interaction)]'
                '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Interaction)][,N]'
                '\n{6}'\
                .format(datatablename, column,\
                        breakl, breakh, breakdx, labelbreakh,\
                        outputtablename)\
                )
        elif interaction_pair != 'all' and nchol is not None:
            mydf = ro.r(\
                '\n{0}$bin{1} <- cut({0}${1},'
                                'breaks=seq({2}, {3}, {4}),'
                                'labels=seq({2}, {5}, {4}))'
                '\n{0} <- {0}[Interaction=="{6}"][,Interaction:=NULL][NChol<5]'
                '\nbckp_data <- {0}'
                '\n{7} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(bin{1}, NChol)]'
                '\n{7}$weight <- bckp_data[,.N,by=list(bin{1}, NChol)][,N]'
                '\n{7}'\
                .format(datatablename, column,\
                        breakl, breakh, breakdx, labelbreakh,\
                        interaction_pair, outputtablename)\
                )
            #print(ro.r('max({}$NChol)'.format(datatablename)))
#===============================================================================
#             mydf = ro.r(\
#                         '\n{0}$bin{1} <- cut({0}${1},'
#                                         'breaks=seq({2}, {3}, {4}),'
#                                         'labels=seq({2}, {5}, {4}))'
#                         '\nbckp_data <- {0}[NChol<5]'
#                         #'\n{0}intermed1 <- {0}[Interaction=="{6}"][,Interaction:=NULL]'
#     #                    \n'{0}intermed2 <- {0}[,temperature:=NULL]\n'
#     #                    \n'{0}intermed2'\
#                         '\n{7} <- {0}[,lapply(.SD, mean), by=list(bin{1}, NChol)]'
# 
#                         '\n{7}$weight <- {0}[,.N,by=list(bin{1}, NChol)][,N]'
#                         '\n{7}'\
#                         .format(datatablename, column, breakl, breakh, breakdx,\
#                                 labelbreakh, interaction_pair, outputtablename))
#===============================================================================
            print("Max CHOL neighbors:", ro.r('max({}$NChol)'.format(outputtablename)))
        elif interaction_pair != 'all' and nchol is None:
            mydf = ro.r(\
                        '\n{0}$bin{1} <- cut({0}${1},'
                                        'breaks=seq({2}, {3}, {4}),'
                                        'labels=seq({2}, {5}, {4}))'
                        '\n{0}intermed1 <- {0}[Interaction=="{6}"][,Interaction:=NULL]'
#                        '\n{0}intermed1 <- {0}[Interaction=="{6}"][,Interaction:=NULL]'
    #                    \n'{0}intermed2 <- {0}[,temperature:=NULL]\n'
    #                    \n'{0}intermed2'\
#                        '\n{0}fin <- {0}intermed1[,lapply(.SD, mean), by=bin{1}]'
                        '\n{7} <- {0}intermed1[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=bin{1}]'
                        '\n{7}$weight <- {0}intermed1[,.N,by=bin{1}][,N]'
                        '\n{7}'\
                        .format(datatablename, column, breakl, breakh, breakdx,\
                                labelbreakh, interaction_pair, outputtablename, )\
                        )
        #print("Thats my df:", mydf)
        print("Finished binning")
        return mydf

    def read_EofScd(self, dataname, datafile):
        df = ro.r('{} <- fread("{}")'.format(dataname, datafile))
        return df   # SysSource type AvgScd Etot Occurrence/weight NChol
    
    def temperature_avg(self, outputtable, systemname, pair, interaction_pair, nchol=None):
        ''' ___ Calculate average over temperature range ___
        Variables are:
            mysys: Systemname (like dppc_chol)
            Tlow['sys'], Thigh['sys']: Temperature range
            pair: lipidtype pair (DPPC_DPPC)
            interaction: part of interaction (complete/headtail/headtailhalfs/...)
            interaction_pair: (w_w/h_t/...)
            nchol: None | not None
        Final data.frame structure:
            "binAvgScd Time Host Host_Scd Neib Neib_Scd DeltaScd AvgScd Etot Evdw Ecoul NChol weight"
        '''
        print("INPUT:", outputtable, systemname, pair, interaction_pair,)
        framenames = []
        for temp in range(self.Tlow[systemname], self.Thigh[systemname]+1, 10):
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            print(system)
            if self.interaction == 'complete':
                file_to_read = '{}/Eofscd{}.dat'.format(system, pair)
            else:
                file_to_read = '{}/Eofscd{}{}.dat'.format(system, pair, self.interaction)
            if os.path.exists(file_to_read):
                self.read_EofScd(system, file_to_read)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
        #Averaging before everything's collected=============
            df = self.bin_raw_data(system, system, interaction_pair, nchol=nchol)
            #print("binned df:", df)
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dt = merge_data(outputtable, 'temperature', *framenames)
        #print("MERGED DATA:", dt)
        if nchol is not None and interaction_pair == 'all':
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                         '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binAvgScd, Interaction, NChol)]'
                         #'\n{1}$NChol <- factor({1}$NChol)'
                         .format(outputtable))
        elif nchol is not None:
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                         '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binAvgScd, NChol)]'
                         #'\n{1}$NChol <- factor({1}$NChol)'
                         .format(outputtable))
        elif interaction_pair is 'all':
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                         '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binAvgScd, Interaction)]'
                         #'\n{1}$NChol <- factor({1}$NChol)'
                         .format(outputtable))
        elif nchol is None and interaction_pair != 'all':
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
                         '{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binAvgScd]'\
                         .format(outputtable))
        #print("FINALFRAME:", df_final)
        return df_final
        #======================================================================
        #     framenames.append(system)
        # print(framenames)
        # framenames = [(i, i[:-3]) for i in framenames]
        # print(framenames)
        # finalframe = system+'final'
        # dmerge = self.merge_data(finalframe, 'temperature', *framenames)
        # print(dmerge)
        # df_final = self.bin_raw_data(finalframe, interaction_pair)
        # return df_final
        #======================================================================

def add_axis(dim, name, breaks, limits=None, pair=None):
    ''' Axes scales and labels
    dim: x/y/z axis
    breaks: (start, end, bin)
    limits: (start, end)
    '''
    if name == 'avgscd':
        name = ro.r('expression(bar("S"[CD]))')
    elif name == 'enthalpy':
        #name = ro.r('expression(paste("H"[{}],"(S"[CD],") / kJ/mol"))'.format(pair))
        name = ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))')
    elif name == 'deltascd':
        name = ro.r('expression(paste(Delta,"S"[CD]))')
    
    if breaks == 'auto':
        raise NotImplementedError("not implemented")
        breaks = 'NULL'
        limits = 'NULL'
    else:
        if limits is None:
            limits = ro.FloatVector([breaks[0], breaks[1]])
        else:
            limits = ro.FloatVector(limits)
        breaks = ro.FloatVector(np.around(np.arange(breaks[0], breaks[1]+1, breaks[2]), decimals=1))
    expand = ro.FloatVector([0,0])
    if dim == 'x':
        axis = ggplot2.scale_x_continuous(name, limits=limits, breaks=breaks, expand=expand)
    if dim == 'y':
        axis = ggplot2.scale_y_continuous(name, limits=limits, breaks=breaks, expand=expand)
    return axis


def my_theme():
    ''' pass theme here '''
    theme = ggplot2.theme(**{\
                'axis.text':        ggplot2.ggplot2.element_text(size=14, colour="black"),\
                'axis.title':       ggplot2.ggplot2.element_text(size=20),\
                'axis.ticks':       ggplot2.ggplot2.element_line(size=.5),\
                'panel.background': ggplot2.ggplot2.element_rect(fill="white"),\
                'panel.border':     ggplot2.ggplot2.element_rect(colour="black", size=0.5),\
                'panel.grid.major': ggplot2.ggplot2.element_line(colour="grey90", size=0.5),\
                'panel.grid.minor': ggplot2.ggplot2.element_line(colour="grey90", size=0.5),\
                'legend.title':     ggplot2.ggplot2.element_text(size=20),\
                'legend.text':      ggplot2.ggplot2.element_text(size =14),\
                'legend.key':       ggplot2.ggplot2.element_blank(),\
                'plot.title':       ggplot2.ggplot2.element_text(size=20, hjust = 0.5),\
                'legend.key.width': ro.r.unit(1,"lines"),\
                'legend.key.height':ro.r.unit(1.5,"lines")\
                })
    return theme

def save_plot(filename, plot, path='.', device='pdf', width=15, height=20, unit='cm', dpi=360):
    ''' Save plot commands '''
    ro.r('ggsave("{}", plot = {}, device = "{}", path = "{}",'\
         'scale = 1, width = {}, height = {}, units = "{}",'\
         'dpi = {}, )'.format(filename, plot, device, path, width, height, unit, dpi))

def merge_data(tablename, groupname, *args):
    ''' merges data.tables
        args must be tuple: (data.table-name, frame-specifier)
        Something like: groupname = Temperature, arg=(df1, DPPC_CHOL) ...
     '''
    print('Merging data . . .')
    for arg in args:
        #print(arg)
        #print("Mergedata:", ro.r(arg[0]))
        ro.r('{}${}="{}"'.format(arg[0], groupname, arg[1]))
    all_frames = ','.join([arg[0] for arg in args])
    #print("All frames:", all_frames)
    df = ro.r('{} <- rbind({})'.format(tablename, all_frames))
    return df

