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
           'dppc_chol':[290, 300],\
#           'dppc_dupc':[290, 330],\
           'dppc_dupc_chol':[290, 330],\
           'dppc_dupc_chol05':[290, 330],\
           'dppc_dupc_chol40':[290, 330],\
#           'dupc_chol':[270, 330],\
           }
interactions = {\
                'complete':['w_w'],\
                'head-tail': ['h_h', 'h_t', 't_t',],\
                'head-tailhalfs':['h_h', 'h_t12', 'h_t22', 'h_w',\
                                    't12_t12', 't12_t22', 't12_w',\
                                    't22_t22', 't22_w',\
                                    'w_w'],\
                'carbons':carbonlist,\
                }

#def rdfplot(datafile, outputfile_name, ncholrange=None, temperature=None):
def rdfplot(system, interaction, outputfile_name, ncholrange=None, temperature=None):
    ''' creates an rdf-plot
        with ncholrange the range for datafiles with nchol specification can be given with [start, end]
        
        datafilename path: <system>_<temperature>/datafiles/rdf/rdf_<inter1>_<inter2>.xvg
    '''
    datafile=''
    if temperature is not None:
        for temp in Tranges[system]:
            pass
    else:
        if ncholrange is not None:
            frame_regex = re.compile(r'(.*)(\..*)$')
            dataframes = {}
            framespecifiers = []
            for nchol in range(ncholrange[0], ncholrange[1]+1):
                regmatch = re.match(frame_regex, datafile)
                filesuffix = regmatch.group(2)
                dtfile_name = regmatch.group(1)+str(nchol)+filesuffix
                dataframes[nchol] = ro.r('dt{0} <- fread("{1}")'
                      '\ndt{0}$V2 <- dt{0}$V2/max(dt{0}$V2)'
                      '\ndt{0}'.format(nchol, dtfile_name))
                framespecifiers.append(['dt'+str(nchol), nchol])
            dt = merge_data('dt', 'NChol', *framespecifiers)
            aes = ggplot2.aes_string(x='V1', y='V2', color='NChol')
        else:
            dt = ro.r('{0} <- fread("{1}")'
                          '\n{0}$V2 <- {0}$V2/max({0}$V2)'
                          '\n{0}'.format('dt_raw', datafile))
            aes = ggplot2.aes_string(x='V1', y='V2')
    plot = ggplot2.ggplot(dt)\
            + aes\
            + ggplot2.geom_point(alpha='0.6')\
            + ggplot2.geom_line()\
            + add_axis('x', 'distance r / nm', breaks=[0.0, 4.0, 0.5])\
            + add_axis('y', 'g(r)', breaks=[0.0, 1.1, 0.1])\
            + ggplot2.theme_light() + my_theme()
#            + ggplot2.geom_smooth(se='False', span='0.01', n=900)\
    gr.pdf(file=outputfile_name)
    print(plot)
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
    def parts_together(self, pair, sysname, filename):
        '''
        Plot all interactionpairs together
        '''
        tables = {}
        tablelist = []
        for interaction_pair in interactions[self.interaction]:
            tables[sysname] = self.temperature_avg(sysname, pair, interaction_pair)
            tablelist.append((sysname, interaction_pair))
            print(tables[sysname])
        final_table = merge_data('final_table', 'Interaction', *tablelist)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', size='weight', color='Interaction', weight='weight')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth()\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-50, 0, 5], pair=pair)\
                + ggplot2.theme_light() + my_theme()
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
            tables[sysname] = self.temperature_avg(sysname, pair, interaction_pair, nchol=nchol)
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
                + ggplot2.scale_size("rel. sample size")\
                + ggplot2.theme_light() + my_theme()
#                + ggplot2.geom_smooth(se='False')\
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()
    def energy_components(self, pair, interaction_pair, filename):
        '''
        Plot Vdw and Coulomb separately
        '''
        tables = {}
        tablelist = []
        for sysname in self.systemnames:
            tables[sysname] = self.temperature_avg(sysname, pair, interaction_pair)
            tablelist.append((sysname, sysname))
        final_table = merge_data('final_table', 'System', *tablelist)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(size='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Evdw', color='Red'))\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Ecoul', color='Blue'))\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Evdw', color='Red'))\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Ecoul', color='Blue'))\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-100, -45, 5], pair=pair)\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()
    def standard_plot(self, systemname, pair, interaction_pair, filename, nchol=None):
        ''' Creates a standard EofScd plot of one system '''
        tavg = self.temperature_avg(systemname, pair, interaction_pair, nchol=nchol)
        if tavg is None:
            return
        if nchol is None:
            plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x='AvgScd', y='Etot')\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se='False')\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                + ggplot2.theme_light() + my_theme()
        else:
            plot = ggplot2.ggplot(tavg)\
                    + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
                    + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                    + ggplot2.geom_smooth(se='False')\
                    + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                    + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                    + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()

    def bin_raw_data(self, datatablename, interactions_pair, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05, nchol=None):
        ''' Discretize data by aggregating specified column '''
        print("Binning: {} : {}".format(datatablename, interactions_pair))
        labelbreakh = breakh-breakdx
        if nchol is None:
            mydf = ro.r(\
                        '\n{0}$bin{1} <- cut({0}${1},'
                                        'breaks=seq({2}, {3}, {4}),'
                                        'labels=seq({2}, {5}, {4}))'
                        '\n{0}intermed1 <- {0}[Interaction=="{6}"]'
#                        '\n{0}intermed1 <- {0}[Interaction=="{6}"][,Interaction:=NULL]'
    #                    \n'{0}intermed2 <- {0}[,temperature:=NULL]\n'
    #                    \n'{0}intermed2'\
#                        '\n{0}fin <- {0}intermed1[,lapply(.SD, mean), by=bin{1}]'
                        '\n{0}fin <- {0}intermed1[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=bin{1}]'
                        '\n{0}fin'\
                        '\n{0}fin$weight <- {0}[,.N,by=bin{1}][,N]'
                        '\n{0}fin'\
                        .format(datatablename, column, breakl, breakh, breakdx,\
                                labelbreakh, interactions_pair)\
                        )
        else:
            #print(ro.r('max({}$NChol)'.format(datatablename)))
            mydf = ro.r(\
                        '\n{0}$bin{1} <- cut({0}${1},'
                                        'breaks=seq({2}, {3}, {4}),'
                                        'labels=seq({2}, {5}, {4}))'
                        '\n{0}intermed1 <- {0}[Interaction=="{6}"][,Interaction:=NULL]'
    #                    \n'{0}intermed2 <- {0}[,temperature:=NULL]\n'
    #                    \n'{0}intermed2'\
                        '\n{0}fin <- {0}intermed1[,lapply(.SD, mean), by=list(bin{1}, NChol)]'
                        '\n{0}fin'
                        '\n{0}fin$weight <- {0}[,.N,by=list(bin{1}, NChol)][,N]'
                        '\n{0}fin <- {0}fin[NChol<5]'
                        '\n{0}fin'\
                        .format(datatablename, column, breakl, breakh, breakdx,\
                                labelbreakh, interactions_pair))
            #print(ro.r('max({}fin$NChol)'.format(datatablename)))
        #print(mydf)
        print("Finished binning")
        return mydf

    def read_EofScd(self, dataname, datafile):
        df = ro.r('{} <- fread("{}")'.format(dataname, datafile))
        print(ro.r('max({}$NChol)'.format(dataname)))
        return df   # SysSource type AvgScd Etot Occurrence/weight NChol
    
    def temperature_avg(self, systemname, pair, interaction_pair, nchol=None):
        ''' ___ Calculate average over temperature range ___
        Variables are:
            mysys: Systemname (like dppc_chol)
            Tlow['sys'], Thigh['sys']: Temperature range
            pair: lipidtype pair (DPPC_DPPC)
            interaction: part of interaction (complete/headtail/headtailhalfs/...)
            interaction_pair: (w_w/h_t/...)
        Final data.frame structure:
            "binAvgScd Time Host Host_Scd Neib Neib_Scd DeltaScd AvgScd Etot Evdw Ecoul NChol weight"
        '''
        framenames = []
        for temp in range(self.Tlow[systemname], self.Thigh[systemname]+1, 10):
            system = '{}_{}'.format(systemname, temp)
            print(system)
            if self.interaction == 'complete':
                file_to_read = '{}/Eofscd{}.dat'.format(system, pair)
            else:
                file_to_read = '{}/Eofscd{}{}.dat'.format(system, pair, self.interaction)
            if os.path.exists(file_to_read):
                dt = self.read_EofScd(system, file_to_read)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
            #print(dt)
        #Averaging before everything's collected=============
            self.bin_raw_data(system, interaction_pair, nchol=nchol)
            framenames.append(system+'fin')
        print(framenames)
        framenames = [(i, i[:-3]) for i in framenames]
        print(framenames)
        finalframe = 'df'
        merge_data(finalframe, 'temperature', *framenames)
        if nchol is not None:
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                         '\n{1} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binAvgScd, NChol)]'
                         #'\n{1}$NChol <- factor({1}$NChol)'
                         .format(finalframe, systemname)) 
        else:
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
                         '{1} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binAvgScd]'.format(finalframe, systemname))
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
        name = ro.r('expression(paste("H"[{}],"(S"[CD],") / kJ/mol"))'.format(pair))
    
    if breaks == 'auto':
        raise NotImplementedError("not implemented")
        breaks = 'NULL'
        limits = 'NULL'
    else:
        if limits is None:
            limits = ro.FloatVector([breaks[0], breaks[1]])
        breaks = ro.FloatVector(np.around(np.arange(breaks[0], breaks[1], breaks[2]), decimals=1))
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
                'legend.key.width': ro.r.unit(1,"lines"),\
                'legend.key.height':ro.r.unit(1.5,"lines")\
                })
    return theme

def save_plot(filename, plot, path='NULL', device='pdf', width=15, height=20, unit='mm', dpi=360):
    ''' Save plot commands '''
    ro.r('ggsave({}, plot = {}, device = {}, path = {},'\
         'scale = 1, width = {}, height = {}, units = {},'\
         'dpi = {}, )'.format(filename, plot, device, path, width, height, unit, dpi))
