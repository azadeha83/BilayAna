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
import numpy as np
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from code import interact
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

Tranges = {\
           'dppc':[290, 360],\
           'dppc_chol':[290, 350],\
           'dppc_dupc':[290, 330],\
           'dppc_dupc_chol':[290, 330],\
           'dppc_dupc_chol05':[290, 330],\
           'dppc_dupc_chol40':[290, 330],\
           'dupc_chol':[270, 330],\
           }
interactions = {\
                'complete':['w_w'],\
                'head-tail': ['h_h', 'h_t', 'h_w', 't_t', 't_w', 'w_w'],\
                'head-tailhalfs':['h_h', 'h_t12', 'h_t22', 'h_w',\
                                    't12_t12', 't12_t22', 't12_w',\
                                    't22_t22', 't22_w',\
                                    'w_w'],\
                'carbons':carbonlist,\
                }



class EofScdplots():
    ''' Controls plotting of EofScd plots
        Needs raw_input files with following structure and header: 
            "Time Host Host_Scd Neib Neib_Scd Interaction DeltaScd AvgScd Etot Evdw Ecoul Nchol"
            
        --> Bug1: binAvgScd seems not to be numeric!?
    '''
    def __init__(self, interaction, systemnames=Tranges.keys()):
        self.systemnames = systemnames
        self.Tlow = {}
        self.Thigh = {}
        for sys in systemnames:
            self.Tlow[sys], self.Thigh[sys] = Tranges[sys]
        self.interaction = interaction
    def parts_together(self, pair, sysname, filename):
        '''
        Plot all interactionpairs together
        '''
        tables = {}
        tablelist = []
        for interaction_pair in interactions[self.interaction]:
            tables[sysname] = self.temperature_avg(sysname, pair, interaction_pair)
            tablelist.append((interaction_pair, interaction_pair))
        final_table = self.merge_data('final_table', 'Interaction', *tablelist)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', size='weight', color='Interaction')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth()\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'etot', breaks=[-100, 0, 20], pair=pair)\
                + ggplot2.theme_light() + my_theme()
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
        final_table = self.merge_data('final_table', 'System', *tablelist)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', size='weight', color='System')\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Etot', size='weight', color='System'))\
                + ggplot2.geom_smooth()\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'etot', breaks=[-100, 0, 20], pair=pair)\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()
    def all_systems(self, pair, interaction_pair, filename):
        '''
        Plot all systems together 
        '''
        tables = {}
        tablelist = []
        for sysname in self.systemnames:
            tables[sysname] = self.temperature_avg(sysname, pair, interaction_pair)
            tablelist.append((sysname, sysname))
        final_table = self.merge_data('final_table', 'System', *tablelist)
        plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(size='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Evdw', color='Red'))\
                + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Ecoul', color='Blue'))\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Evdw', color='Red'))\
                + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Ecoul', color='Blue'))\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-100, 0, 20], pair=pair)\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()
    def standard_plot(self, systemname, pair, interaction_pair, filename):
        ''' Creates a standard EofScd plot of one system '''
        tavg = self.temperature_avg(systemname, pair, interaction_pair)
        plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', size='weight')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth()\
                + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'etot', breaks=[-100, 0, 20], pair=pair)\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename)
        print(plot)
        gr.dev_off()

    def bin_raw_data(self, datatablename, interactions_pair, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05):
        ''' Discretize data by aggregating specified column '''
        labelbreakh = breakh-breakdx
        mydf = ro.r(\
                    '{0}$bin{1} <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))\n'
                    '{0}intermed <- {0}[Interaction=="{6}"][,Interaction:=NULL]\n'
                    '{0}fin <- {0}intermed[,lapply(.SD, mean), by=bin{1}]\n'
                    '{0}fin$weight <- {0}[,.N,by=bin{1}][,N]\n'
                    '{0}fin'\
                    .format(datatablename, column, breakl, breakh, breakdx,\
                            labelbreakh, interactions_pair)\
                    )
        return mydf
    def merge_data(self, tablename, groupname, *args):
        ''' merges data.tables
            args must be tuple: (data.table-name, framename)
            Something like: groupname = Temperature, arg=(df1, DPPC_CHOL) ...
         '''
        print('Merging data.')
        for arg in args:
            print(arg)
            ro.r('{}${}="{}"'\
                 .format(arg[0], groupname, arg[1]))
        df = ro.r('{} <- rbind({})'.format(tablename, ','.join([arg[0] for arg in args])))
        return df
    def read_EofScd(self, dataname, datafile):
        df = ro.r('{} <- fread("{}")'.format(dataname, datafile))
        return df   # SysSource type AvgScd Etot Occurrence/weight Nchol
    
    def temperature_avg(self, systemname, pair, interaction_pair):
        ''' ___ Calculate average over temperature range ___
        Variables are:
            mysys: Systemname (like dppc_chol)
            Tlow['sys'], Thigh['sys']: Temperature range
            pair: lipidtype pair (DPPC_DPPC)
            interaction: part of interaction (complete/headtail/headtailhalfs/...)
            interaction_pair: (w_w/h_t/...)
        Final data.frame structure:
            "binAvgScd Time Host Host_Scd Neib Neib_Scd DeltaScd AvgScd Etot Evdw Ecoul Nchol weight"
        '''
        framenames = []
        for temp in range(self.Tlow[systemname], self.Thigh[systemname]+1, 10):
            system = '{}_{}'.format(systemname, temp)
            file_to_read = '{}/Eofscd{}{}.dat'.format(system, pair, self.interaction)
            self.read_EofScd(system, file_to_read)
            self.bin_raw_data(system, interaction_pair)
            framenames.append(system+'fin')
            print(system+'fin')
        print(framenames)
        framenames = [(i, i[:-3]) for i in framenames]
        print(framenames)
        finalframe = 'df'
        self.merge_data(finalframe, 'temperature', *framenames)
        df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
                        '{1} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binAvgScd]'.format(finalframe, systemname)) 
        return df_final

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




