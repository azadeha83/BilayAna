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

import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
base = importr('base')
datatable = importr('data.table')
gr = importr('grDevices')

#mtcars = data(datasets).fetch('mtcars')['mtcars']
def standard_plot():
    ''' Creates a standard EofScd plot '''

class EofScdplots():
    ''' Controls plotting of EofScd plots
        Needs raw_input files with following structure and header: 
            "Time Host Host_Scd Neib Neib_Scd Interaction DeltaScd AvgScd Etot Evdw Ecoul Nchol"
    '''
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
    
    def merge_data(self, groupname, *args):
        ''' merges data.tables
            args must be tuple: (data.table-name, framename)
            Something like: groupname = Temperature, arg=(df1, DPPC_CHOL) ...
         '''
        print('Merging data.')
        for arg in args:
            print(arg)
            ro.r('{}${}="{}"'\
                 .format(arg[0], groupname, arg[1]))
        df = ro.r('df <- rbind({})'.format(','.join([arg[0] for arg in args])))
        return df
    
    def read_EofScd(self, dataname, datafile):
        df = ro.r('{} <- fread("{}")'.format(dataname, datafile))
        return df   # SysSource type AvgScd Etot Occurrence/weight Nchol
    
    def temperature_avg(self, mysys, Tlow, Thigh, pair, interaction, interaction_pair):
        ''' Average over temperature range
            mysys: Systemname (like dppc_chol)
            Tlow, Thigh: Temperature range
            pair: lipidtype pair (DPPC_DPPC)
            interaction: part of interaction (complete/headtail/headtailhalfs/...)
            interaction_pair: (w_w/h_t/...)
        '''
        framenames = []
        for temp in range(Tlow, Thigh+1, 10):
            systemname = '{}_{}'.format(mysys, temp)
            file_to_read = '{}/Eofscd{}{}.dat'.format(systemname, pair, interaction)
            self.read_EofScd(systemname, file_to_read)
            self.bin_raw_data(systemname, interaction_pair)
            framenames.append(systemname+'fin')
            print(systemname+'fin')
        print(framenames)
        framenames = [(i, i[:-3]) for i in framenames]
        print(framenames)
        #self.merge_data('temperature', *framenames)
        df = self.bin_raw_data('df', interaction, column="temperature")
        return df

    def add_specific_geoms(self, plottype):
        ''' Add geoms related to a name??? '''
    
    def add_axislabels(self, scalex=(-0.5, 1, 0.1), scaley=(-100, 0, 10), **kwargs):
        ''' Axes scales and labels
            With kwargs  x='time', y='energy' ...
                    or scalex=(0,100,dt), scaley= ...
        '''
        axisdescription = ''
        axisdescription += ggplot2.scale_x_continuous(kwargs['scalex'])
        axisdescription += ggplot2.scale_y_continuous(breaks=kwargs['scaley'])
        return axisdescription

def add_theme():
    ''' pass theme here '''
    theme = ggplot2.theme_light()\
            + ggplot2.theme(**{\
                'axis.text':        ggplot2.ggplot2.element_text(size=14, colour="black"),\
                'axis.title':       ggplot2.ggplot2.element_text(size=20),\
                'axis.ticks':       ggplot2.ggplot2.element_line(size=.5),\
                'panel.background': ggplot2.ggplot2.element_rect(fill="white"),\
                'panel.border':     ggplot2.ggplot2.element_rect(colour="black", size=1),\
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




