import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
base = importr('base')
datatable = importr('data.table')
gr = importr('grDevices')

#mtcars = data(datasets).fetch('mtcars')['mtcars']

def bin_data(datatablename, interaction, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05):
    ''' Discretize data by aggregating specified column '''
    labelbreakh = breakh-breakdx
    mydf = ro.r(\
                '{0}$bin{1} <- cut(dat${1},'
                                'breaks=seq({2}, {3}, {4}),'
                                'labels=seq({2}, {5}, {4}))'
                '{0}[Interaction=="{6}"][,Interaction:=NULL]'
                'mydf <- {0}[,lapply(.SD, mean), by=bin{1}]'
                'mydf$weight <- {0}[,.N,by=bin{1}][N]'
                'mydf'\
                .format(datatablename, column, breakl, breakh, breakdx,\
                        labelbreakh, interaction)\
                )
    return mydf

def merge_data(*args):
    ''' merges data.tables. args must be tuple (tablename, grpname) '''
    for arg in args:
        ro.r('{}$grp={}'\
             .format(arg[0], arg[1]))
    df = ro.r('df <- rbind({})'.format(','.join(args)))
    return df

def temperature_avg(System, Thigh, Tlow, interaction):
    ''' Average over temperature range '''
    files_to_read = []
    for sys in System:
        for temp in range(Tlow, Thigh+1, 10):
            files_to_read.append('{}_{}'.format(sys, temp))
    datafiles = read_EofScd(interaction, *files_to_read)
    for data in datafiles:
        df_tavg = ro.r('{0}'
                       ''\
                    .format(data)) # How to parse the dataframe names??
    return df_tavg
        
def read_EofScd(interaction, *args):
    datafiles = []
    for i, file in enumerate(args):
        ro.r('dat{0} <- fread("{1}")'.format(i, file))
        datafiles.append(bin_data('dat'+str(i), interaction))
    return datafiles   # SysSource type AvgScd Etot Occurrence/weight Nchol

def add_specific_geoms():
    ''' Add geoms related to a name??? '''

def add_axislabels():
    ''' Axes scales and labels '''

def add_theme():
    ''' pass theme here '''
    theme = ggplot2.theme_light()\
            + ggplot2.theme(**{\
                'axis.text':        ggplot2.ggplot2.element_text(size=14, colour = "black"),\
                'axis.title':       ggplot2.ggplot2.element_text(size=20),\
                'axis.ticks':       ggplot2.ggplot2.element_line(size=.5),\
                'panel.background': ggplot2.ggplot2.element_rect(fill = "white"),\
                'panel.border':     ggplot2.ggplot2.element_rect(colour="black",size=1),\
                'panel.grid.major': ggplot2.ggplot2.element_line(colour="grey90",size=0.5),\
                'panel.grid.minor': ggplot2.ggplot2.element_line(colour="grey90",size=0.5),\
                'legend.title':     ggplot2.ggplot2.element_text(size=20),\
                'legend.text':      ggplot2.ggplot2.element_text(size =14),\
                'legend.key':       ggplot2.ggplot2.element_blank(),\
                'legend.key.width': ro.r.unit(1,"lines"),\
                'legend.key.height':ro.r.unit(1.5,"lines")\
                })
    return theme


#def open_device(devicename, geometry):
#    ''' Open a device with specific geometry? '''
def save_plot(filename, plot, path='NULL', device='pdf', width=15, height=20, unit='mm', dpi=360):
    ''' Save plot commands '''
    ro.r('ggsave({}, plot = {}, device = {}, path = {},'\
         'scale = 1, width = {}, height = {}, units = {},'\
         'dpi = {}, )'.format(filename, plot, device, path, width, height, unit, dpi))
    