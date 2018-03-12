'''
Contains all plotting schemes:

1st part: Definitions of constants and useful lists
2nd part: Definitions of useful functions for plotting
3rd part: Classes containing the actual plot routines

'''
import os#, sys
#import re
#import time
import numpy as np
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2 import rinterface as rint
from rpy2.robjects.lib.ggplot2 import aes_string



base = importr('base')
datatable = importr('data.table')
zoo = importr('zoo')
gr = importr('grDevices')

carbonlist = []
for i in range(7):
    for j in range(7):
        if j > i:
            continue
        else:
            carbonlist.append('C{0}_C{1}'.format(i,j))

carbonlist = ['C0_C0', 'C1_C1', 'C2_C2', 'C3_C3', 'C4_C4', 'C5_C5', 'C6_C6']

Tranges_whole = {\
           #====================================================================
           # 'dppc':[290, 360],
           # 'dupc':[270, 330],
           # 'dppc_chol20':[290, 350],
           # 'dppc_chol10':[290, 350],
           # 'dppc_chol30':[290, 350],
           # 'dppc_dupc':[290, 330],
           # 'dppc_dupc_chol25':[290, 330],
           # 'dppc_dupc_chol05':[290, 330],
           # 'dppc_dupc_chol40':[290, 330],
           # 'dupc_chol10':[290,330],
           # 'dupc_chol20':[270, 330],
           # 'dupc_chol30':[290,330],
           #====================================================================
           'dppc':[i for i in range(290, 360+1, 10)],
           'dupc':[i for i in range(270, 330+1, 10)],
           'dppc_chol10':[i for i in range(290, 350+1, 10)],
           'dppc_chol20':[i for i in range(290, 350+1, 10)],
           'dppc_chol30':[i for i in range(290, 350+1, 10)],
           'dppc_dupc':[i for i in range(290, 330+1, 10)],
           'dppc_dupc_chol25':[i for i in range(290, 330+1, 10)],
           'dppc_dupc_chol05':[i for i in range(290, 330+1, 10)],
           'dppc_dupc_chol40':[i for i in range(290, 330+1, 10)],
           'dupc_chol10':[290,],
           'dupc_chol20':[i for i in range(290, 330+1, 10)],
           'dupc_chol30':[290,],
           }

Tranges_high = {\
           'dppc':[i for i in range(330, 360+1, 10)],
           'dupc':[i for i in range(290, 330+1, 10)],
           'dppc_chol10':[i for i in range(330, 350+1, 10)],
           'dppc_chol20':[i for i in range(330, 350+1, 10)],
           'dppc_chol30':[i for i in range(330, 350+1, 10)],
           'dppc_dupc':[i for i in range(300, 330+1, 10)],
           'dppc_dupc_chol25':[i for i in range(310, 330+1, 10)],
           'dppc_dupc_chol05':[i for i in range(310, 330+1, 10)],
           'dppc_dupc_chol40':[i for i in range(310, 330+1, 10)],
           'dupc_chol10':[330,],
           'dupc_chol20':[i for i in range(300, 330+1, 10)],
           'dupc_chol30':[330],
           #====================================================================
           # 'dppc':[330, 360],
           # 'dupc':[290, 330],
           # 'dppc_chol20':[330, 350],
           # 'dppc_chol10':[330, 350],
           # 'dppc_chol30':[330, 350],
           # 'dppc_dupc':[300, 330],
           # 'dppc_dupc_chol25':[310, 330],
           # 'dppc_dupc_chol05':[310, 330],
           # 'dppc_dupc_chol40':[310, 330],
           # 'dupc_chol10':[330,],
           # 'dupc_chol20':[300, 330],
           # 'dupc_chol30':[330,],
           #====================================================================
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


def read_datafile(dataname, datafile):
    df = ro.r('{} <- fread("{}")'.format(dataname, datafile))
    return df   # SysSource type AvgScd Etot Occurrence/weight NChol


def add_axis(dim, name, breaks, limits=None, pair=None):
    ''' Axes scales and labels
    dim: x/y/z axis
    breaks: (start, end, bin)
    limits: (start, end)
    '''
    if name == 'avgscd':
        name = ro.r('expression(bar("S"[CD]))')
    elif name == 'scd':
        name = ro.r('expression("S"[CD])')
    elif name == 'enthalpy':
        #name = ro.r('expression(paste("H"[{}],"(S"[CD],") / kJ/mol"))'.format(pair))
        name = ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))')
    elif name == 'deltascd':
        name = ro.r('expression(paste(Delta,"S"[CD]))')
    elif name == 'neighbors':
        name = ro.r('expression(paste("N","(S"[CD],")"))')
    
    if breaks == 'auto':
        breaks = rint.NULL
        limits = rint.NULL
    else:
        if limits is None:
            limits = ro.FloatVector([breaks[0], breaks[1]])
        else:
            limits = ro.FloatVector(limits)
        breaks = ro.FloatVector(np.around(np.arange(breaks[0], breaks[1]+1, breaks[2]), decimals=2))
    expand = ro.FloatVector([0,0])
    if dim == 'x':
        axis = ggplot2.scale_x_continuous(name, limits=limits, breaks=breaks, expand=expand)
    if dim == 'y':
        axis = ggplot2.scale_y_continuous(name, limits=limits, breaks=breaks, expand=expand)
    return axis


def my_theme():
    ''' pass theme here '''
    theme = ggplot2.theme(**{
                'axis.text':        ggplot2.ggplot2.element_text(size=25, colour="black"),
                'axis.title':       ggplot2.ggplot2.element_text(size=30),
                'axis.ticks':       ggplot2.ggplot2.element_line(size=1),
                'panel.background': ggplot2.ggplot2.element_rect(fill="white"),
                'panel.border':     ggplot2.ggplot2.element_rect(colour="black", size=0.5),
                'panel.grid.major': ggplot2.ggplot2.element_line(colour="grey90", size=0.5),
                'panel.grid.minor': ggplot2.ggplot2.element_line(colour="grey90", size=0.5),
                'legend.title':     ggplot2.ggplot2.element_text(size=30),
                'legend.text':      ggplot2.ggplot2.element_text(size =25),
                'legend.key':       ggplot2.ggplot2.element_blank(),
                'plot.title':       ggplot2.ggplot2.element_text(size=30, hjust = 0.5),
                'legend.key.width': ro.r.unit(1,"lines"),
                'legend.key.height':ro.r.unit(1.5,"lines"),
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
    print('Merging data . . .', tablename, groupname, args)
    for arg in args:
        #print(arg)
        print("Mergedata:", ro.r(arg[0]))
        data = ro.r('{0}${1}="{2}"'.format(arg[0], groupname, arg[1]))
        #print(data)
    all_frames = ','.join([arg[0] for arg in args])
    #print("All frames:", all_frames)
    df = ro.r('{} <- rbind({})'.format(tablename, all_frames))
    return df


def poly_fit_data(datatable, xcol, ycol, subset=None, deg=6):
    ''' Fit data of inputtable '''
    if subset is not None:
        fit = ro.r('form <- {0}${2}~poly({0}${1},{3},raw=TRUE)'
                   '\nm <- lm(form, data={0}, weight={0}$weight)'
                   'summary(m)'
                   .format(datatable, xcol, ycol, deg))
    print(fit)
    return fit





def p_of_N(datafilename, lipidtype, outputfilename, count_type='Ntot'):
    ''' '''
    df = read_datafile('df', datafilename)
    print(df)
    aes = ggplot2.aes_string(x='Ntot', fill='factor(Scd)', y=lipidtype)
    plot = ggplot2.ggplot(df)\
            + aes\
            + ggplot2.ggtitle(lipidtype)\
            + ggplot2.geom_bar(stat="identity", position="dodge")\
            + ggplot2.theme_light() + my_theme()\
            + ggplot2.scale_y_continuous('Frequency')\
            + add_axis('x', 'Number of '+count_type, breaks=[0, 10, 1])
    gr.pdf(file=outputfilename, width=10, height=6)
    print(plot)
    gr.dev_off()
    ro.r('rm(list=ls(all=TRUE))')
#            + add_axis('y', 'Frequency', breaks=[0.0, 25000, 2000])

#===============================================================================
# def Nrplot(system, temperature, interaction):
#     if isinstance(temperature, list):
#         tempframes = []
#         for temp in range(temperature[0], temperature[1]+1, temperature[2]):
#             outputfile_name = 'nr_{0}_{2}{3}.pdf'.format(system, temp, *interaction)
#             datafile_name = '{0}_{1}/datafiles/rdf/nr_{2}-{3}.xvg'\
#                                 .format(system, temp, interaction[0], interaction[1])
#             dataframes = {}
#             dataframes[temp] = ro.r('dt{0} <- fread("{1}")'
#                   '\ndt{0}$temp <- {2}'.format(temp, datafile_name, temp))
#             tempframes.append(['dt'+str(temp), temp])
#         aes = ggplot2.aes_string(x='V1', y='V2')
#         dt = merge_data('dt_tavg', 'Temperature', *tempframes)
#     else:
#         datafile_name = '{0}_{1}/datafiles/rdf/nr_{2}-{3}.xvg'\
#                                 .format(system, temp, interaction[0], interaction[1])
#         outputfile_name = 'nr_{0}_{1}_{2}_{3}.pdf'.format(system, temp, *interaction)
#         dt = ro.r('{0} <- fread("{1}")'
#                       '\n{0}<-{0}[V1>0.1]'
#                       '\n{0}'.format('dt_raw', datafile_name))
#         aes = ggplot2.aes_string(x='V1', y='V2')
#     if interaction == ('P', 'P'):
#         ybreaks = [0.0, 6.5, 2.0]
#         yline = 4
#     else:
#         ybreaks = [0.0, 12.0, 3.0]
#         if system == 'dppc_chol20':
#             yline = 4.0
#         elif system == 'dppc_chol30':
#             yline = 4.0
#         elif system == 'dppc_chol10':
#             yline = 4.0
#         else:
#             yline = 4.0
#     
#     plot = ggplot2.ggplot(dt)\
#         + aes\
#         + ggplot2.geom_point(alpha='0.6')\
#         + ggplot2.geom_hline(yintercept=yline)\
#         + add_axis('x', 'distance r / nm', breaks=[0.4, 1.2, 0.2])\
#         + add_axis('y', "N(r)", breaks=ybreaks)\
#         + ggplot2.ggtitle("{}: {}-{}".format(system, *interaction))\
#         + ggplot2.theme_light() + my_theme()\
# #            + ggplot2.geom_smooth(se='False', span='0.01', n=900)\
#     if isinstance(temperature, list):
#         plot += ggplot2.facet_grid('Temperature ~ .')
#     gr.pdf(file=outputfile_name, width=10, height=6)
#     print(plot)
#     #save_plot(outputfile_name, plot, width=25, height=15)
#     gr.dev_off()
#     ro.r('rm(list=ls(all=TRUE))')
#===============================================================================

class Noft():
    ''' All around Noft data '''
    @staticmethod
    def plot(systemname, outputfilename, lipid, neighbor):
        ''' Plot N(t) '''
        framenames = []
        temperatures = {}
        temperatures[systemname] = Tranges_whole[systemname] 
        for temp in temperatures[systemname]:
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            print(system)
            file_to_read='{}/Nofscd.dat'.format(system)
            if os.path.exists(file_to_read):
                df = read_datafile(system, file_to_read)
                print(df)   
                df = ro.r('\n{0}$Time = {0}$Time/1000'
                          '\n{0} <- {0}[Lipid_type == "{1}"]'
                          '\n{0}[,Lipid_type:=NULL]'
                          '\n{0} <- {0}[,lapply(.SD, mean),by=Time]'
                          '\n{0} <- data.table(rollmean({0}, 150))'
                          '\n{0}'.format(system, lipid))
                print(df)
               # if lipid:
               #     df = ro.r('{0}[Lipid_type == {1}]'
               #         '\n{0}[,Lipid_type:=NULL]'
               #         '\n{0} <- data.table(rollmean({0}, 100))'
               #         '\n{0}'.format(system, lipid))
               #     print(df)
               # else:
               #     raise NotImplementedError('Specify a lipid please')
            else:
                print("File {} does not exist".format(file_to_read))
                return None
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dat = merge_data('dat', 'temperature', *framenames)
        print(dat)
        plot = ggplot2.ggplot(dat)\
                + ggplot2.aes_string(x='Time', y='{}/Ntot'.format(neighbor), color='temperature')\
                + ggplot2.geom_point()\
                + ggplot2.guides(color=ggplot2.guide_legend('T / K'))\
                + add_axis('x', 'Time / ns', breaks=[0, 2600, 500])\
                + add_axis('y', 'N_{}/Ntot'.format(neighbor), breaks=[0.3, 0.55, 0.05])\
                + ggplot2.ggtitle('Segregation process in {}'.format(systemname))\
                + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(outputfilename)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

class RDF():
    def my_rdfplot(self, system, temperature, particlespair):
        ''' Plots data generated by selfimplementation of radial distribution function '''
        ncholranges = {'dppc_chol20':2,
                       'dppc_chol10':2,
                       'dppc_chol30':3,
                       'dupc_chol20':2,
                       'dppc_dupc_chol':2,
                       'dppc_dupc_chol05':1,
                       'dppc_dupc_chol40':4,
                        }
        temp = temperature
        print("PLOTTING:", system, temperature, particlespair)
        if isinstance(temperature, list):
            tframes = []
            for T in temperature:
                print("T is ", T)
                outputfile_name = '{}_rdf_allT_{}.pdf'.format(system, particlespair)
                #frame_regex = re.compile(r'(.*)(\..*)$')
                datafile_name = '{0}_{1}/rdf_{2}.dat'.format(system, T, particlespair)
                print("READING", datafile_name)
                dt = ro.r('dtcomplete <- fread("{0}")'\
                                              #'\ndtcomplete$temperature <- {1}'
                                              .format(datafile_name,))
                print(dt)
                framespecifiers = [['dtcomplete', 'complete'], ]
                if 'chol' in [i[:4] for i in system.split('_')] and ncholranges[system] is not None and particlespair != 'CHOL_CHOL':
                    for nchol in range(ncholranges[system]+1):
                        datafile_name = '{0}_{1}/rdf_{2}_{3}.dat'.format(system, T, particlespair, nchol)
                        print("READING", datafile_name)
                        dt = ro.r('dt{0} <- fread("{1}")'\
                                                 #'\ndt{0}$temperature <- {2}'
                                                 .format(nchol, datafile_name))
                        print(dt)
                        framespecifiers.append(['dt'+str(nchol), nchol])
                    dt = merge_data('dt'+str(T), 'NChol', *framespecifiers)
                    print("MERGED", dt)
                else:
                    dt = ro.r('dt{} <- dtcomplete'.format(T))
                tframes.append(['dt'+str(T), T])
            dt = merge_data('dt', 'temperature', *tframes)
            print("MERGED TEMPS", dt)
        else:
            outputfile_name = '{0}_rdf_{1}_{2}.pdf'.format(system, temp, particlespair)
            #frame_regex = re.compile(r'(.*)(\..*)$')
            datafile_name = '{0}_{1}/rdf_{2}.dat'.format(system, temp, particlespair)
            print("READING", datafile_name)
            dt = ro.r('dtcomplete <- fread("{0}")'.format(datafile_name)
                                          )
            print(dt)
            framespecifiers = [['dtcomplete', 'complete']]
            if 'chol' in [i[:4] for i in system.split('_')] and ncholranges[system] is not None and particlespair != 'CHOL_CHOL':
                for nchol in range(ncholranges[system]+1):
                    datafile_name = '{0}_{1}/rdf_{2}_{3}.dat'.format(system, temp, particlespair, nchol)
                    print("READING", datafile_name)
                    dt = ro.r('dt{0} <- fread("{1}")'.format(nchol, datafile_name))
                    framespecifiers.append(['dt'+str(nchol), nchol])
                dt = merge_data('dt', 'NChol', *framespecifiers)
        if 'chol' in [i[:4] for i in system.split('_')] and ncholranges[system] is not None and particlespair != 'CHOL_CHOL':
            aes = ggplot2.aes_string(x='V1', y='V2', color='NChol')
        else:
            aes = ggplot2.aes_string(x='V1', y='V2')
        plot = ggplot2.ggplot(dt)\
                + aes\
                + ggplot2.geom_point(alpha='0.4')\
                + add_axis('x', 'distance r / nm', breaks=[0.0, 2.5, 0.2])\
                + add_axis('y', 'g(r)', breaks=[0.0, 4, 0.5])\
                + ggplot2.theme_light() + my_theme()
                #===================================================================
                # + ggplot2.geom_line()\
                #===================================================================
                #+ ggplot2.geom_smooth(se='False', span='0.01', n=900)
        if isinstance(temperature, list):
            plot += ggplot2.facet_grid('temperature ~ .')
        gr.pdf(file=outputfile_name, width=10, height=6)
        print(plot)
        #save_plot(outputfile_name, plot, width=25, height=15)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')


    #def rdfplot(datafile, outputfile_name, ncholrange=None, temperature=None):
    def rdfplot(self, system, temperatures, interaction, lipid=None, ncholrange=None):
        ''' creates an rdf-plot
            datafilename path: <system>_<temperature>/datafiles/rdf/rdf_<inter1>_<inter2>.xvg
            system:         systemname
            temperature:    can be either a list with [Tstart, Tend] or an integer
            interaction:    tuple ('HEAD', 'HEAD') or ('HEAD', 'TAIL')
            ncholrange:     range for datafiles with nchol specification can be given with [start, end, step]
            allsys:         list of systems
        '''
        print("inputparameters:", system, temperatures, interaction, lipid, ncholrange)
        if isinstance(temperatures, list):
            tempframes = []
            for temp in temperatures:
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
                                                         #'\ndt{2}{0}$V2 <- dt{2}{0}$V2/max(dt{2}{0}$V2)'
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
                                                         #'\ndt{2}{0}$V2 <- dt{2}{0}$V2/max(dt{2}{0}$V2)'
                                                         '\ndt{2}{0}$temp <- {2}'
                                                         '\ndt{2}{0}'.format(syst, datafile_name, temp))
                        print("dataframe:", dataframes[(temp, syst)])
                        framespecifiers.append(['dt'+str(temp)+str(syst), syst])
                    dt = merge_data('dt'+str(system)+str(temp), 'System', *framespecifiers)
                    print("THIS IS DT at", temp, "\n", dt)
                    tempframes.append(['dt'+str(system)+str(temp), temp])
                    aes = ggplot2.aes_string(x='V1', y='V2', color='System')
                else:
                    if lipid is None:
                        lipid = ''
                    outputfile_name = 'rdf_{0}_tavg_{2}_{3}{4}.pdf'.format(system, temp, lipid, *interaction)
                    datafile_name = '{0}_{1}/datafiles/rdf/rdf_{2}-{3}.xvg'\
                                    .format(system, temp, interaction[0], interaction[1])
                    dataframes = {}
                    dataframes[temp] = ro.r('dt{0} <- fread("{1}")'
                                            #'\ndt{0}$V2 <- dt{0}$V2/max(dt{0}$V2)'
                                            '\ndt{0}$temp <- {2}'.format(temp, datafile_name, temp))
                    tempframes.append(['dt'+str(temp), temp])
                    aes = ggplot2.aes_string(x='V1', y='V2')
            dt = merge_data('dt_tavg', 'Temperature', *tempframes)
            print(dt)
        else:
            temp = temperatures
            if ncholrange is not None:
                outputfile_name = '{0}_rdf_{1}_{2}_{3}{4}_nchol.pdf'.format(system, temp, lipid, *interaction)
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
                                              #'\ndt{0}$V2 <- dt{0}$V2/max(dt{0}$V2)'
                                              '\ndt{0}'.format(nchol, datafile_name))
                    framespecifiers.append(['dt'+str(nchol), nchol])
                dt = merge_data('dt', 'NChol', *framespecifiers)
                aes = ggplot2.aes_string(x='V1', y='V2', color='NChol')
            else:
                datafile_name = '{0}_{1}/datafiles/rdf/rdf_{2}-{3}.xvg'\
                                    .format(system, temp, interaction[0], interaction[1])
                outputfile_name = 'rdf_{0}_{1}_{2}{3}.pdf'.format(system, temp, *interaction)
                dt = ro.r('{0} <- fread("{1}")'
                              '\n{0}<-{0}[V1>0.1]'
                              #'\n{0}$V2 <- {0}$V2/max({0}$V2)'
                              '\n{0}'.format('dt_raw', datafile_name))
                aes = ggplot2.aes_string(x='V1', y='V2')
        plot = ggplot2.ggplot(dt)\
                + aes\
                + ggplot2.geom_point(alpha='0.6')\
                + add_axis('x', 'distance r / nm', breaks=[0.0, 3.5, 0.5], limits=[0.0, 3.7])\
                + ggplot2.scale_y_continuous(name='g(r)')\
                + ggplot2.theme_light() + my_theme()\
    #            + ggplot2.geom_smooth(se='False', span='0.01', n=900)\
        if isinstance(temperatures, list):
            plot += ggplot2.facet_grid('Temperature ~ .')
        gr.pdf(file=outputfile_name, width=10, height=6)
        print(plot)
        #save_plot(outputfile_name, plot, width=25, height=15)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')


class autocorrelation():
    ''' Calculate and plot the autocorrelation function of a time series '''
    def __init__(self, systemname, object_with_order, orderparameter_name, lipidtype, temperature='high'):
        '''
        object_with_order: Colum name of object of interest (e.g. Residue or Pair)
        orderparameter_name: Column name of order parameter (e.g. Scd, AvgOrd, ...)
        lipidtype: Type of lipid
        '''

        self.systemname = systemname
        self.object_with_order = object_with_order
        self.orderparameter_name = orderparameter_name
        self.lipidtype = lipidtype
        self.temperature_set = temperature

    def collect_datatables(self):
        '''  '''
        #systems = [self.systemname+str(T) for T in Tranges[self.systemname]]
        if self.temperature_set == 'whole':
            temperatures = Tranges_whole[self.systemname]
        elif self.temperature_set == 'high':
            temperatures = Tranges_high[self.systemname]
        framenames = []
        for temp in temperatures:
            print("On temperature", temp)
            datafilename = self.systemname+'_'+str(temp)+'/scd_distribution.dat'
            self.calc_acf_for_T('df'+str(temp), datafilename)
            framenames.append(('df'+str(temp), temp))
        dt = merge_data('final', 'temperature', *framenames)
        print("binned df:", dt)
        return dt

    def calc_acf_for_T(self, tablename, datafilename):
        read_datafile(tablename, datafilename)
        ro.r('total <- data.table(lag=numeric(), acf=numeric())'
                 '\n{3} = {3}[Type=="{1}"][Time>max({3}$Time)*0.8]'
                 #'\n{3} = {3}[Type=="{1}"][Time>1000000]'
                 '\nfor(i in min({3}${0}):max({3}${0})){{'
                 '\nacf <- acf({3}[{0}==i]${2}, lag.max={3}[{0}==i][,.N], plot=FALSE)'
                 #'\nacf <- acf({3}[{0}==i]${2}, lag.max=150)'
                 '\ntotal <- rbind(total, list(acf$lag, acf$acf))'
                 '\n}}'
                 '\ntotal'
                 .format(self.object_with_order, self.lipidtype, self.orderparameter_name, tablename,))
        dat = ro.r('{} <- total[,lapply(.SD, mean), by=lag]'.format(tablename))
        print("Systems acf\n", dat)
        return dat

    def plot_acf(self, filename):
        ''' Function for plotting the calculated acf data '''
        dat = self.collect_datatables()
        print('final', dat)
        aes = ggplot2.aes_string(x='lag', y='acf', color='factor(temperature)')
        guide = ggplot2.guides(color=ggplot2.guide_legend('T / K'))
        plot = ggplot2.ggplot(dat)\
                + aes\
                + guide\
                + ggplot2.scale_x_continuous("lag / ns")\
                + add_axis('y', "ACF", limits=[-0.1, 0.55], breaks=[-0.1, 0.5, 0.1])\
                + ggplot2.ggtitle(self.systemname)\
                + ggplot2.geom_point(alpha=0.7)\
                + ggplot2.theme_light() + my_theme()
                #+ ggplot2.geom_smooth(ggplot2.aes_string(linetype='factor(temperature)'), se=False)\
        print(plot)
        #time.sleep(5)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')


class Cholesterol_inputplots():
    def __init__(self, systemnames=Tranges_high.keys(), averaging=True, temperatures='whole'):
        self.lipidranges = {'DPPC':(3, 6, 1), 'CHL1':(2, 5, 0.5), 'DUPC':(3, 8, 1)}
        #self.enthalpyranges = {}
        self.temperature_set = temperatures
        self.systemnames = systemnames
        self.temperatures = {}
        for syst in systemnames:
            if self.temperature_set == 'whole':
                self.temperatures[syst] = Tranges_whole[syst]
            elif self.temperature_set == 'high':
                self.temperatures[syst] = Tranges_high[syst]
        self.enthalpyranges = {}
        self.systemnames = systemnames

    def plot_H_cholchol(self, systemname, filename, averaging=True):
        ''' H(nchol) Chol-Chol interaction plots '''
        if isinstance(systemname, list):
            tablelist = []
            for syst in systemname:
                tavg = self.temperature_avg(syst, syst, datatype='H', averaging=averaging)
                print(tavg)
                tablelist.append((syst, syst))
            tavg = merge_data('final_table', 'System', *tablelist)
            aes =  ggplot2.aes_string(x='Chol', y='Ntot', color='System')
            #yname ="N CHOL neighbors"
            #yname = 'N total'
        else:
            tavg = self.temperature_avg(systemname, systemname, datatype='H', averaging=averaging)
            print("TAVG", tavg)
        if tavg is None:
            return
        if averaging:
            aes = ggplot2.aes_string(x='NChol', y='Etot')
            guide = ggplot2.guides(size=False)
        else:
            aes = ggplot2.aes_string(x='NChol', y='Etot', color='factor(temperature)')
            guide = ggplot2.guides(size=False, color=ggplot2.guide_legend('T'))
        if isinstance(systemname, list):
            aes = ggplot2.aes_string(x='NChol', y='Etot', color='factor(System)')
            guide = ggplot2.guides(size=False, color=ggplot2.guide_legend('System'))
        plot = ggplot2.ggplot(tavg)\
                + aes\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + add_axis('x', "Number of cholesterol neighbors", breaks=[0, 4, 1], limits=[-0.2, 4.2])\
                + add_axis('y', 'enthalpy', breaks=[-30, -10, 5])\
                + ggplot2.geom_smooth(se=False, method='lm', formula='y ~ x ')\
                + guide\
                + ggplot2.ggtitle('CHOL-CHOL interaction')\
                + ggplot2.scale_size("weight")\
                + ggplot2.theme_light() + my_theme()
                #===============================================================
                # + add_axis('y', 'enthalpy', breaks=[-24, -14, 2])\
                #===============================================================
        #print(plot)
        #=======================================================================
        # xaxis = add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair='CHL1_CHL1')
        #=======================================================================
        #time.sleep(5)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
    def fit_H(self, systemname, filename, fit_degree=1):
        ''' 
            Outputs fit coefficient summary for CHOL-CHOL interaction
            and plot with respective fit
            
        '''
        tavg = self.temperature_avg(systemname, systemname, datatype='H')
        print("TAVG", tavg)
        if len(tavg[1]) == 0:
            return
        fit = ro.r('fit <- lm({0}$Etot ~ poly({0}$NChol, degree={1}, raw=TRUE))'
                   '\nsummary(fit)'.format(systemname, fit_degree)
                   )
        with open('{}_H_{}_deg{}_polyfit.dat'.format(systemname, 'CHL1_CHL1', fit_degree), "w") as fitf:
            print(fit, file=fitf)
        plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x='NChol', y='Etot')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth(method='lm', se=False, fullrange=True, formula='y~poly(x, degree={})'.format(fit_degree))\
                + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(filename)
        print(plot)
        gr.dev_off()
    def fit_N(self, systemname, filename, fit_degree=1):
        ''' 
            Outputs fit coefficient summary
            for N_CHOL and plot with respective fit 
        '''
        data = self.temperature_avg(systemname, systemname, datatype='N', averaging=False)
        print("TAVG", data)
        if len(data[1]) == 0:
            return
        data = ro.r('{0} <- {0}[,lapply(.SD, mean), by=list(Chol, temperature)]'.format(systemname))
        fits = []
        for Nc in range(5):
            fit = ro.r('{0}$temperature <- as.numeric({0}$temperature)'
                       '\n{0}{2} <- {0}[Chol == {2}]' 
                       '\nfit <- lm({0}{2}$Ntot ~ poly({0}{2}$temperature, degree={1}, raw=TRUE))'
                       '\nsummary(fit)'.format(systemname, fit_degree, Nc)
                       )
            fits.append(fit)
            with open('{}_N_{}_Nc{}_deg{}_polyfit.dat'.format(systemname, 'CHL1', Nc, fit_degree,), "w") as fitf:
                print(fit, file=fitf)
        plot = ggplot2.ggplot(data)\
                + ggplot2.aes_string(x='as.numeric(temperature)', y='Ntot', color='factor(Chol)')\
                + ggplot2.geom_point()\
                + ggplot2.guides(color=ggplot2.guide_legend('NChol'))\
                + ggplot2.geom_smooth(method='lm', se=False, fullrange=True, formula='y~poly(x, degree={})'.format(fit_degree))\
                + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(filename)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
    def plot_N_chol(self, systemname, filename, neib='Ntot', averaging=False):
        ''' N(nchol) for host cholesterol '''
        if isinstance(systemname, list):
            tablelist = []
            for syst in systemname:
                tavg = self.temperature_avg(syst, syst, datatype='N', averaging=averaging)
                print(tavg)
                tablelist.append((syst, syst))
            tavg = merge_data('final_table', 'System', *tablelist)
            aes =  ggplot2.aes_string(x='Chol', y=neib, color='factor(temperature)', shape='System')
            guide = ggplot2.guides(color=ggplot2.guide_legend("T / K"))
            #yname ="N CHOL neighbors"
        else:
            tavg = self.temperature_avg(systemname, systemname, datatype='N', averaging=averaging)
            print(tavg)
            aes =  ggplot2.aes_string(x='Chol', y=neib, color="factor(temperature)")
            guide = ggplot2.guides(color=ggplot2.guide_legend("T / K"))
        if tavg is None:
            return
        if neib == 'Ntot':
            yname = 'N'
            axis = add_axis('y', yname, breaks=[3.0, 7.5, 0.5], limits=[2.6, 7.7])
        else:
            yname = '{} neighbors'.format(neib)
            #axis = add_axis('y', yname, breaks=[3.0, 4.0, 0.5], limits=[2.6, 4.2])
            axis = add_axis('y', yname, breaks=[3.0, 7.5, 0.5], limits=[2.6, 7.7])
        #guide =  ggplot2.guides()
        plot = ggplot2.ggplot(tavg)\
            + aes\
            + ggplot2.scale_x_continuous('CHOL neighbors')\
            + axis\
            + guide\
            + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.8)\
            + ggplot2.guides(size=False)\
            + ggplot2.scale_size()\
            + ggplot2.ggtitle('Neighbors of cholesterol')\
            + ggplot2.theme_light() + my_theme()
            #===================================================================
            # + add_axis('x', 'scd', breaks='auto')\
            #===================================================================
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    def temperature_avg(self, outputtable, systemname, datatype, averaging=True):
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
        print("INPUT:", outputtable, systemname, 'Chol')
        framenames = []
        for temp in self.temperatures[systemname]:
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            print(system)
            if datatype == 'N':
                file_to_read = '{}/Nofscd.dat'.format(system)
            elif datatype == 'H':
                file_to_read = '{}/EofscdCHL1_CHL1.dat'.format(system)
            if os.path.exists(file_to_read):
                df = read_datafile(system, file_to_read)
                print(df)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
            #Averaging before everything's collected=============
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dt = merge_data(outputtable, 'temperature', *framenames)
        print("MERGED DATA:", dt)
        dt = ro.r('{0}$temperature <- as.numeric({0}$temperature)'.format(systemname))
        #=======================================================================
        # if averaging == False:
        #     print("NO AVERAGE")
        #     dt = ro.r('{0} <- {0}[Lipid_type=="{1}"]\n'.format(outputtable, 'CHL1'))
        #     print(dt)
        #     return dt
        #=======================================================================
        if datatype == 'N':
            if not averaging:
                df_final = ro.r(
                                '\n{0} <- {0}[Lipid_type=="{1}"]\n'
                                '\n{0} <- {0}[,Lipid_type:=NULL]\n'
                                '\nweightdata <- {0}'
                                '\n{0} <- {0}[,lapply(.SD,mean),by=list(Chol, temperature)]'
                                '\n{0}$weight  <- weightdata[,.N,by=list(Chol, temperature)][,N]'
                                '\n{0}$weight <- {0}$weight/sum({0}$weight)'
                                '\ninsignificant <- 0.01*max({0}$weight)\n'
                                '\n{0} <- {0}[weight>insignificant][Chol<5]\n'
                                '\n{0}'
                                .format(outputtable, 'CHL1'))
            else:
                df_final = ro.r('\n{0} <- subset({0},select=-c(temperature))\n'
                                '\n{0} <- {0}[Lipid_type=="{1}"]\n'
                                '\n{0} <- {0}[,Lipid_type:=NULL]\n'
                                '\nweightdata <- {0}'
                                '\n{0} <- {0}[,lapply(.SD,mean,),by=Chol]'
                                '\n{0}$weight  <- weightdata[,.N,by=Chol][,N]'
                                '\n{0}$weight <- {0}$weight/sum({0}$weight)'
                                '\ninsignificant <- 0.01*max({0}$weight)\n'
                                '\n{0} <- {0}[weight>insignificant][Chol<5]\n'
                                '\n{0}'
                                .format(outputtable, 'CHL1'))
        elif datatype == 'H':
            df_final = ro.r('\n{0} <- subset({0},select=-c(temperature))\n'
                            '\nweightdata <- {0}'
                            '\n{0} <- {0}[,lapply(.SD,mean,),by=NChol]'
                            '\n{0}$weight  <- weightdata[,.N,by=NChol][,N]'
                            '\n{0}$weight <- {0}$weight/sum({0}$weight)'
                            '\ninsignificant <- 0.001*max({0}$weight)\n'
                            '\n{0} <- {0}[weight>insignificant]\n'
                            '\n{0}'
                            .format(outputtable))
        print("FINALFRAME:", df_final)
        return df_final
    
    #===========================================================================
    # def bin_raw_data(self, outputtablename, datatablename, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05):
    #     ''' Discretize data by aggregating specified column '''
    #     print("Binning... {}".format(datatablename))
    #     labelbreakh = breakh-breakdx
    #     mydf = ro.r(\
    #         '\n{0}$bin{1} <- cut({0}${1},'
    #                         'breaks=seq({2}, {3}, {4}),'
    #                         'labels=seq({2}, {5}, {4}))'
    #         '\nbckp_data <- {0}'
    #         '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=bin{1}]'
    #         '\n{6}$weight <- bckp_data[,.N,by=bin{1}][,N]'
    #         '\n{6}$weight <- {6}$weight/sum({6}$weight)'
    #         '\n{6}'\
    #         .format(datatablename, column,\
    #                 breakl, breakh, breakdx, labelbreakh,\
    #                 outputtablename)\
    #         )
    #     print("Finished binning")
    #     return mydf
    #===========================================================================


class PofScdplots:
    def __init__(self, systemnames=Tranges_whole.keys(), temperature='high'):
        self.lipidranges = {'DPPC':(3, 6, 1), 'CHL1':(2, 5, 0.5), 'DUPC':(3, 8, 1)}
        #self.enthalpyranges = {}
        self.systemnames = systemnames
        self.temperature_set = temperature
        self.temperatures = {}
        for syst in systemnames:
            if self.temperature_set == 'whole':
                self.temperatures[syst] = Tranges_whole[syst]
            elif self.temperature_set == 'high':
                self.temperatures[syst] = Tranges_high[syst]
    def fitfunction(self, systemname, lipid, filename, fit_degree, temperature, binsize=0.05):
        ''' Outputs fit functions '''
        self.read_data(systemname, systemname, lipid, temperature=temperature)
        dat = ro.r('{0} = {0}[Time>0.8*max({0}$Time)][Type=="{2}"]'
                   '\n{0}$bin = cut({0}$Scd, breaks=seq(-0.5, 1.0, 0.05), labels=seq(-0.5, 0.95, {1}))'
                   '\n{0}[,Type:=NULL]'
                   '\nbckp = {0}'
                   '\ndat2 = {0}[, lapply(.SD, mean), by=bin]'
                   '\ndat2$freq = bckp[, .N, by=bin][,N]'
                   '\ndat2$freq = dat2$freq/(sum(dat2$freq)*{1})'
                   '\ndat2$lnfreq = log(dat2$freq)'
                   '\ndat2'.format(systemname, binsize, lipid)
                   )
        if len(dat[1]) == 0:
            return
        fit = ro.r('fit <- lm(dat2$freq ~ log(poly(dat2$Scd, degree={0}, raw=TRUE)))'
                   '\nsummary(fit)'.format(fit_degree)
                   )
        with open('{}_{}_{}_deg{}_polyfit.dat'.format(systemname, temperature, lipid, fit_degree), "w") as fitf:
            print(fit, file=fitf)
        plot = ggplot2.ggplot(dat)\
                + ggplot2.aes_string(x='Scd', y='freq')\
                + ggplot2.geom_point(color='Red')\
                + ggplot2.geom_point(aes_string(y='lnfreq'), color='Green')\
                + add_axis('x', 'scd', breaks=[-0.3, 1.0, 0.2])\
                + ggplot2.geom_smooth(method='lm', se=False, fullrange=True, formula='log(y)~poly(x, degree={})'.format(fit_degree))\
                + ggplot2.geom_smooth(method='lm', se=False, fullrange=True, formula='y~poly(x, degree={})'.format(fit_degree))\
                + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(filename)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
        
    def plot_standard(self, systemname, lipid, filename, temperature=None, timerange=None):
        ''' '''
        if timerange is not None and temperature is None:
            raise ValueError('If timerange is specified, specify a temperature as well.')
        dat = self.read_data(systemname, systemname, lipid, temperature=temperature)
        if timerange is not None:
            framenames = []
            for time in range(timerange[0], timerange[1], 100):
                dt = ro.r('{0}{1} <- {0}[Time>{2}][Time<{3}]'.format(systemname, time, time*1000, (time+100)*1000))
                print("TIMES:{} - {}\n".format(time, time+100), dt)
                framenames.append((systemname+str(time), str(time)))
            dat = merge_data(systemname, 'timeranges', *framenames)
            aes = ggplot2.aes_string(x='Scd', fill='timeranges')
            guide =  ggplot2.guides(fill=ggplot2.guide_legend("[time,+100ns["))
        else:
            aes =  ggplot2.aes_string(x='Scd', fill='temperature')
            guide =  ggplot2.guides(color=ggplot2.guide_legend("Temperature"))
        plot = ggplot2.ggplot(dat)\
            + aes\
            + ggplot2.geom_density(alpha=0.5)\
            + add_axis('x', 'scd', breaks=[-0.5, 1.0, 0.2])\
            + ggplot2.scale_y_continuous(name='frequency')\
            + ggplot2.scale_size()\
            + guide\
            + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
    def read_data(self, outputtable, systemname, lipid, temperature=None):
        print("INPUT:", outputtable, systemname, lipid)
        framenames = []
        if temperature is None:
            for temp in self.temperatures[systemname]:
                print("At temperature:", temp)
                system = '{}_{}'.format(systemname, temp)
                print(system)
                file_to_read = '{}/scd_distribution.dat'.format(system)
                if os.path.exists(file_to_read):
                    df = read_datafile(system, file_to_read)
                    print(df)
                else:
                    print("File {} does not exist".format(file_to_read))
                    return None
                #Averaging before everything's collected=============
                framenames.append(system)
            framenames = [(i, i[-3:]) for i in framenames]
            print("Framenames:", framenames)
            merge_data(outputtable, 'temperature', *framenames)
        elif temperature is not None:
            system = '{}_{}'.format(systemname, temperature)
            file_to_read = '{}/scd_distribution.dat'.format(system)
            read_datafile(outputtable, file_to_read)
            ro.r('{0}$temperature <- {1}'.format(outputtable, temperature))
        dt = ro.r('{0} <- {0}[Type=="{1}"]\n'.format(outputtable, lipid))
        print("LAST", dt)
        return dt
    def avg_scd(self, outputfilename, lipid, systems):
        ''' Calculate rolling average of scd values over time for systems and all temperatures '''
        frames_syst = []
        if not isinstance(systems, list):
            systems = [systems,]
        print("INPUT:", outputfilename, lipid, systems)
        for syst in systems:
            print("CALCULATING FOR:", syst)
            frames_temp = []
            for temp in self.temperatures[syst]:
                print("At temperature:", temp)
                system = '{}_{}'.format(syst, temp)
                print(system)
                file_to_read = '{}/scd_distribution.dat'.format(system)
                if os.path.exists(file_to_read):
                    df = read_datafile(system, file_to_read)
                    ro.r('{0} <- {0}[Type == "{1}"]'
                         '\n{0}[, Residue:=NULL]'
                         '\n{0}[, Type:=NULL]'
                         '\n{0}$Time = {0}$Time/1000'
                         '\n{0} <- {0}[, lapply(.SD, mean), by=Time]'
                         .format(system, lipid))
                    print(df)
                else:
                    print("File {} does not exist".format(file_to_read))
                    return None
                #Averaging before everything's collected=============
                frames_temp.append(system)
            frames_temp = [(i, i[-3:]) for i in frames_temp]
            print("Framenames:", frames_temp)
            merge_data(syst, 'temperature', *frames_temp)
            frames_syst.append((syst, syst))
        dt = merge_data('final', 'system', *frames_syst)
        print(dt)
        if not isinstance(systems, list) or len(systems) == 1:
            aes = ggplot2.aes_string(x="Time", y="Scd", color="temperature")
            guide = ggplot2.guides(color=ggplot2.guide_legend("T / K"))
        else:
            aes = ggplot2.aes_string(x="Time", y="Scd", color="temperature", shape="system", linetype="system")
            guide = ggplot2.guides(color=ggplot2.guide_legend("T / K"), shape=ggplot2.guide_legend("system"))
        plot = ggplot2.ggplot(dt)\
            + aes\
            + ggplot2.geom_point(alpha=0.4)\
            + ggplot2.geom_smooth(se=False, size=2)\
            + ggplot2.scale_x_continuous(name='time / ns')\
            + ggplot2.scale_y_continuous(name='average order')\
            + ggplot2.ggtitle(lipid)\
            + ggplot2.scale_size()\
            + guide\
            + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(file=outputfilename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
        
class NofScdplots():
    def __init__(self, systemnames=Tranges_whole.keys(), temperature='whole'):
        self.lipidranges = {'DPPC':(3, 6, 1), 'CHL1':(2, 5, 0.5), 'DUPC':(3, 8, 1)}
        #self.enthalpyranges = {}
        self.systemnames = systemnames
        self.temperature_set = temperature
        self.temperatures = {}
        for syst in systemnames:
            if self.temperature_set == 'whole':
                self.temperatures[syst] = Tranges_whole[syst]
            elif self.temperature_set == 'high':
                self.temperatures[syst] = Tranges_high[syst]

    def fit(self, systemname, lipid, fit_degree=1, nchol=None):
        ''' Returns fit parameters '''
        dat = self.temperature_avg(systemname, systemname, lipid, averaging=False)
        if nchol is not None:
            for Nc in range(0, 5):
                filename = '{}_nofscd_c{}_{}_deg{}fit.pdf'.format(systemname, Nc, lipid, fit_degree)
                fit = ro.r('dat{2} <- {0}[Chol=={2}]'
                           'fit <- lm(Ntot ~ poly(Host_Scd, as.numeric(temperature), degree={1}, raw=TRUE), data=dat{2})'
                           '\nsummary(fit)'.format(systemname, fit_degree, Nc))
                with open("{}_nofscd_chol{}_{}_deg{}_fit.dat".format(systemname, Nc, lipid, fit_degree), "w") as fitf:
                    print(fit, file=fitf)
                plot = ggplot2.ggplot(dat)\
                    + ggplot2.aes_string(x='Host_Scd', y='Ntot', z='temperature', color="factor(Chol")\
                    + ggplot2.geom_point()\
                    + ggplot2.geom_smooth(se=False, fullrange=True, method='lm', formula='y ~ poly(x, z, degree={}, raw=True'.format(fit_degree))\
                    + add_axis('x', 'scd', breaks=[-0.3, 1.0, 0.2])\
                    + ggplot2.ggtitle('Neighbors of '+ str(lipid))\
                    + ggplot2.theme_light() + my_theme()
                    #+ add_axis('y', yname, breaks=[3.0, 4.5, 0.5])\
                gr.pdf(file=filename, width=10, height=7)
                print(plot)
                gr.dev_off()
                ro.r('rm(list=ls(all=TRUE))')
        else:
            filename = '{}_nofscd_{}_deg{}fit.pdf'.format(systemname, lipid, fit_degree)
            fit = ro.r(
                       'fit <- lm(Ntot ~ poly(Host_Scd, as.numeric(temperature), degree={1}, raw=TRUE), data=dat{0})'
                       '\nsummary(fit)'.format(systemname, fit_degree))
            with open("{}_nofscd_{}_deg{}_fit.dat".format(systemname, lipid, fit_degree), "w") as fitf:
                print(fit, file=fitf)
            plot = ggplot2.ggplot(dat)\
                + ggplot2.aes_string(x='Host_Scd', y='Ntot', z='temperature')\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth(se=False, fullrange=True, method='lm', formula='y ~ poly(x, z, degree={}, raw=True'.format(fit_degree))\
                + add_axis('x', 'scd', breaks=[-0.3, 1.0, 0.2])\
                + ggplot2.ggtitle('Neighbors of '+ str(lipid))\
                + ggplot2.theme_light() + my_theme()
                #+ add_axis('y', yname, breaks=[3.0, 4.5, 0.5])\
            gr.pdf(file=filename, width=10, height=7)
            print(plot)
            gr.dev_off()

    def plot_standard(self, systemname, lipid, filename, y_type='Ntot', neib_spec=None, averaging=True, plotrel=False):
        yaxis = None
        tavg = self.temperature_avg(systemname, systemname, lipid, neib_spec=neib_spec, averaging=averaging, plotrel=plotrel)
        if tavg is None:
            return
        print(tavg)
        if lipid =='CHL1':
            lipid = 'Chol'
        if y_type != 'Ntot':
            yname = 'N {}'.format(y_type)
        else:
            if isinstance(lipid, list):
                yname = 'N'
            else:
                yname = 'N'
        #if isinstance(neib_spec, list):
        #    aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, )
        #    guide =  ggplot2.guides(color=ggplot2.guide_legend(neib_spec), shape=ggplot2.guide_legend("Lipid"))
        if neib_spec is None:
            if isinstance(lipid, list):
                aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, color='Lipid_type')
                guide =  ggplot2.guides(color=ggplot2.guide_legend("Lipid"), size=False)
            elif averaging == False:
                #aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, color='temperature')
                aes =  ggplot2.aes_string(x='as.numeric(temperature)', y=y_type, )
                guide =  ggplot2.guides(color=ggplot2.guide_legend("Temperature"), size=False)
            else:
                neib_spec = ''
                aes =  ggplot2.aes_string(x='Host_Scd', y=y_type)
                guide = ggplot2.guides(shape=ggplot2.guide_legend("Dummy"), size=False)
        elif neib_spec is not None:
            if averaging == False:
                #aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, color='factor('+neib_spec+')',)# shape='temperature')
                aes =  ggplot2.aes_string(x='as.numeric(temperature)', y=y_type, color='factor('+neib_spec+')',)# shape='temperature')
                guide =  ggplot2.guides(shape=ggplot2.guide_legend("Temperature"), color=ggplot2.guide_legend('N'+neib_spec), size=False)
            else:
                aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, color='factor('+neib_spec+')')
                guide =  ggplot2.guides(color=ggplot2.guide_legend('N'+neib_spec), size=False)
                neib_spec = ''
        #elif isinstance(lipid, list) and neib_spec is not None:
        #    neib_spec = ''
        #    aes =  ggplot2.aes_string(x='Host_Scd', y=y_type, shape='Lipid_type', color='factor('+neib_spec+')')
        #    guide =  ggplot2.guides(color=ggplot2.guide_legend(neib_spec), shape=ggplot2.guide_legend("Lipid"))
        #dupcfunc = ggplot2.stat_function(ggplot2.aes_string(linetype='"DLiPC"'), fun=ro.r("function(x){3.52807504+0.12268938*x-0.09141714*x**2}"), color='black', size=2)
        #dppcfunc = ggplot2.stat_function(ggplot2.aes_string(linetype='"DPPC"'), fun=ro.r("function(x){(0.94123979/(1 + exp(-18.29302015 * (x - 0.63623191))) + 3.88006658)}"), color='black', size=2, )
        if lipid == 'DPPC' and (systemname == 'dppc_chol20' or systemname == 'dppc_chol10' or systemname == 'dppc_chol30'):
            if neib_spec == 'Chol':
                yaxis = add_axis('y', yname, breaks=[3, 10, 1])
            else:
                if y_type == 'Ntot':
                    yaxis = add_axis('y', yname, breaks=[3.5, 10, 1])
                else:
                    if y_type == 'DPPC':
                        yaxis = add_axis('y', yname, breaks=[3, 7, 1])
                    elif y_type == 'CHL1':
                        yaxis = add_axis('y', yname, breaks=[0, 2, 0.5])
        elif lipid == 'DUPC' and (systemname == 'dupc_chol20'):
            if neib_spec == 'Chol':
                yaxis = add_axis('y', yname, breaks=[3, 6.5, 1])
            else:
                if y_type == 'Ntot':
                    yaxis = add_axis('y', yname, breaks=[3, 7, 1])
                else:
                    if y_type == 'DUPC':
                        yaxis = add_axis('y', yname, breaks=[2.5 , 5, 1])
                    elif y_type == 'CHL1':
                        yaxis = add_axis('y', yname, breaks=[0, 2, 0.5])
        elif lipid == 'CHL1':
            yaxis = add_axis('y', yname, breaks=[2.5, 8.0, 1])
        plot = ggplot2.ggplot(tavg)\
            + aes\
            + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
            + ggplot2.geom_smooth(se=False, weight='weight')\
            + ggplot2.scale_size()\
            + ggplot2.ggtitle('Neighbors of '+ str(lipid))\
            + guide\
            + ggplot2.theme_light() + my_theme()
            #+ add_axis('x', 'scd', breaks=[-0.3, 1.0, 0.2])\
            #+ add_axis('y', yname, breaks=[3.0, 4.5, 0.5])\
        if yaxis is not None:
            plot += yaxis
        if averaging == False and neib_spec is not None:
            plot += ggplot2.scale_shape_manual(values=ro.r("1:"+str(len(self.temperatures[systemname]))))
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    def temperature_avg(self, outputtable, systemname, lipid, averaging=True, neib_spec=None, plotrel=False):
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
        print("INPUT:", outputtable, systemname, lipid)
        framenames = []
        for temp in self.temperatures[systemname]:
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            print(system)
            file_to_read = '{}/Nofscd.dat'.format(system)
            if os.path.exists(file_to_read):
                df = read_datafile(system, file_to_read)
                print(df)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
            #Averaging before everything's collected=============
            df = self.bin_raw_data(system, system, neib_spec=neib_spec, plotrel=plotrel)
            print("binned df:", df)
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dt = merge_data(outputtable, 'temperature', *framenames)
        print("MERGED DATA:", dt)
        if averaging == False:
            print("NO AVERAGE")
            if neib_spec is None:
                dt = ro.r('{0} <- {0}[Lipid_type=="{1}"]\n'.format(outputtable, lipid))
            else:
                dt = ro.r('{0} <- {0}[Lipid_type=="{1}"][Chol<5]\n'.format(outputtable, lipid))
            print(dt)
            return dt
        if isinstance(lipid, list):
            df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
                            '{0} <- {0}[,lapply(.SD, weighted.mean,w=weight),by=list(Lipid_type, binHost_Scd)]'
                            .format(outputtable))
        elif isinstance(neib_spec, list):
            df_final = ro.r('\n{0} <- subset({0},select=-c(temperature))\n'
                            '\n{0} <- {0}[Lipid_type=="{1}"]\n'
                            '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binHost_Scd, {2})]'
                            '\ninsignificant <- 0.01*max({0}$weight)\n'
                            '\n{0} <- {0}[,Lipid_type:=NULL][weight>insignificant][Chol<4]\n'
                            .format(outputtable, lipid, ','.join(neib_spec)))
        else:
            if neib_spec is None:
                print("NEIB_SPEC {}".format(neib_spec))
                df_final = ro.r('\n{0} <- subset({0},select=-c(temperature))\n'
                                '\n{0} <- {0}[Lipid_type=="{1}"]\n'
                                '\n{0} <- {0}[,Lipid_type:=NULL]\n'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binHost_Scd]'
                                '\ninsignificant <- 0.01*max({0}$weight)\n'
                                '\n{0} <- {0}[weight>insignificant][Chol<4]\n'
                                .format(outputtable, lipid))
            else:
                print("NEIB_SPEC {}".format(neib_spec))
                #if plotrel == False:
                df_final = ro.r('\n{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[Lipid_type=="{1}"]'
                                '\n{0} <- {0}[,Lipid_type:=NULL]'
                                #'{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binHost_Scd, {2})]'\
                                '\n{0} <- {0}[,lapply(.SD, weighted.mean, w=weight),by=list(binHost_Scd, {2})]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant][Chol<5]'
                                .format(outputtable, lipid, neib_spec))
                #===============================================================
                # elif plotrel == True:
                #     quotients = '/'.split(neib_spec)
                #     df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
                #                     '{0} <- {0}[Lipid_type=="{1}"]\n'
                #                     'insignificant <- 0.01*max({0}$weight)\n'
                #                     '{0} <- {0}[weight>insignificant]\n'
                #                     '{0} <- {0}[,Lipid_type:=NULL]\n'
                #                     '{0}${2} <- {0}[{3}]/{0}[{4}]\n'
                #                     '{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binHost_Scd]'\
                #                     .format(outputtable, lipid, neib_spec, quotients[0], quotients[1]))
                #===============================================================
        print("FINALFRAME:", df_final)
        if len(df_final[0]) == 0:
            print("DONT")
            df_final = None
        return df_final
    def bin_raw_data(self, outputtablename, datatablename, column='Host_Scd', neib_spec=None, plotrel=False, breakl=-0.5, breakh=1.0, breakdx=0.05):
        ''' Discretize data by aggregating specified column '''
        print("Binning... {} :".format(datatablename))
        labelbreakh = breakh-breakdx
        if neib_spec is None:
            mydf = ro.r(\
                '\n{0}$bin{1} <- cut({0}${1},'
                                'breaks=seq({2}, {3}, {4}),'
                                'labels=seq({2}, {5}, {4}))'
                '\nbckp_data <- {0}'
                '\n{6} <- {0}[,lapply(.SD, mean), by=list(bin{1}, Lipid_type)]'
                '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Lipid_type)][,N]'
                '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                '\n{6}'\
                .format(datatablename, column,
                        breakl, breakh, breakdx, labelbreakh,
                        outputtablename),)
        #=======================================================================
        # elif isinstance(neib_spec, list):
        #     mydf = ro.r(\
        #         '\n{0}$bin{1} <- cut({0}${1},'
        #                         'breaks=seq({2}, {3}, {4}),'
        #                         'labels=seq({2}, {5}, {4}))'
        #         '\nbckp_data <- {0}'
        #         '\n{6} <- {0}[,lapply(.SD, mean), by=list(bin{1}, Lipid_type, {7})]'
        #         '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Lipid_type, {7})][,N]'
        #         '\n{6}'\
        #         .format(datatablename, column,
        #                 breakl, breakh, breakdx, labelbreakh,
        #                 outputtablename, ','.join(neib_spec)),)
        #=======================================================================
        else:
            if plotrel == False:
                mydf = ro.r(\
                    '\n{0}$bin{1} <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))'
                    '\nbckp_data <- {0}'
                    '\n{6} <- {0}[,lapply(.SD, mean), by=list(bin{1}, Lipid_type, {7})]'
                    '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Lipid_type, {7})][,N]'
                    '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                    '\n{6}'\
                    .format(datatablename, column,
                            breakl, breakh, breakdx, labelbreakh,
                            outputtablename, neib_spec),)
            elif plotrel == True:
                quotients = neib_spec.split("_")
                print(quotients)
                mydf = ro.r(\
                    '\n{0}$bin{1} <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))'
                    '\n{0}${7} <- {0}${8}/{0}${9}'
                    '\n{0}${7} <- cut({0}${7},'
                                    'breaks=seq(-0.125, 1.125, 0.25),'
                                    'labels=seq(0., 1.0, 0.25))'
                    '\nbckp_data <- {0}'
                    '\n{6} <- {0}[,lapply(.SD, mean), by=list(bin{1}, Lipid_type, {7})]'
                    '\n{6}$weight <- bckp_data[,.N,by=list(bin{1}, Lipid_type, {7})][,N]'
                    '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                    '\n{6}'\
                    .format(datatablename, column,
                            breakl, breakh, breakdx, labelbreakh,
                            outputtablename, neib_spec, quotients[0], quotients[1])
                            )
        print(mydf)
        return mydf

class Selfinteraction():
    '''
    Plot and fit 
    '''
    def __init__(self, interaction, systemnames=Tranges_whole.keys(), averaging=True, temperature='whole'):
        self.interaction = interaction
        self.systemnames = systemnames
        self.temperature_set = temperature
        self.temperatures = {}
        for syst in systemnames:
            if self.temperature_set == 'whole':
                self.temperatures[syst] = Tranges_whole[syst]
            elif self.temperature_set == 'high':
                self.temperatures[syst] = Tranges_high[syst]
    def fit(self, systemname, lipid, filename, fit_degree):
        ''' '''
        tavg = self.temperature_avg(systemname, systemname, lipid, averaging=True)
        if len(tavg[1]) == 0:
            return
        #tavg = ro.r('\n{0}$Etot = {0}$Etot-max({0}$Etot)\n{0}'.format(systemname))
        fit = ro.r(
                   '\nfit <- lm({0}$Etot ~ poly({0}$Host_Scd, degree={1}, raw=TRUE))'
                   '\nsummary(fit, digits=10)'.format(systemname, fit_degree)
                   )
        with open('{}_Hself_{}_deg{}_polyfit.dat'.format(systemname, lipid, fit_degree), "w") as fitf:
            print(fit, file=fitf)
        plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x='Host_Scd', y='Etot')\
                + add_axis('x', 'Scd', breaks=[-0.3, 0.8, 0.2])\
                + add_axis('y', 'Etot', breaks=[-300, -255, 10])\
                + ggplot2.geom_point()\
                + ggplot2.geom_smooth(method='lm', se=False, fullrange=True, formula='y~poly(x, degree={})'.format(fit_degree))\
                + ggplot2.theme_light() + my_theme()
        print(plot)
        gr.pdf(filename)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    def plot_standard(self, systemname, lipid, filename, averaging=True, nchol=None):
        ''' Creates a standard EofScd plot of one system '''
        tavg = self.temperature_avg(systemname, systemname, lipid, averaging=averaging)
        print("TAVG", tavg)
        title = '{} self interaction'.format(lipid)
        if tavg is None:
            return
        if averaging:
            aes = ggplot2.aes_string(x='Host_Scd', y='Etot')
        else:
            aes = ggplot2.aes_string(x='Host_Scd', y='Etot', color='factor(temperature)')
        if nchol is not None:
            aes = ggplot2.aes_string(x='Host_Scd', y='Etot', color='factor(NChol)')
        xaxis = add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2])
        #yaxis = add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair])
        guides = ggplot2.guides(color=ggplot2.guide_legend("T / K"), linetype=ggplot2.guide_legend("T / K"), size=False)
        plot = ggplot2.ggplot(tavg)\
                + aes\
                + xaxis\
                + guides\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.7)\
                + ggplot2.ggtitle(title)\
                + ggplot2.scale_size("weight")\
                + ggplot2.theme_light() + my_theme()
                #+ ggplot2.geom_smooth(ggplot2.aes_string(linetype='factor(temperature)'), se=False)\
                                #+ ggplot2.scale_y_continuous(breaks=ro.FloatVector([i for i in range(-230, -150+1, 10)]), limits=ro.FloatVector([-230, -150]))\
        print(plot)
        #time.sleep(5)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    def bin_raw_data(self, outputtablename, datatablename, column='Host_Scd', breakl=-0.5, breakh=1.0, breakdx=0.05):
        ''' Discretize data by aggregating specified column '''
        print("Binning... {}".format(datatablename,))
        labelbreakh = breakh-breakdx
        mydf = ro.r(\
            '\n{0}$binnedScd <- cut({0}${1},'
                            'breaks=seq({2}, {3}, {4}),'
                            'labels=seq({2}, {5}, {4}))'
            '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
            '\nbckp_data <- na.omit({0})'
            '\n{0} <- na.omit({0})'
            '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=binnedScd]'
            '\n{6}$weight <- bckp_data[,.N,by=binnedScd][,N]'
            '\n{6}$weight <- {6}$weight/sum({6}$weight)'
            '\n{6}'\
            .format(datatablename, column,\
                    breakl, breakh, breakdx, labelbreakh,\
                    outputtablename)\
            )
        return mydf
    def temperature_avg(self, outputtable, systemname, lipid, column='Host_Scd', neib_spec=None, averaging=True):
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
        print("INPUT:", outputtable, systemname, lipid, )
        framenames = []
        for temp in self.temperatures[systemname]:
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            print(system)
            if self.interaction == 'complete':
                file_to_read = '{}/Eofscd_self{}.dat'.format(system, lipid)
            else:
                file_to_read = '{}/Eofscd_self{}{}.dat'.format(system, lipid, self.interaction)
            if os.path.exists(file_to_read):
                dat = read_datafile(system, file_to_read)
                print(dat)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
        #Averaging before everything's collected=============
            df = self.bin_raw_data(system, system, column=column)
            print("binned df:", df)
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dt = merge_data(outputtable, 'temperature', *framenames)
        print("MERGED DATA:", dt)
        if neib_spec is None:
            if averaging:
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binnedScd]'
                                '\ninsignificant <- 0.005*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                '\n{0}'
                                .format(outputtable))#column, ))
            else:
                df_final = dt
        else:
            raise NotImplementedError('Not yet implemented')
        return df_final

class EofScdplots():
    ''' Controls plotting of EofScd plots
        Needs raw_input files with following structure and header: 
            "Time Host Host_Scd Neib Neib_Scd Interaction DeltaScd AvgScd Etot Evdw Ecoul NChol Ntot"
    '''
    def __init__(self, interaction, systemnames=Tranges_whole.keys(), averaging=True, temperature='whole'):
        pairs = [('DPPC_DPPC', (-90, -40, 5)), ('DUPC_DUPC', (-90, -40, 5)), 
                 ('DPPC_DUPC', (-90, -40, 5)), ('DPPC_CHL1', (-50, -20, 5)),
                 ('DUPC_CHL1', (-50, -20, 5)), ('CHL1_CHL1', (-30, -10, 5))]
        self.enthalpyranges = {}
        self.systemnames = systemnames
        self.temperature_set = temperature
        self.temperatures = {}
        for syst in systemnames:
            if self.temperature_set == 'whole':
                self.temperatures[syst] = Tranges_whole[syst]
            elif self.temperature_set == 'high':
                self.temperatures[syst] = Tranges_high[syst]
        self.interaction = interaction
        for pair in pairs:
            self.enthalpyranges[pair[0]] = pair[1]

    def fit(self, systemname, pair, degree=6, nchol=None):
        '''
            Creates a fit for lipid-lipid interaction functions
        '''
    def fit_lipid_chol(self, systemname, pair, filename, mdpath, l_orderparam, fit_degree=3, neib_spec=None):
        '''
            Fits lipid chol interaction
        '''
        pairlipids = pair.split('_')
        lipidsequence = self.lipid_sequence_gro(mdpath+'/initial_coords/'+systemname+'.gro')
        for lip in lipidsequence:
            if lip in pairlipids:
                if lip == l_orderparam:
                    xcol = 'Host_Scd'
                else:
                    xcol = 'Neib_Scd'
                print("Taking:", xcol)
                break
        tavg = self.temperature_avg(systemname, systemname, pair, interaction_pair='lipid-chol',
                                    column=xcol, neib_spec=neib_spec, averaging=True)
        for Nc in range(5):
            fitdata = ro.r('my.formula={1}~poly({2},degree={3},raw=TRUE)'
                         '\nm <- lm(my.formula, data={0}[NChol=={4}])'
                         '\nsummary(m)'.format(systemname, 'Etot', 'AvgScd', fit_degree, Nc))
        #print("Fitdata: \n", fitdata, 3*"\n")
            with open('{}_eofscd_{}_deg{}_nchol{}_polyfit.dat'.format(systemname, pair, fit_degree, Nc), "w") as fitf:
                print(fitdata, file=fitf)
        plot = ggplot2.ggplot(tavg)\
                + ggplot2.aes_string(x=xcol, y='Etot', color='factor({})'.format(neib_spec))\
                + add_axis('x', 'scd', breaks=[-0.3, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=[-50, 20, 10], pair=pair)\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se=False, method='lm', fullrange=True, formula='y ~ poly(x, degree={}, raw=TRUE)'.format(fit_degree))\
                + ggplot2.ggtitle(pair+' interaction')\
                + ggplot2.scale_size("weight")\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
    def lipid_sequence_gro(self, gropath):
        ''' Determine in which the lipid species appear in the initial structure file '''
        lipidsequence = []
        with open(gropath,"r") as fgro:
            # get rid of first 2 lines:
            fgro.readline()
            fgro.readline() 
            for lines in fgro:
                lipid = lines[5:9]
                if lipid not in lipidsequence:
                    lipidsequence.append(lipid)
        return lipidsequence

    def plot_H_lipidchol(self, systemname, pair, filename, mdpath,  l_orderparam=None, averaging=True,
                         neib_spec=None,):
        ''' H(S_orderparam), H(S, nchol) Lipid-Chol interaction plots
            
            orderparam: <Lipid species> which order to take
            formula = ['poly', order]
        '''
        def lipidtypes_in_sys(sysname):
            lipids_in_sys = sysname.split('_')
            lipids_in_sys = [i[:4].upper() for i in lipids_in_sys]
            if len(lipids_in_sys) > 2:
                raise NotImplementedError('Not implemented for ternary systems')
            lipids_in_sys.remove('CHOL')
            return lipids_in_sys[0]
        titlepair = []
        pairlipids = pair.split('_')
        for lip in pairlipids:
            if lip == 'CHL1':
                titlepair.append('CHOL')
            else:
                titlepair.append(lip)
        if isinstance(systemname, list):
            pair_copy = pair
            tablelist = []
            for syst in systemname:
                if pair_copy == 'LIPID_CHL1':
                    print(pair, pairlipids)
                    pairlipids = [lipidtypes_in_sys(syst)] + ['CHL1']
                    pair = lipidtypes_in_sys(syst)+'_'+'CHL1'
                elif not self.lipid_in_sys(pair, syst):
                    continue
                lipidsequence = self.lipid_sequence_gro(mdpath+'/initial_coords/'+syst+'.gro')
                for lip in lipidsequence:
                    if lip in pairlipids:
                        if lip != 'CHL1':
                            xcol = 'Host_Scd'
                        else:
                            xcol = 'Neib_Scd'
                        break
                tavg = self.temperature_avg(syst, syst, pair, interaction_pair='lipid-chol',
                                            column=xcol, neib_spec=neib_spec, averaging=averaging)
                #print(syst, tavg)
                tablelist.append((syst, syst))
                print("TABLES", tablelist)
            tavg = merge_data('final_table', 'System', *tablelist)
            print("TAVG SYSTEMS", tavg)
            title = '-'.join(titlepair)
            if tavg is None:
                return
            if neib_spec is None:
                if averaging:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', color='System')
                else:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', shape='System', color='factor(temperature)')
                xaxis = add_axis('x', 'scd', breaks=[-0.2, 1.0, 0.2], pair=pair)
                yaxis = add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)
                guides = ggplot2.guides(color=ggplot2.guide_legend("System"), size=False)
            else:
                if averaging:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', shape='System', color='factor('+neib_spec+')',
                                            weight='weight')
                else:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', color='factor('+neib_spec+')',
                                            weight='weight', shape='factor(temperature)')
                xaxis = add_axis('x', 'scd', breaks=[-0.2, 1.0, 0.2], pair=pair)
                yaxis = ggplot2.scale_y_continuous(name= ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))'))
                guides = ggplot2.guides(color=ggplot2.guide_legend(neib_spec), size=False)
        else:
            lipidsequence = self.lipid_sequence_gro(mdpath+'/initial_coords/'+systemname+'.gro')
            for lip in lipidsequence:
                if lip in pairlipids:
                    if lip == l_orderparam:
                        xcol = 'Host_Scd'
                    else:
                        xcol = 'Neib_Scd'
                    print("Taking:", xcol)
                    break
            tavg = self.temperature_avg(systemname, systemname, pair, interaction_pair='lipid-chol',
                                        column=xcol, neib_spec=neib_spec, averaging=averaging)
            print("TAVG", tavg)
            title = '-'.join(titlepair)
            if tavg is None:
                return
            if neib_spec is None:
                if averaging:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot')
                else:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', color='factor(temperature)')
                xaxis = add_axis('x', 'scd', breaks=[-0.2, 1.0, 0.2], pair=pair)
                yaxis = add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)
                guides = ggplot2.guides(color=ggplot2.guide_legend("T / K"), size=False)
            else:
                if averaging:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', color='factor('+neib_spec+')', weight='weight')
                else:
                    aes = ggplot2.aes_string(x='binnedScd', y='Etot', color='factor('+neib_spec+')', weight='weight',
                                             shape='factor(temperature)')
                xaxis = add_axis('x', 'scd', breaks=[-0.2, 1.0, 0.2], pair=pair)
                yaxis = ggplot2.scale_y_continuous(name= ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))'))
                guides = ggplot2.guides(color=ggplot2.guide_legend(neib_spec), size=False)
        plot = ggplot2.ggplot(tavg)\
                + aes\
                + guides\
                + xaxis\
                + yaxis\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + ggplot2.geom_smooth(se=False)\
                + ggplot2.ggtitle(title+' interaction')\
                + ggplot2.scale_size("weight")\
                + ggplot2.theme_light() + my_theme()
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    def delta_scd(self, pair, systemnames, interaction_pair, fileoutputname, neib_spec=None):
        ''' '''
        tables = {}
        tablelist = []
        for sysname in systemnames:
            tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair, neib_spec=neib_spec)
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
    
    #===========================================================================
    # def parts_together(self, pair, sysname, filename, neib_spec=None):
    #     '''
    #     Plot all interactionpairs together
    #     '''
    #     if self.interaction == 'complete':
    #         return
    #     breaks_inter = {'head-tail':[-60, 10, 20],\
    #                     'head-tailhalfs':[-20, 0, 5],\
    #                     'carbons':[-2.5, 0, 0.1],\
    #                     }
    #     dt = self.temperature_avg('dataframe', sysname, pair, neib_spec=neib_spec, interaction_pair='all')
    #     if self.interaction == 'carbons' and neib_spec is None:
    #         dt = ro.r(\
    #                   'dataframe[,paste0("int",1:2) := tstrsplit(Interaction,"_")]'
    #                   '\ndataframe'
    #                   )
    #         print(dt)
    #         plot = ggplot2.ggplot(dt)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', weight='weight')\
    #             + ggplot2.geom_point()\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[0.2, 1.0, 0.4], limits=[-0.2, 1.0], pair=pair)\
    #             + add_axis('y', 'enthalpy', breaks=[-2.5, 0, 1], pair=pair)\
    #             + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #             + ggplot2.theme_light() + my_theme()\
    #             + ggplot2.facet_grid('int1 ~ int2')
    #         #save_plot(filename, plot, width=25, height=20)
    #         gr.pdf(file=filename,  width=10, height=7)
    #         print(plot)
    #         return gr.dev_off()
    #     elif self.interaction == 'carbons' and neib_spec is not None:
    #         dt = ro.r(\
    #                   'dataframe[,paste0("int",1:2) := tstrsplit(Interaction,"_")]'
    #                   '\ndataframe'
    #                   )
    #         print(dt)
    #         plot = ggplot2.ggplot(dt)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
    #             + ggplot2.geom_point(alpha=0.3)\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[0.2, 1.0, 0.4], limits=[-0.2, 1.0], pair=pair)\
    #             + add_axis('y', 'enthalpy', breaks=[-2.5, 0, 1], pair=pair)\
    #             + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
    #             + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #             + ggplot2.facet_grid('int1 ~ int2')\
    #             + ggplot2.theme_light() + my_theme()\
    #         #save_plot(filename, plot, width=25, height=20)
    #         gr.pdf(file=filename, width=10, height=7)
    #         print(plot)
    #         return gr.dev_off()
    #     elif self.interaction == 'head-tailhalfs' and neib_spec is not None:
    #         dt = ro.r(\
    #                   'dataframe <- dataframe[Interaction %in% c("t12_t12", "t22_t22")]'
    #                   '\ndataframe'
    #                   )
    #         print(dt)
    #         plot = ggplot2.ggplot(dt)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
    #             + ggplot2.geom_point(alpha=0.3)\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.4], pair=pair)\
    #             + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], limits=[-25, 0], pair=pair)\
    #             + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
    #             + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #             + ggplot2.facet_grid('Interaction ~ .')\
    #             + ggplot2.theme_light() + my_theme()\
    #         #save_plot(filename, plot, width=25, height=20)
    #         gr.pdf(file=filename, width=10, height=7)
    #         print(plot)
    #         return gr.dev_off()
    #     elif neib_spec is None:
    #         plot = ggplot2.ggplot(dt)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', color='Interaction', weight='weight')\
    #             + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
    #             + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
    #             + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #             + ggplot2.theme_light() + my_theme()
    #     else:
    #         plot = ggplot2.ggplot(dt)\
    #                 + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(NChol)', weight='weight')\
    #                 + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.6)\
    #                 + ggplot2.geom_smooth(ggplot2.aes_string(), size=1., se='False')\
    #                 + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
    #                 + add_axis('y', 'enthalpy', breaks=breaks_inter[self.interaction], pair=pair)\
    #                 + ggplot2.theme_light() + my_theme()\
    #                 + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #                 + ggplot2.guides(color=ggplot2.guide_legend("NChol"))\
    #                 + ggplot2.facet_grid('Interaction ~ .')
    #     gr.pdf(file=filename)
    #     print(plot)
    #     gr.dev_off()
    #     ro.r('rm(list=ls(all=TRUE))')
    #===========================================================================

    @staticmethod
    def lipid_in_sys(lipidpair, sysname):
        lipids_in_sys = sysname.split('_')
        #print(lipids_in_sys)
        lipids_in_sys = [i[:4].upper() for i in lipids_in_sys]
        #print(lipids_in_sys)
        if 'CHOL' in lipids_in_sys:
            lipids_in_sys.append('CHL1')
            #print(lipids_in_sys)
        for lipid in lipidpair.split('_'):
            #print(lipid)
            if lipid not in lipids_in_sys:
                return 0
        return 1

    def all_systems(self, pair, interaction_pair, filename, systems=None, neib_spec=None):
        '''
        Plot all systems together 
        '''
        tables = {}
        tablelist = []
        if systems is None:
            for sysname in self.systemnames:
            #for sysname in ['dppc_dupc_chol', 'dppc_chol']:
                if not self.lipid_in_sys(pair, sysname):
                    continue
                tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair, neib_spec=neib_spec)
                tablelist.append((sysname, sysname))
        else:
            for sysname in systems:
            #for sysname in ['dppc_dupc_chol', 'dppc_chol']:
                if not self.lipid_in_sys(pair, sysname):
                    print("Lipid not in Sys:", pair, sysname)
                    continue
                tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair, neib_spec=neib_spec)
                tablelist.append((sysname, sysname))
        #'\n{0}fin$weight <- {0}fin$weight/sum({0}fin$weight)'
        print(tablelist)
        if not tablelist:
            return
        final_table = merge_data('final_table', 'System', *tablelist)
        print(final_table)
        titlepair = []
        for lip in pair.split('_'):
            if lip == 'CHL1':
                titlepair.append('CHOL')
            else:
                titlepair.append(lip)
        title = '-'.join(titlepair)
        if neib_spec is None:
            plot = ggplot2.ggplot(final_table)\
                + ggplot2.aes_string(x='AvgScd', y='Etot', color='System', weight='weight')\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
                + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)\
                + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                + ggplot2.guides(size=False)\
                + ggplot2.scale_size("weight")\
                + ggplot2.ggtitle(title+' interaction')\
                + ggplot2.theme_light() + my_theme()
        elif neib_spec == 'Ntot':
            print("PLOTTING NTOT")
            plot = ggplot2.ggplot(final_table)\
                    + ggplot2.aes_string(x='AvgScd', y='Etot', shape='factor(System)', color='factor(bin'+neib_spec+')', weight='weight')\
                    + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.7)\
                    + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)\
                    + ggplot2.scale_y_continuous(name= ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))'))\
                    + ggplot2.guides(color=ggplot2.guide_legend(neib_spec), shape=ggplot2.guide_legend("System"), size=False)\
                    + ggplot2.scale_size("weight")\
                    + ggplot2.ggtitle(title+' interaction')\
                    + ggplot2.theme_light() + my_theme() + ggplot2.theme(**{'legend.key.width':ro.r.unit(1,"cm")})
                    #===========================================================
                    # + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
                    #==========================================================
        else:
            #final_table = ro.r("final_table <- final_table[NChol==3]")
            #+ ggplot2.scale_y_continuous(name= ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))'))\
            plot = ggplot2.ggplot(final_table)\
                    + ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(System)', shape='factor('+neib_spec+')', weight='weight')\
                    + add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)\
                    + add_axis('y', 'enthalpy', limits=[-90, -40], breaks=[-90, -40, 5])\
                    + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.7)\
                    + ggplot2.guides(shape=ggplot2.guide_legend(neib_spec), color=ggplot2.guide_legend("System"), linetype=ggplot2.guide_legend("System"),size=False)\
                    + ggplot2.scale_size("weight")\
                    + ggplot2.ggtitle(title+' interaction')\
                    + ggplot2.theme_light() + my_theme() + ggplot2.theme(**{'legend.key.width':ro.r.unit(1,"cm")})
                    #===========================================================
                    # + ggplot2.geom_smooth(ggplot2.aes_string(linetype="factor(System)"), se=False, size=1)\
                    # + ggplot2.scale_linetype_manual(values=ro.r('c("solid", "dashed", "dotted", "dotdash")'))\
                    #===========================================================
        #print(plot)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')

    #===========================================================================
    # def energy_components(self, pair, systemname, interaction_pair, filename):
    #     '''
    #     Plot Vdw and Coulomb separately
    #     '''
    #     #=======================================================================
    #     # tables = {}
    #     # tablelist = []
    #     # for sysname in self.systemnames:
    #     #     tables[sysname] = self.temperature_avg(sysname, sysname, pair, interaction_pair)
    #     #     tablelist.append((sysname, sysname))
    #     # final_table = merge_data('final_table', 'System', *tablelist)
    #     #=======================================================================
    #     final_table = self.temperature_avg(systemname, systemname, pair, interaction_pair)
    #     plot = ggplot2.ggplot(final_table)\
    #             + ggplot2.aes_string(size='weight')\
    #             + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Evdw'), color='red')\
    #             + ggplot2.geom_point(ggplot2.aes_string(x='AvgScd', y='Ecoul'), color='blue')\
    #             + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Evdw'), color='red')\
    #             + ggplot2.geom_smooth(ggplot2.aes_string(x='AvgScd', y='Ecoul'),color='blue')\
    #             + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2], pair=pair)\
    #             + ggplot2.ggtitle('-'.join(pair.split("_"))+' interaction')\
    #             + ggplot2.theme_light() + my_theme()
    #             #+ add_axis('y', 'enthalpy', breaks=[-50, 0, 5], pair=pair)\
    #     gr.pdf(file=filename)
    #     print(plot)
    #     gr.dev_off()
    #     ro.r('rm(list=ls(all=TRUE))')
    #===========================================================================

    #===========================================================================
    # def pairs_together(self, datasources, filename, interaction_pair='w_W', neib_spec=None): 
    #     tables = {}
    #     tablelist = []
    #     for sysname, pair in datasources:
    #         tables[sysname] = self.temperature_avg(pair, sysname, pair, interaction_pair, neib_spec=neib_spec)
    #         tablelist.append((pair, pair))
    #         #'\n{0}fin$weight <- {0}fin$weight/sum({0}fin$weight)'
    #     final_table = merge_data('final_table', 'Pair', *tablelist)
    #     print(final_table)
    #     plot = ggplot2.ggplot(final_table)\
    #             + ggplot2.aes_string(x='AvgScd', y='Etot', color='Pair', weight='weight')\
    #             + ggplot2.geom_point(ggplot2.aes_string(size='weight'))\
    #             + ggplot2.geom_smooth(se='False')\
    #             + add_axis('x', 'avgscd', breaks=[-0.4, 1.0, 0.2])\
    #             + add_axis('y', 'enthalpy', breaks=[-50, -10, 10])\
    #             + ggplot2.scale_size("weight")\
    #             + ggplot2.theme_light() + my_theme()
    #             #                + ggplot2.geom_smooth(se='False')\
    #     gr.pdf(file=filename, width=10, height=7)
    #     print(plot)
    #     gr.dev_off()
    #===========================================================================

    def standard_plot(self, systemname, pair, interaction_pair, filename, averaging=True, neib_spec=None):
        ''' Creates a standard EofScd plot of one system '''
        tavg = self.temperature_avg(systemname, systemname, pair, interaction_pair, neib_spec=neib_spec, averaging=averaging)
        print("TAVG", tavg)
        titlepair = []
        for lip in pair.split('_'):
            if lip == 'CHL1':
                titlepair.append('CHOL')
            else:
                titlepair.append(lip)
        title = '-'.join(titlepair)
        if tavg is None:
            return
        fitdata = ro.r('my.formula={0}${1}~poly({0}${2},4,raw=TRUE)'
           '\nm <- lm(my.formula, data={0})'
           '\nsummary(m)'.format(systemname, 'Etot', 'AvgScd', ))
        print(fitdata)
        if neib_spec is None:
            if averaging:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot')
            else:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(temperature)')
            xaxis = add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)
            yaxis = add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)
            guides = ggplot2.guides(color=ggplot2.guide_legend("T / K"), linetype=ggplot2.guide_legend("T / K"), size=False)
        elif neib_spec == 'Ntot':
            if averaging:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(bin'+neib_spec+')', weight='weight')
            else:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot', color='factor(bin'+neib_spec+')', weight='weight',
                                        shape='factor(temperature)')
            xaxis = add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)
            yaxis = ggplot2.scale_y_continuous(name= ro.r('expression(paste("H","(S"[CD],") / kJ/mol"))'))
            guides = ggplot2.guides(color=ggplot2.guide_legend(neib_spec), size=False)
            #===========================================================
            # + add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)\
            #===========================================================
        else:
            if averaging:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot', color='factor('+neib_spec+')', weight='weight')
            else:
                aes = ggplot2.aes_string(x='AvgScd', y='Etot', color='factor('+neib_spec+')', weight='weight',
                                         shape='factor(temperature)')
            xaxis = add_axis('x', 'avgscd', breaks=[-0.2, 1.0, 0.2], pair=pair)
            yaxis = add_axis('y', 'enthalpy', breaks=self.enthalpyranges[pair], pair=pair)
            guides = ggplot2.guides(color=ggplot2.guide_legend(neib_spec), shape=ggplot2.guide_legend("T / K"), 
                                    linetype=ggplot2.guide_legend("T / K"), size=False)
        plot = ggplot2.ggplot(tavg)\
                + aes\
                + guides\
                + xaxis\
                + yaxis\
                + ggplot2.geom_point(ggplot2.aes_string(size='weight'), alpha=0.7)\
                + ggplot2.ggtitle(title+' interaction')\
                + ggplot2.scale_size("weight")\
                + ggplot2.theme_light() + my_theme()
                #+ ggplot2.geom_smooth(ggplot2.aes_string(linetype='factor(temperature)'), se=False)\
        #print(plot)
        #time.sleep(5)
        gr.pdf(file=filename, width=10, height=7)
        print(plot)
        gr.dev_off()
        ro.r('rm(list=ls(all=TRUE))')
    def bin_raw_data(self, outputtablename, datatablename, interaction_pair, column='AvgScd', breakl=-0.5, breakh=1.0, breakdx=0.05, neib_spec=None):
        ''' Discretize data by aggregating specified column '''
        print("Binning... {} : {}".format(datatablename, interaction_pair))
        labelbreakh = breakh-breakdx
        print("interaction pair is:", interaction_pair)
        if interaction_pair == 'complete':
            if neib_spec == 'Ntot':
                print("NTOT")
                mydf = ro.r(\
                    '\n{0}$binnedScd <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))'
                    '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
                    '\n{0}$binNtot <- cut({0}$Ntot,'
                                    'breaks=seq(0, 14, 2),)'
                                    #'labels=seq(1, 1, 2))'
                    '\nbckp_data <- na.omit({0})'
                    '\n{0} <- na.omit({0})'
                    '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction, binNtot)]'
                    '\n{6}$weight <- bckp_data[,.N,by=list(binnedScd, Interaction, binNtot)][,N]'
                    '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                    '\n{6}'\
                    .format(datatablename, column,\
                            breakl, breakh, breakdx, labelbreakh,\
                            outputtablename)\
                    )
            elif neib_spec is not None:
                mydf = ro.r(\
                    '\n{0}$binnedScd <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))'
                    '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
                    '\nbckp_data <- {0}'
                    #'\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction, Ntot, {7})]'
                    '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction, {7})]'
                    '\n{6}$weight <- bckp_data[,.N,by=list(binnedScd, Interaction, {7})][,N]'
                    '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                    '\n{6}'\
                    .format(datatablename, column,\
                            breakl, breakh, breakdx, labelbreakh,\
                            outputtablename, neib_spec)\
                    )
            elif neib_spec is None:
                mydf = ro.r(\
                    '\n{0}$binnedScd <- cut({0}${1},'
                                    'breaks=seq({2}, {3}, {4}),'
                                    'labels=seq({2}, {5}, {4}))'
                    '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
                    '\nbckp_data <- {0}'
                    '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction)]'
                    '\n{6}$weight <- bckp_data[,.N,by=list(binnedScd, Interaction)][,N]'
                    '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                    '\n{6}'\
                    .format(datatablename, column,\
                            breakl, breakh, breakdx, labelbreakh,\
                            outputtablename)\
                    )

        elif interaction_pair == 'lipid-chol':
            if neib_spec is not None:
                mydf = ro.r(\
                            '\n{0}$binnedScd <- cut({0}${1},'
                                            'breaks=seq({2}, {3}, {4}),'
                                            'labels=seq({2}, {5}, {4}))'
                            '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
                            '\nbckp_data <- {0}'
                            '\n{7} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction, {6})]'
                            '\n{7}$weight <- bckp_data[,.N,by=list(binnedScd, Interaction, {6})][,N]'
                            '\n{7}$weight <- {7}$weight/sum({7}$weight)'
                            '\n{7}'\
                            .format(datatablename, column,
                                    breakl, breakh, breakdx, labelbreakh,
                                    neib_spec, outputtablename),
                            )
            elif neib_spec is None:
                mydf = ro.r(\
                            '\n{0}$binnedScd <- cut({0}${1},'
                                            'breaks=seq({2}, {3}, {4}),'
                                            'labels=seq({2}, {5}, {4}))'
                            '\n{0}$binnedScd <- as.numeric(levels({0}$binnedScd))[{0}$binnedScd]'
                            '\nbckp_data <- {0}'
                            '\n{6} <- {0}[,lapply(.SD, function(x) if (!is.numeric(x)) x else mean(x)), by=list(binnedScd, Interaction)]'
                            '\n{6}$weight <- bckp_data[,.N,by=list(binnedScd, Interaction)][,N]'
                            '\n{6}$weight <- {6}$weight/sum({6}$weight)'
                            '\n{6}'\
                            .format(datatablename, column, breakl, breakh, breakdx,
                                    labelbreakh, outputtablename, ),
                            )
        else:
            raise NotImplementedError("Sorry, not implemented anymore...")
        #print("Thats my df:", mydf)
        print("Finished binning")
        return mydf



    def temperature_avg(self, outputtable, systemname, pair, interaction_pair, column='AvgScd', neib_spec=None, averaging=True):
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
        for temp in self.temperatures[systemname]:
            print("At temperature:", temp)
            system = '{}_{}'.format(systemname, temp)
            #print(system)
            if self.interaction == 'complete':
                file_to_read = '{}/Eofscd{}.dat'.format(system, pair)
            elif self.interaction == 'lipid-chol':
                file_to_read = '{}/Eofscd{}.dat'.format(system, pair)
            else:
                file_to_read = '{}/Eofscd{}{}.dat'.format(system, pair, self.interaction)
            if os.path.exists(file_to_read):
                read_datafile(system, file_to_read)
            else:
                print("File {} does not exist".format(file_to_read))
                return None
        #Averaging before everything's collected=============
            df = self.bin_raw_data(system, system, interaction_pair, column=column, neib_spec=neib_spec)
            #print("binned df:", df)
            framenames.append(system)
        framenames = [(i, i[-3:]) for i in framenames]
        print("Framenames:", framenames)
        dt = merge_data(outputtable, 'temperature', *framenames)
        #print("MERGED DATA:", dt)
        if not averaging:
            #dt = ro.r('{0} <- {0}[Lipid_type=="{1}"]\n'.format(outputtable))
            if neib_spec is not None:
                df_final = ro.r('\n{0}Tlow <- {0}[temperature<325][,temperature:=NULL]'
                                '\n{0}Thigh <- {0}[temperature>325][,temperature:=NULL]'
                                '\n{0}compl <- {0}[,temperature:=NULL]'
                                '\n{0}Tlow <- {0}Tlow[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                '\n{0}Thigh <- {0}Thigh[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                '\n{0}compl <- {0}compl[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                '\n{0}Tlow$temperature <- "290-320 K"'
                                '\n{0}Thigh$temperature <- "330-350 K"'
                                '\n{0}compl$temperature <- "complete"'
                                '\ninsignificantlow <- 0.015*max({0}Tlow$weight)'
                                '\ninsignificanthigh <- 0.015*max({0}Thigh$weight)'
                                '\ninsignificantcompl <- 0.015*max({0}compl$weight)'
                                '\n{0}Tlow <- {0}Tlow[weight>insignificantlow]'
                                '\n{0}Thigh <- {0}Thigh[weight>insignificanthigh]'
                                '\n{0}compl <- {0}compl[weight>insignificantcompl]'
                                '\n{0} <- rbind({0}Tlow, {0}Thigh, {0}compl)'
                                '\n{0}'.format(outputtable, neib_spec))
            else:
                df_final = ro.r('\n{0}Tlow <- {0}[temperature<325][,temperature:=NULL]'
                                '\n{0}Thigh <- {0}[temperature>325][,temperature:=NULL]'
                                '\n{0}compl <- {0}[,temperature:=NULL]'
                                '\n{0}Tlow <- {0}Tlow[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction)]'
                                '\n{0}Thigh <- {0}Thigh[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction)]'
                                '\n{0}compl <- {0}compl[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction)]'
                                '\n{0}Tlow$temperature <- "290-320 K"'
                                '\n{0}Thigh$temperature <- "330-350 K"'
                                '\n{0}compl$temperature <- "complete"'
                                '\ninsignificantlow <- 0.015*max({0}Tlow$weight)'
                                '\ninsignificanthigh <- 0.015*max({0}Thigh$weight)'
                                '\ninsignificantcompl <- 0.015*max({0}compl$weight)'
                                '\n{0}Tlow <- {0}Tlow[weight>insignificantlow]'
                                '\n{0}Thigh <- {0}Thigh[weight>insignificanthigh]'
                                '\n{0}compl <- {0}compl[weight>insignificantcompl]'
                                '\n{0} <- rbind({0}Tlow, {0}Thigh, {0}compl)'
                                '\n{0}'.format(outputtable))
            return df_final
        if interaction_pair == 'complete':
            if neib_spec == 'Ntot':
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, binNtot)]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                '\n{0}'
                                    #'{0} <- {0}[,Lipid_type:=NULL][weight>insignificant]\n'
                             #'\n{1}$NChol <- factor({1}$NChol)'
                             .format(outputtable))#, column))
                print("FINAL:", df_final)
            elif neib_spec is not None:
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                #'\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                #'\n{0} <- {0}[Ntot<8][Ntot>4]'
                                '\n{0} <- {0}[NChol<5]'
                                #'\n{0} <- {0}[Ntot==6]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                '\n{0}'
                                .format(outputtable, neib_spec))#column, ))
                #print("this is what remains", df_final)
            elif neib_spec is None:
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction)]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                .format(outputtable))#, column))
        elif interaction_pair == 'lipid-chol':
            if neib_spec is None:
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction)]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                .format(outputtable))#, column))
            elif neib_spec is not None:
                df_final = ro.r('{0} <- subset({0},select=-c(temperature))'
                                '\n{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=list(binnedScd, Interaction, {1})]'
                                '\ninsignificant <- 0.01*max({0}$weight)'
                                '\n{0} <- {0}[weight>insignificant]'
                                '\n{0} <- {0}[NChol<7]'
                                '\n{0}'
                                .format(outputtable, neib_spec))# column))
        else:
            raise NotImplementedError("Not working anymore...")
            #===================================================================
            # if neib_spec is None:
            #     df_final = ro.r('{0} <- subset({0},select=-c(temperature))\n'
            #                     '{0} <- {0}[,lapply(.SD,weighted.mean,w=weight),by=binAvgScd]'\
            #                     '\ninsignificant <- 0.01*max({0}$weight)'
            #                     '\n{0} <- {0}[weight>insignificant]'
            #                  .format(outputtable))
            #===================================================================
        #print("FINALFRAME:", df_final)
        return df_final




