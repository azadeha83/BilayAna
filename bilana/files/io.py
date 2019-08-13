
'''
    Functions that read files that were created with Bilana and translates data to dict.
    NOTE: This is probably inefficient and could be done using pandas or the like
'''
def read_energyinput(energyfile):
    timetoenergy = {}
    with open(energyfile, "r") as efile:
        efile.readline()
        for line in efile:
            #cols = [x.strip() for x in line.split(' ')]
            cols = line.split()
            if int(cols[1]) > int(cols[2]):
                continue
            time = float(cols[0])
            respair = (int(cols[1]), int(cols[2]))
            #print(time, respair, end='\r')
            Etot = float(cols[6])
            VDW = float(cols[5])
            COUL = float(cols[4])
            inttype = cols[3]
            timetoenergy.update({(time, respair, inttype):(Etot, VDW, COUL)})
    return timetoenergy, time

def read_scdinput(scdfile):
    timetoscd = {}
    with open(scdfile, "r") as sfile:
        header = sfile.readline().split()
        for line in sfile:
            cols = line.split()
            time = float(cols[header.index("Time")])
            res = int(cols[header.index("Residue")])
            scd = float(cols[header.index("Scd")])
            timetoscd.update({(time, res):scd})
    return timetoscd, time

def read_neighborinput(neighborfile):
    neighbors_of_host = {}
    hosts_without_neib = []
    with open(neighborfile, "r") as nfile:
        nfile.readline()
        for line in nfile:
            cols = line.split()
            time = float(cols[1])
            host = int(cols[0])
            if int(cols[2]) == 0:
                #print(host, "has no neighbors at time", time)
                hosts_without_neib += [(time, host)]
                continue
            neighbors = cols[3]
            neighbors_of_host.update({(time, host):neighbors})
    return neighbors_of_host, hosts_without_neib
