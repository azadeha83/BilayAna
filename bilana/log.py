'''
    BilAnas logger
    Creates a global logging object that should be imported by all other
'''
import logging
LOGGER = logging.getLogger("BilAna")
LOGGER.setLevel(logging.INFO)

# Creating Streamhandler -- This handler will write to console
SH = logging.StreamHandler()
SH.setLevel(logging.DEBUG)
SH.setFormatter(logging.Formatter('%(asctime)s %(funcName)s - %(levelname)s - %(message)s'))
LOGGER.addHandler(SH)


def create_filehandler(logname, loginstance):
    ''' Creating Filehandler -- Writing to Bilana.log '''
    FH = logging.FileHandler(logname)
    FH.setLevel(logging.INFO)
    FH.setFormatter(logging.Formatter('%(asctime)s in %(funcName)s:  %(message)s\n'))
    loginstance.addHandler(FH)
    return loginstance

def set_verbosity(slevel):
    ''' Sets level of of LOGGER instance
        slevels:
            0 Warning
            1 Info
            2 Debug
    '''
    if isinstance(slevel, int):
        if slevel is None:
            slevel = 0
        elif slevel > 2:
            slevel = 2
        slevel = 30 -slevel*10
        LOGGER.setLevel(slevel)
    elif isinstance(slevel, str):
        LOGGER.setLevel(slevel)
