"""Inspired by https://stackoverflow.com/questions/8951787/defining-python-decorators-for-a-complete-module
(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""

import functools
import datetime,os

from extract_ev_for_clustering import parse_results as extract_ev_for_clustering
from  extract_subfastas import main as extract_subfastas, extract_seq_from_hhm as extract_seq_from_hhm
from get_representatives_from_clusters import main as get_representatives_from_clusters
from make_M_N import main as make_M_N
from my_little_merger_error_cor import main as my_little_merger_error_cor
from  my_little_replacer import main as my_little_replacer

class Logger(object):
    def __new__(cls,*args,**kwargs):
        if not hasattr(cls,"_my_little_logger"):
            cls._my_little_logger = super(Logger,cls).__new__(cls,*args,**kwargs)
        return cls._my_little_logger
    def __init__(self, outputname=''):
        if not outputname:
            outputname = "{}/debug_{}.log".format(os.getcwd(),datetime.datetime.now().strftime("%d-%m-%y_%H:%M%S"))
        self.outputfn = outputname
        self.lastlog = None
    def log(self,fn,fnm,what):
        out = ''
        now = datetime.datetime.now()
        timestamp = now.strftime("%d-%m-%y_%H:%M%S")
        if self.lastlog is not None:
            out += "ELAPSED for {} ({}): {}\n".format(fn,fnm, now-self.lastlog)
        self.lastlog = now
        out += "{ts}:\t{fn} ({fnm})\t{w}\n".format(fn=fn,fnm=fnm,ts=timestamp,w=what)
        self.write(out)
    def write(self,data):
        with open(self.outputfn,"a",0) as out:
            out.write(data)

def tracer(func):
    logger = Logger()
    @functools.wraps(func)
    def logwrap(*args,**kwargs):
        logger.log(fn=func.__name__,fnm=func.__module__,what="calling with {ar} {kw}\n".format(ar=args,kw=kwargs))
        now = datetime.datetime.now()
        result = func(*args,**kwargs)
        logger.log(fn=func.__name__,fnm=func.__module__,what="returned: {res} after : {elaps}".format(res=result,elaps=datetime.datetime.now()-now))
        return result
    return logwrap

funcs_to_wrap = ['extract_ev_for_clustering', 'extract_subfastas', 'extract_seq_from_hhm', 'get_representatives_from_clusters',
                 'make_M_N', 'my_little_merger_error_cor', 'my_little_replacer']

for name in funcs_to_wrap:
    globals()[name] = tracer(globals()[name])
if __name__ == "__main__":
    print vars(globals()['make_M_N'])
    #print globals()