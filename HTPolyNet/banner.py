"""

.. module:: banner
   :synopsis: Defines the banner method for printing the banner to a logging channel
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from HTPolyNet.stringthings import my_logger
from HTPolyNet.__init__ import HTPOLYNET_VERSION

banner_message="""
    HTPolyNet {:s}
    https://abramsgroup.github.io/HTPolyNet/

    Ming Huang
    mh3429@dragons.drexel.edu

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants W911NF-17-2-0227 
    and W911NF-12-R-0011 from the US Army Research Lab

    Please cite the HTPolyNet paper:
    
    Ming Huang and Cameron F. Abrams, HTPolyNet: A general 
    system generator for all-atom molecular simulations of 
    amorphous crosslinked polymers, SoftwareX, vol. 21, 
    pp. 101303, 2023 (doi:10.1016/j.softx.2022.101303)
    """.format(HTPOLYNET_VERSION)
def banner(logf):
    my_logger(banner_message,logf,fill=' ',just='<')