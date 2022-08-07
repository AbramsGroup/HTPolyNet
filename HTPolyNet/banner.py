from HTPolyNet.stringthings import my_logger
def banner(logf):
    banner_message="""
    HTPolyNet
    https://abramsgroup.github.io/HTPolyNet/

    Ming Huang
    mh3429@dragons.drexel.edu

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants W911NF-17-2-0227 
    and W911NF-12-R-0011 from the US Army Research Lab
    """
    my_logger(banner_message,logf,fill=' ',just='<')