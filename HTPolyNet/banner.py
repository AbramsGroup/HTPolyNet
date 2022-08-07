from HTPolyNet.stringthings import my_logger
def banner(logf):
    banner_message="""
    HTPolyNet
    https://abramsgroup.github.io/HTPolyNet/

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants W911NF-17-2-0227, W911NF-12-R-0011 
    """
    my_logger(banner_message,logf,pad=' ',just='<')