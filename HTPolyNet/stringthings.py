import pandas as pd

def my_logger(msg,logf,width=67,fill='*',sep=', ',just='^'):
    fmt=r'{'+r':'+fill+just+f'{width}'+r'}'
    ll=' ' if just in ['^','>'] else ''
    rr=' ' if just in ['^','<'] else ''
    if type(msg)==list:
        rr=' ' if ' ' not in sep else ''
        lnlst=[]
        for tok in msg:
            ll_sofar=sum([len(x) for x in lnlst])
            test_ll=ll_sofar+len(tok)+len(sep)*(len(lnlst)+1)
            if test_ll>width-2*(1+len(sep)):
                outstr=ll+sep.join(lnlst)+sep+rr
                logf(fmt.format(outstr))
                lnlst=[tok]
            else:
                lnlst.append(tok)
        outstr=ll+sep.join(lnlst)+' '
        logf(fmt.format(outstr))
    elif type(msg)==pd.DataFrame:
        for ln in msg.to_string().split('\n'):
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
    else:
        lns=msg.split('\n')
        for ln in lns:
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
        
