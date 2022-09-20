"""

.. module:: dataframetools
   :synopsis: Some convenient tools for handling pandas dataframes in the context of htpolynet coordinates
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import pandas as pd
import logging

logger=logging.getLogger(__name__)

def get_row_attribute(df:pd.DataFrame,name,attributes):
    """get_row_attribute returns a scalar value of attribute "name" in row
        expected to be uniquely defined by attributes dict

    :param df: dataframe to search
    :type df: pandas.DataFrame
    :param name: name of attribute whose value you want
    :type name: str
    :param attributes: dictionary of attribute:value pairs that defines target set or row
    :type attributes: dict
    :return: value of attribute name
    :rtype: scalar
    """
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find row with attributes {attributes}'
    if type(name)==list:
        name_in_df=all([n in df for n in name])
    else:
        name_in_df= name in df
    assert name_in_df,f'Attribute(s) {name} not found'
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
    return df[list(l)][name].values[0]

def get_row_as_string(df:pd.DataFrame,attributes):
    """get_row_as_string returns a scalar value of attribute "name" in row
        expected to be uniquely defined by attributes dict

    :param df: a pandas dataframe
    :type df: pd.DataFrame
    :param attributes: dictionary of column names (keys) and values that specify set of rows to be returned
    :type attributes: dict(str,obj)
    :return: selected dataframe converted to a string
    :rtype: str
    """
    ga={k:v for k,v in attributes.items() if k in df}
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
    return df[list(l)].to_string()

def get_rows_w_attribute(df:pd.DataFrame,name,attributes:dict):
    ''' returns a series of values of attribute "name" from
        all rows matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find any rows with attributes {attributes}'
    if type(name)==list:
        name_in_df=all([n in df for n in name])
    else:
        name_in_df= name in df
    assert name_in_df,f'Attribute(s) {name} not found'
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
    return df[list(l)][name].values

def set_row_attribute(df:pd.DataFrame,name,value,attributes):
    ''' set value of attribute name to value in all rows matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    exla={k:v for k,v in attributes.items() if not k in df}
    if len(exla)>0:
        logger.warning(f'Caller attempts to use unrecognized attributes to refer to row: {exla}')
    if name in df and len(ga)>0:
        c=[df[k] for k in ga]
        V=list(ga.values())
        l=[True]*df.shape[0]
        for i in range(len(c)):
            l = (l) & (c[i]==V[i])
        cidx=[c==name for c in df.columns]
        df.loc[list(l),cidx]=value

def set_rows_attributes_from_dict(df:pd.DataFrame,valdict,attributes):
    ''' set values of attributes in valdict dict of all rows
        matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    exla={k:v for k,v in attributes.items() if not k in df}
    if len(exla)>0:
        logger.warning(f'using unknown attributes to refer to atom: {exla}')
    if all([x in df for x in valdict]) and len(ga)>0:
        c=[df[k] for k in ga]
        V=list(ga.values())
        l=[True]*df.shape[0]
        for i in range(len(c)):
            l = (l) & (c[i]==V[i])
        for k,v in valdict.items():
            cidx=[c==k for c in df.columns]
            df.loc[list(l),cidx]=v 
