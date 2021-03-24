# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:05:00 2020
Revised on Thur Mar 18 16:04:00 2021

@author: Starlitnightly

New Version 1.2.3
"""
import itertools
import numpy as np
import pandas as pd
from upsetplot import from_memberships
from upsetplot import plot


def FindERG(data, depth=2, sort_num=20, verbose=False):
    '''
    Find out endogenous reference gene

    Parameters
    ----------
    data:pandas.DataFrmae
        DataFrame of data points with each entry in the form:['gene_id','sample1',...]
    depth:int
        Accuracy of endogenous reference gene,must be larger that 2
        The larger the number, the fewer genes are screened out,Accuracy improvement
    sort_num:int
    	The size of the  peendogenous reference gener filter
    	When the sample is large, it is recommended to increase the value
    verbose: bool
        Make the function noisy, writing times and results.
    Returns
    -------
    result:list
        a list of endogenous reference gene
    '''
    lp=[]
    if verbose:
        import time,datetime
        start = time.time()
    if depth==1:
        print('the depth must larger than 2')
        return
    if len(data.columns)<=2:
        print('the number of samples must larger than 2')
        return
    if depth>(len(data.columns)-1):
        print('depth larger than samples')
        return
    count=0
    result=[]#result
    bucket_size = 1000
    for i in itertools.combinations(data.columns[0:depth], 2):
        start = time.time()
        count=count+1
        test=data.replace(0,np.nan).dropna()

        last_std=pd.DataFrame()
        for k in range(0 ,len(data), bucket_size):
            test1=test[i[0]].iloc[k:k + bucket_size]
            test2=test[i[1]].iloc[k:k + bucket_size]
            data_len=len(test1.values)
            table1=np.array(test1.values.tolist()*data_len).reshape(data_len,data_len)
            table2=pd.DataFrame(table1.T/table1)
            table2.index=test1.index

            table4=np.array(test2.values.tolist()*data_len).reshape(data_len,data_len)
            table5=pd.DataFrame(table4.T/table4)
            table5.index=test1.index

            table6=(table2-table5).std()
            table6.index=test1.index
            l_std=table6.sort_values()[0:sort_num]
            if(k==0):
                last_std=l_std
            else:
                last_std=pd.concat([last_std,l_std])

        last_std=last_std.sort_values()[0:sort_num]

        testlist=list(last_std.index)
        #print(testlist)
        lp.append(testlist)
        #print(lllll)
        if(count==1):
            result=testlist
        if(count>1):
            result=list(set(testlist).intersection(set(result))) #Venn
    example = from_memberships(lp,data=range(len(lp)))
    if verbose:
        end = time.time()
        print("calculate time:%.2fs"%(end-start))
        print(result)
    if depth>2:
        plot(example)
    return result


def normalizationdata( data, ERGname):
    '''
    Get a DataFrame disposed by standardization

    Parameters
    ----------
    data:pandas.DataFrmae
        DataFrame of data points with each entry in the form:['gene_id','sample1',...]
    ERGname:str
        the gene use to normalization

    Returns
    -------
    result:pandas.DataFrame
        A DataFrame disposed by standardization
    '''
    l=data[data.index.isin([ERGname])]
    da=data.iloc[:,0:len(data.columns)]/np.asarray(l.iloc[0,0:len(data.columns)])
    da=da*l.iloc[0,0]
    return da
