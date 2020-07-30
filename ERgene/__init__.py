# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:05:00 2020

@author: Starlitnightly

New Version 1.1.1
"""

import itertools
import numpy as np
import pandas as pd
from upsetplot import from_memberships
from upsetplot import plot



def FindERG( data, depth=2):
    '''
    Find out endogenous reference gene
    
    Parameters
    ----------
    data:pandas.DataFrmae
        DataFrame of data points with each entry in the form:['gene_id','sample1',...]
    depth:int
        Accuracy of endogenous reference gene,must be larger that 2
        The larger the number, the fewer genes are screened out,Accuracy improvement
          
    Returns
    -------
    result:list
        a list of endogenous reference gene
    '''
    lp=[]
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
    datana=pd.DataFrame()
    data1=pd.DataFrame()
    data2=pd.DataFrame()
    for i in itertools.combinations(range(1,depth+1), 2): 
        count=count+1 #calculate circle times
        data=data.dropna() 
        data=data.drop_duplicates(data.columns[0])
        data.reset_index(drop=True, inplace=True)
        
        last_std=pd.DataFrame()
        length=len(data)//1000
        remain=len(data)-1000*length
        for k in range(1,length+1):            
            datana=data.iloc[1000*(k-1):1000*k,0:1]
            data1=data.iloc[1000*(k-1):1000*k,i[0]:i[0]+1]
            data2=data.iloc[1000*(k-1):1000*k,i[1]:i[1]+1]
            l1=pd.DataFrame()
            l2=pd.DataFrame()
            for j in range(1000*(k-1),1000*k):
                l1[datana.loc[j]]=data1
                l2[datana.loc[j]]=data2
            l1=l1.div(np.asarray(data.iloc[1000*(k-1):1000*k,i[0]]))
            l2=l2.div(np.asarray(data.iloc[1000*(k-1):1000*k,i[1]]))           
            l=l1-l2
            l_std=l.std(axis=0)
            l_std=l_std.sort_values()[0:20]
            if(k==1):
                last_std=l_std
            else:
                last_std=pd.concat([last_std,l_std])
        datana=data.iloc[length*1000:length*1000+remain,0:1]
        data1=data.iloc[length*1000:length*1000+remain,i[0]:i[0]+1]
        data2=data.iloc[length*1000:length*1000+remain,i[1]:i[1]+1]
        l1=pd.DataFrame()
        l2=pd.DataFrame()
        for j in range(length*1000,length*1000+remain):
            l1[datana.loc[j]]=data1
            l2[datana.loc[j]]=data2
        l1=l1.div(np.asarray(data.iloc[length*1000:length*1000+remain,i[0]]))
        l2=l2.div(np.asarray(data.iloc[length*1000:length*1000+remain,i[1]]))       
        l=l1-l2
        l_std=l.std(axis=0)
        l_std=l_std.sort_values()[0:20]
        if(length>0):
            last_std=pd.concat([last_std,l_std])
        else:
            last_std=l_std            
        last_std=last_std.sort_values()[0:20]
        
        testlist=list(last_std.index)
        lp.append(testlist)
        #print(lllll)
        if(count==1):
            result=testlist
        if(count>1):
            result=list(set(testlist).intersection(set(result))) #Venn
    example = from_memberships(lp,data=range(len(lp)))
    end = time.time()
    print("calculate time:%.2fs"%(end-start))
    print(result)
    if depth>2:
        plot(example)
    


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
    l=data[data.iloc[:,0].isin([ERGname])]
    da=data.iloc[:,1:len(data.columns)]/np.asarray(l.iloc[0,1:len(data.columns)])
    da=da*l.iloc[0,1]
    dataname=data.iloc[:,0:1]
    da.insert(0,dataname.columns[0],dataname)
    return da

            
        
        
        
        

