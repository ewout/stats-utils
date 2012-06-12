#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache 

import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from statlib import convert_fff, sasinput

DATAFILE = '/home/ewout/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.txt'
DICFILE =  '/home/ewout/enem/Microdados ENEM 2010/Input_SAS/INPUT_SAS_ENEM_2010.SAS'
CSVFILE =  '/home/ewout/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.csv'

gabarito2010 = {89:'BACEADCAECABDDACBABEBDEDAEECDBCDBDEEABDACEDDC',
                90:'ACAEBCCEADADDBACBEBABDEEADEDDBCBCDDBAEEACDDEC',
                91:'CABEAACDCEBADADCEBBBADDEAEEBDDCCBDBDAEEADCDCE',
                92:'ABAECECCDADADABCABBBEDEEAEDBDCCBDDEEADBAECDCD',
                105:'DBACCDADCCCEEEECAADCAEAACEADADDCADBEDDBBEAEAE',
                106:'CDBCACADDCECEEECAADACEAACAECDADDDDABDEEBBAEEA',
                107:'ACDCBCACDDEEECECCADAAEEACAADCDDAEDDADBEBABEAE',
                108:'BACCDDACCDCEEEECACDAAEAECAAADCDDBDEDDAEBAAEBE'
}

def stats(acertos):
    ''
    teststats = {}
    itemstats = {}
    N = acertos.shape[0]
    k = acertos.shape[1]
    itemf = acertos.sum(axis=0) / float(N)
    itemv = acertos.var(axis=0)
    scores = acertos.sum(axis=1)

    alpha = 1.0*k/(k-1)*(1-sum(itemv)/scores.var())

    teststats['alpha'] = alpha
    teststats['scores'] = scores
    itemstats['N'] = N
    itemstats['k'] = k
    itemstats['itemf'] = itemf
    itemstats['itemv'] = itemv
    
    return itemstats, teststats 

def resvec(df,rescol,gabcol):
    ''
    res = df[rescol]
    gab = df[gabcol]
    l = []
    for rvec,gvec in zip(res,gab):
        l.append([1 if x==y else 0 for x,y in zip(rvec,gvec)])

    a = np.array(l)
    itemstats, teststats = stats(a)
    df['res'] = list(a)
    df['ressum'] = a.sum(axis=1)
    df['resstd'] = a.std(axis=1)

    for qn,x in enumerate(a.T):
        df['AQ'+str(qn+1)] = x

        
    return df, itemstats, teststats
        
def resvec2(df,rescol='TX_RESPOSTAS_CN'):
    'Transforma o vetore de resolução em colunas do dataframe '
    res = df[rescol]
    l = []
    for rvec in res:
        l.append(list(rvec))
    a = np.array(l)
    for qn,x in enumerate(a.T):
        df['Q'+str(qn+1)] = x
    return df

def questionstats(df,qn):
    ''
    q = df['Q'+str(qn)]
    correct = gabarito2010[idprova][qn-1]
    

def gradebar(df,qn,fig=None,ax=None):
    ''
    idprova = 89
    q = df['Q'+str(qn)]
    correct = gabarito2010[idprova][qn-1]
    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
    c = Counter(q)
    labels = sorted(c.keys())
    correctindex = labels.index(correct)
    values = [c[val] for val in labels]
    N = len(labels)
    x = np.arange(N)
    width = 0.6
    ax.bar(x,values,width,color='r')
    ax.bar([correctindex],values[correctindex],width,color='b')
    #ax.bar(x,values,width,color='b')
    ax.set_xticks(x+0.5*width)
    ax.set_xticklabels(labels)
    return ax

def gradegrid(df,ncols=5,nrows=9):
    ''
    fig = plt.figure()
    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            ax = gradebar(df,qn,fig=fig,ax=ax)
            qn += 1
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.9,wspace=0.4,hspace=0.4)
    return fig

def rawdata2csv(dadosfile = DATAFILE, dicfile = DICFILE, outfile = None):
    'raw inep fixed format to csv'
    dic = sasinput(dicfile,filtercols=['NU_INSCRICAO','ID_PROVA_CN','NU_NT_CN','TX_RESPOSTAS_CN','DS_GABARITO_CN'])
    print dic
    if not outfile:
        outfile = dadosfile[:-3] + 'csv'
    convert_fff(dadosfile,outfile,dic,sample = 0.001)

def csv2df(csvfile = CSVFILE,idprov=89,tipprov='CN'):
    'Import enem csv to dataframe, clean it up.'
    df = pandas.read_table(csvfile)
    df = df[(df['ID_PROVA_'+tipprov] == idprov) & (df['NU_NT_'+tipprov] != '         ')]
    df['nota'] = df['NU_NT_'+tipprov].apply(float)
    df, teststats, itemstats = resvec(df,'TX_RESPOSTAS_'+tipprov,'DS_GABARITO_'+tipprov)
    df = resvec2(df)
    return df, itemstats, teststats


if __name__ == '__main__':
    df, itemstats, teststats = csv2df(CSVFILE)
    
