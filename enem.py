#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache 

import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from statlib import convert_fff, sasinput, stats

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


def resvec(df,rescol,gabcol):
    ''
    res = df[rescol]
    gab = df[gabcol]
    l = []
    for rvec,gvec in zip(res,gab):
        l.append([1 if x==y else 0 for x,y in zip(rvec,gvec)])

    a = np.array(l)
    itemstats, teststats = stats(a,df['nota'])
    df['res'] = list(a)
    df['gab'] = gab
    df['ressum'] = a.sum(axis=1)
    df['resstd'] = a.std(axis=1)

    for qn,x in enumerate(a.T):
        df['AQ'+str(qn+1)] = x
        
    return df, itemstats, teststats, a
        
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

def itemfbar(acertos,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    itemf = itemstats['itemf']
    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Fração dos alunos que acertaram a questão")
    width = 0.8
    ax.bar(np.arange(1,len(itemf)+1),itemf,width=width,align='center')
    ax.plot([0,len(itemf)+1],[0.2,0.2],'b',label="Chute",linewidth=2)
    ax.plot([0,len(itemf)],[itemf.mean(),itemf.mean()],'k',label=u"Média",linewidth=2)
    ax.legend()
    ax.set_xlabel(u"Questão")
    ax.set_ylim(0,1)
    ax.set_xlim(0,48)
    ax.set_xticks([1,5,10,15,20,25,30,35,40,45])
    return fig,ax

def idbar(acertos,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    id50 = itemstats['id50']
    id27 = itemstats['id27']

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Índice de Discriminação")
    width = 0.8
    ax.bar(np.arange(1,len(id27)+1),id27,width=width,align='center',color='b',label=u'diferença de acertos entre os piores e melhores 27%')
    ax.bar(np.arange(1,len(id50)+1),id50,width=width,align='center',color='r',alpha=0.5,label='50%')
    ax.set_xlabel(u"Questão")
    ax.set_ylim(-0.05,1)
    ax.set_xlim(0,48)
    ax.set_xticks([1,5,10,15,20,25,30,35,40,45])
    ax.legend()
    return fig,ax



def gradebar(df,qn,fig=None,ax=None):
    ''
    #idprova = 89
    q = df['Q'+str(qn)]
    #correct = gabarito2010[idprova][qn-1]
    correct = df['gab'].values[0][qn-1]
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
    ax.set_xticks(x+0.5*width)
    ax.set_xticklabels(labels)
    ax.text(0.02,0.8,'Q'+str(qn),transform = ax.transAxes)
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

def iccgraph(df,acertos,qn,fig=None,ax=None):
    ''
    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Curva Característica Questão "+str(qn))

    itemstats, teststats = stats(acertos,hscale = df['nota'])
    icc = itemstats['icc'][qn-1]
    hbin = icc[:,0]
    nbin = icc[:,1]
    acertos_no_bin = icc[:,2]
    prob = icc[:,3]
    err = icc[:,4]
    ax.errorbar(hbin,prob,yerr=err,fmt='o')
    ax.set_ylim(0,1)
    ax.set_yticks([0,0.5,1])
    ax.text(0.03,0.85,'Q'+str(qn),transform = ax.transAxes)
    return fig, ax

def iccgrid(df,acertos,ncols=5,nrows=9):
    ''
    fig = plt.figure()
    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            ax = iccgraph(df,acertos,qn,fig=fig,ax=ax)
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

def csv2df(idprov=89,tipprov='CN'):
    'Import enem csv to dataframe, clean it up.'
    #if tipprov == 'CN':
    #    csvfile = CSVFILE[:-4] + '-CN.csv'
    #elif tipprov == 'CH':
    #    csvfile = CSVFILE[:-4] + '-CH.csv'
    csvfile = CSVFILE
    df = pandas.read_table(csvfile)
    df = df[(df['ID_PROVA_'+tipprov] == idprov) & (df['NU_NT_'+tipprov] != '         ')]
    df['nota'] = df['NU_NT_'+tipprov].apply(float)
    df, itemstats, teststats, acertos = resvec(df,'TX_RESPOSTAS_'+tipprov,'DS_GABARITO_'+tipprov)
    df = resvec2(df,rescol='TX_RESPOSTAS_'+tipprov)
    return df, itemstats, teststats, acertos


if __name__ == '__main__':
    df, itemstats, teststats, acertos = csv2df(CSVFILE)
    
