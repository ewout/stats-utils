#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache 

import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from statlib import convert_fff, sasinput, stats,invlogit
import statsmodels.api as sm

#DATAFILE = '/home/ewout/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.txt'
#DICFILE =  '/home/ewout/enem/Microdados ENEM 2010/Input_SAS/INPUT_SAS_ENEM_2010.SAS'
#CSVFILE =  '/home/ewout/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.csv'
CSVFILE =  '/home/ewout/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.csv'

#gabarito2010 = {89:'BACEADCAECABDDACBABEBDEDAEECDBCDBDEEABDACEDDC',
#                90:'ACAEBCCEADADDBACBEBABDEEADEDDBCBCDDBAEEACDDEC',
#                91:'CABEAACDCEBADADCEBBBADDEAEEBDDCCBDBDAEEADCDCE',
#                92:'ABAECECCDADADABCABBBEDEEAEDBDCCBDDEEADBAECDCD',
#                105:'DBACCDADCCCEEEECAADCAEAACEADADDCADBEDDBBEAEAE',
#                106:'CDBCACADDCECEEECAADACEAACAECDADDDDABDEEBBAEEA',
#                107:'ACDCBCACDDEEECECCADAAEEACAADCDDAEDDADBEBABEAE',
#                108:'BACCDDACCDCEEEECACDAAEAECAAADCDDBDEDDAEBAAEBE'
#}


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

    #for qn,x in enumerate(a.T):
    #    df['AQ'+str(qn+1)] = x
        
    return df, itemstats, teststats, a
        
def resvec2(df,rescol='TX_RESPOSTAS_CN'):
    'Transforma o vetor de resolução em colunas do dataframe '
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
    #X = sm.add_constant(df['nota'],prepend=True)
    #lf = itemstats['iccfit'][qn-1].predict(X)
    const,sconst,nota,snota,itemd,sitemd = itemstats['iccfitsparam'][qn-1]
    ax.errorbar(hbin,prob,yerr=err,fmt='o')
    x = np.linspace(0.9*min(df['nota']),1.1*max(df['nota']),200)
    p = invlogit(const+nota*x)
    ax.plot(x,p,'g-')  
    ax.set_ylim(0,1)
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([200,400,600,800,1000])
    ax.text(0.03,0.85,'Q'+str(qn),transform = ax.transAxes)
    return fig, ax

def iccgrid(df,acertos,ncols=5,nrows=9):
    ''
    fig = plt.figure()
    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            fig, ax = iccgraph(df,acertos,qn,fig=fig,ax=ax)
            qn += 1
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.9,wspace=0.4,hspace=0.4)
    return fig

def iccfitgraph(df,acertos,fig=None):
    ''
    if not fig:
        fig = plt.figure()
        #fig.suptitle(u"Parámetros dos fits logisticos")

    ax1 = fig.add_subplot(211)
    ax1.set_title(u"Dificuldade")
    ax2 = fig.add_subplot(212)
    ax2.set_title(u"Discriminição")

    itemstats, teststats = stats(acertos,df['nota'])
    iccfitsparam = itemstats['iccfitsparam']

    itemds = [p[4] for p in iccfitsparam]
    err = [p[5] for p in iccfitsparam]
    x = np.arange(1,len(itemds)+1)
    ax1.errorbar(x,itemds,yerr=err,fmt='o')
    ax1.set_ylim(0,1400)
    ax1.set_xlim(0,48)
    ax1.set_xticklabels([])

    itemn = [p[2] for p in iccfitsparam]
    err = [p[3] for p in iccfitsparam]
    x = np.arange(1,len(itemn)+1)
    ax2.errorbar(x,itemn,yerr=err,fmt='o')
    ax2.set_xlabel(u"Questão")
    #ax1.set_ylim(0,1000)
    ax2.set_xlim(0,48)
    ax2.set_xticks([1,5,10,15,20,25,30,35,40,45])

    return fig   


def csv2df(idprov=89,tipprov='CN',sexo=None, raca=None):
    'Import enem csv to dataframe, clean it up.'
    csvfile = CSVFILE
    df = pandas.read_table(csvfile)
    if not idprov:
        df = df[df['NU_NT_'+tipprov] != '         ']
    else:
        df = df[(df['ID_PROVA_'+tipprov] == idprov) & (df['NU_NT_'+tipprov] != '         ')]


    if sexo:
        df = df[df['TP_SEXO'] == sexo]

    if raca:
        df = df[df['TP_COR_RACA'] == raca]
        
    df['nota'] = df['NU_NT_'+tipprov].apply(float)
    df, itemstats, teststats, acertos = resvec(df,'TX_RESPOSTAS_'+tipprov,'DS_GABARITO_'+tipprov)
    df = resvec2(df,rescol='TX_RESPOSTAS_'+tipprov)
    return df, itemstats, teststats, acertos


def iccgriddif(idprov=89,tipprov='CN',ncols=5,nrows=9):
    ''
    df, itemstats, teststats, acertos = csv2df(idprov,tipprov)
    dfm, itemstatsm, teststatsm, acertosm = csv2df(idprov,tipprov,sexo='M')
    dff, itemstatsf, teststatsf, acertosf = csv2df(idprov,tipprov,sexo='F')
    #dfm, itemstatsm, teststatsm, acertosm = csv2df(idprov,tipprov,raca='1')
    #dff, itemstatsf, teststatsf, acertosf = csv2df(idprov,tipprov,raca='2')
    fig = plt.figure()
    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            fig, ax = iccgraph(df,acertos,qn,fig=fig,ax=ax)
            fig, ax = iccgraph(dfm,acertosm,qn,fig=fig,ax=ax)
            fig, ax = iccgraph(dff,acertosf,qn,fig=fig,ax=ax)
            qn += 1
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.9,wspace=0.4,hspace=0.4)
    return fig

if __name__ == '__main__':
    df, itemstats, teststats, acertos = csv2df()
    
