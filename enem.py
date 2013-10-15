#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache

import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import Counter
from statlib import convert_fff, sasinput, stats,invlogit
import statsmodels.api as sm
import os

CSVFILE2009 =  '/home/ewout/enem/2009/dados/enem_2009_1.csv'
CSVFILE2010 =  '/home/ewout/enem/2010/dados/enem_2010_1.csv'
CSVFILE2011 =  '/home/ewout/enem/2011/dados/enem_2011_1.csv'

def resvec(df,rescol,gabcol,hscale=None):
    ''

    dmap = {'A':1,
            'B':2,
            'C':3,
            'D':4,
            'E':5,
            '.':'NA',
            '*':'NA'
            }
    tonumbers = lambda x: dmap[x]


    res = df[rescol]
    gab = df[gabcol]
    l = []
    l1 = []
    for rvec,gvec in zip(res,gab):
        l.append([1 if x==y else 0 for x,y in zip(rvec,gvec)])
        l1.append(map(tonumbers,rvec))
    a = np.array(l)
    an = np.array(l1)
    df['res'] = list(a)
    df['gab'] = gab
    df['ressum'] = a.sum(axis=1)
    df['resstd'] = a.std(axis=1)

    itemstats, teststats = stats(a,hscale,df)

    return df, itemstats, teststats, a, an

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

def orderedfig(itemlabels,itemvalues,itemlabels9,itemvalues9,maxitems=10,fig=None,ax=None):
    'Faz um gráfico de itens em ordem (reverso) de qualidade'

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
    width = 0.5
    ax.bar(np.arange(1,maxitems+1),itemvalues,width=width,align='center',color='b',label="2010")
    ax.bar(np.arange(1+width,maxitems+1+width),itemvalues9,width=width,align='center',color="g",label="2009")

    ax.legend(loc=2)
    ax.set_xlim(0,maxitems+1)
    ax.set_xticks([])
    trans = ax.get_xaxis_transform()
    for i,label in enumerate(itemlabels):
        ax.text(i+1,-0.04,label,clip_on=False,color='b',rotation="horizontal",ha='center',weight='bold',transform=trans)

    for i,label in enumerate(itemlabels9):
        ax.text(i+1+width,-0.06,label,clip_on=False,color='g',rotation="horizontal",ha='center',weight='bold',transform=trans)

    return fig,ax


def itemfbar(acertos,acertos2,order=True,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    itemf = itemstats['itemf']

    itemstats2, teststats2 = stats(acertos2)
    itemf2 = itemstats2['itemf']


    if order:
        itemfdf = pandas.DataFrame(itemstats['itemf'],index=range(1,len(itemf)+1))
        itemfdf = itemfdf.sort(columns=0)
        itemfdf2 = pandas.DataFrame(itemstats2['itemf'],index=range(1,len(itemf2)+1))
        itemfdf2 = itemfdf2.sort(columns=0)


    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"ENEM 2010")
    width = 0.5
    if order:
        ax.bar(np.arange(1,len(itemf)+1),itemfdf[0],width=width,align='center',color='b',label="ENEM 2010")
        ax.bar(np.arange(1+width,len(itemf)+1+width),itemfdf2[0],width=width,align='center',color="g",label="ENEM 2009")
    else:
        ax.bar(np.arange(1,len(itemf)+1),itemf,width=width,align='center')

    ax.plot([0,len(itemf)+1],[0.2,0.2],'b',label="Chute",linewidth=2)
    ax.plot([0,len(itemf)+1],[itemf.mean(),itemf.mean()],'k',label=u"Média",linewidth=2)
    #    m10 = itemf.mean()
    #m9 = itemf2.mean()
    #ax.annotate(u"Médias",xy=(0,m10),xytext=(10,m10),color='k',arrowprops=dict(facecolor='b',width=2,shrink=0.05))
    #ax.annotate("",xy=(0,m9),xytext=(10,m9),color='g',arrowprops=dict(facecolor='g',width=2,shrink=0.05))
    #ax.annotate("Chute",xy=(0,0.2),xytext=(10,0.25),color='k',arrowprops=dict(facecolor='grey',width=1,shrink=0.05))
    ax.legend()
    ax.set_ylabel(u"Fração dos alunos que acertaram a questão")
    ax.set_xlabel(u"Item")
    ax.set_ylim(0,1)
    ax.set_xlim(0,len(itemf)+1)
    if order:
        ax.set_xticks([])
    else:
        ax.set_xticks([1,5,10,15,20,25,30,35,40,45])
    return fig,ax

def itemfbar2(acertos,acertos2,maxitems=10,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    itemf = itemstats['itemf']

    itemstats2, teststats2 = stats(acertos2)
    itemf2 = itemstats2['itemf']


    itemfdf = pandas.DataFrame(itemstats['itemf'],index=range(1,len(itemf)+1))
    itemfdf = itemfdf.sort(columns=0)
    itemfdf = itemfdf[0:maxitems]
    itemfdf2 = pandas.DataFrame(itemstats2['itemf'],index=range(1,len(itemf2)+1))
    itemfdf2 = itemfdf2.sort(columns=0)
    itemfdf2 = itemfdf2[0:maxitems]


    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Índice de Dificuldade")

    ax.set_ylabel(u"Fração dos alunos que acertaram a questão")
    ax.set_ylim(0,0.25)
    fig,ax = orderedfig(itemfdf.index,itemfdf[0],itemfdf2.index,itemfdf2[0],maxitems,fig,ax)
    ax.text(0.5,-0.1,u"Itens em ordem de dificuldade",clip_on=False,transform = ax.transAxes,ha='center')

    orderedstats = itemfdf.reset_index()
    orderedstats.columns = ['fQ2010','f2010']

    orderedstats2 = itemfdf2.reset_index()
    orderedstats2.columns = ['fQ2009','f2009']

    return fig,ax, orderedstats.join(orderedstats2)


def idbar(acertos,acertos2,order=True,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    itemstats2, teststats2 = stats(acertos2)

    #id50 = itemstats['id50']
    id27 = itemstats['id27']
    id279 = itemstats2['id27']
    if order:
        iddf = pandas.DataFrame(id27,index=range(1,len(id27)+1))
        iddf = iddf.sort(columns=0)
        iddf9 = pandas.DataFrame(id279,index=range(1,len(id279)+1))
        iddf9 = iddf9.sort(columns=0)



    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Índice de Discriminação")
    width = 0.5
    if order:
        ax.bar(np.arange(1,len(id27)+1),iddf[0],width=width,align='center',color='b',label=u'2010')
        ax.bar(np.arange(1+width,len(id279)+1+width),iddf9[0],width=width,align='center',color='g',label="2009")
    else:
        ax.bar(np.arange(1,len(id27)+1),id27,width=width,align='center',color='b',label=u'diferença de acertos entre os piores e melhores 27%')
    #ax.text(u'diferença de acertos entre os piores e melhores 27%')
    ax.legend(loc="upper left")
    ax.set_xlabel(u"Itens em ordem de discriminação")
    ax.set_ylim(-0.05,0.4)
    ax.set_xlim(0,len(id27)+1)
    if order:
        ax.set_xticks([])
    else:
        ax.set_xticks([1,5,10,15,20,25,30,35,40,45])
    return fig,ax

def idbar2(acertos,acertos2,maxitems=10,fig=None,ax=None):
    ''
    itemstats, teststats = stats(acertos)
    itemstats2, teststats2 = stats(acertos2)

    id27 = itemstats['id25']
    id279 = itemstats2['id25']

    iddf = pandas.DataFrame(id27,index=range(1,len(id27)+1))
    iddf = iddf.sort(columns=0)
    iddf = iddf[0:maxitems]
    iddf9 = pandas.DataFrame(id279,index=range(1,len(id279)+1))
    iddf9 = iddf9.sort(columns=0)
    iddf9 = iddf9[0:maxitems]

    orderedstats = iddf.reset_index()
    orderedstats.columns = ['idQ2010','id2010']

    orderedstats2 = iddf9.reset_index()
    orderedstats2.columns = ['idQ2009','id2009']

    orderedstats = orderedstats.join(orderedstats2)


    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Índice de Discriminação (quartis)")


    fig,ax = orderedfig(iddf.index,iddf[0],iddf9.index,iddf9[0],maxitems,fig,ax)
    ax.set_xlabel(u"")
    ax.set_ylim(-0.01,0.15)

    ax.text(0.5,-0.1,u"Item",clip_on=False,transform = ax.transAxes,ha='center')



    return fig,ax, orderedstats

def biscorr(maxitems=10,fig=None,ax=None):
    ''
    from string import lstrip

    biscorr49 = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc49-biscorr.csv',header=None,names=['Q','biscorr'])
    biscorr89 = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc89-biscorr.csv',header=None,names=['Q','biscorr'])
    biscorr121 = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc121-biscorr.csv',header=None,names=['Q','biscorr'])


    stripV = lambda s: lstrip(s,'V')
    biscorr49['Q'] = biscorr49['Q'].apply(stripV)
    biscorr89['Q'] = biscorr89['Q'].apply(stripV)
    biscorr121['Q'] = biscorr121['Q'].apply(stripV)
    biscorr49 = biscorr49.sort(columns='biscorr')
    biscorr89 = biscorr89.sort(columns='biscorr')
    biscorr121 = biscorr121.sort(columns='biscorr')
    biscorr49 = biscorr49[0:maxitems]
    biscorr89 = biscorr89[0:maxitems]
    biscorr121 = biscorr121[0:maxitems]


    ostats = biscorr89.copy()
    ostats.index = range(maxitems)
    ostats.columns = ['bcQ2010','bc2010']
    ostats2 = biscorr49.copy()
    ostats2.index = range(maxitems)
    ostats2.columns = ['bcQ2009','bc2009']
    ostats3 = biscorr121.copy()
    ostats3.index = range(maxitems)
    ostats3.columns = ['bcQ2011','bc2011']

    ostats = ostats.join(ostats2)
    ostats = ostats.join(ostats3)

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Correlação biserial")

    fig,ax = orderedfig(biscorr89['Q'],biscorr89['biscorr'],biscorr49['Q'],biscorr49['biscorr'],maxitems,fig,ax)

    ax.set_xlabel(u"")
    ax.set_ylim(0,0.2)

    ax.set_xticks([])
    ax.text(0.5,-0.1,u"Item",clip_on=False,transform = ax.transAxes,ha='center')


    return fig,ax,ostats

def logfit(itemstats,itemstats2,maxitems=10,fig=None,ax=None):
    ''
    N = itemstats['k']
    disc = np.array(itemstats['iccfitsparam'])[:,2]
    disc9 = np.array(itemstats2['iccfitsparam'])[:,2]
    disc = pandas.DataFrame(disc,index=range(1,N+1)).sort(columns=0)[0:maxitems]
    disc9 = pandas.DataFrame(disc9,index=range(1,N+1)).sort(columns=0)[0:maxitems]

    ostats = disc.reset_index()
    ostats.columns = ['Q2010','logfit2010']

    ostats2 = disc9.reset_index()
    ostats2.columns = ['Q2009','logfit2009']

    ostats = ostats.join(ostats2)


    fig,ax = orderedfig(disc.index,disc[0],disc9.index,disc9[0])
    ax.set_title(u"Discriminação via CCI empírica")

    return fig,ax, ostats

def gradebar(df,qn,fig=None,ax=None):
    'qn from 1 til 45'
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
    values = [c[val]*1.0/len(df) for val in labels]
    N = len(labels)
    x = np.arange(N)
    width = 0.6
    ax.bar(x,values,width,color='r')
    ax.bar([correctindex],values[correctindex],width,color='b')
    ax.set_xticks(x+0.5*width)
    ax.set_xticklabels(labels)
    ax.set_ylim(0,0.5)
    ax.yaxis.set_major_locator(MaxNLocator(2))
    #ax.yaxis.set_minor_locator(MaxNLocator(2))
    ax.text(0.02,0.8,'Q'+str(qn),transform = ax.transAxes)
    return fig, ax

def gradegrid(df,ncols=5,nrows=9,fig=None):
    ''
    if not fig:
        fig = plt.figure()

    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            ax = gradebar(df,qn,fig=fig,ax=ax)
            qn += 1
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.9,wspace=0.4,hspace=0.4)
    return fig

def iccgraph(df,acertos,qn,hs = "scores",fig=None,ax=None):
    ''
    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(111)
        ax.set_title(u"Curva Característica do Item "+str(qn))

    itemstats, teststats = stats(acertos,hs = hs,df=df)
    hscale = itemstats['hscale']
    icc = itemstats['icc'][qn-1]
    hbin = icc[:,0]
    nbin = icc[:,1]
    acertos_no_bin = icc[:,2]
    prob = icc[:,3]
    err = icc[:,4]

    const,sconst,nota,snota,itemd,sitemd = itemstats['iccfitsparam'][qn-1]
    ax.errorbar(hbin,prob,yerr=err,fmt='o')
    x = np.linspace(0.9*min(hscale),1.1*max(hscale),200)
    p = invlogit(const+nota*x)
    ax.plot(x,p,'g-')
    ax.set_ylim(0,1)
    ax.set_yticks([0,0.5,1])
    if hs == 'scores':
        ax.set_xlabel(u"Acertos")
    else:
        ax.set_xlabel(u"Escore Enem")
    ax.set_ylabel(u"Probabilidade")

    return fig, ax

def ltmfitparams(provid):
    ''

    if provid == 89:
        ltmsummary = pandas.read_table('/home/ewout/Dropbox/RIRT/ltm89summary.csv')
    elif provid == 49:
        ltmsummary = pandas.read_table('/home/ewout/Dropbox/RIRT/ltm49summary.csv')
    elif provid == 121:
        ltmsummary = pandas.read_table('/home/ewout/Dropbox/RIRT/ltm121summary.csv')
    else:
        raise Exception('Só usar com prova 49,89 ou 121')
    ltmdif = ltmsummary[:45]
    ltmdisc = ltmsummary[45:]
    ltmdif.index = range(1,46)
    ltmdisc.index = range(1,46)

    return ltmdif,ltmdisc

def ltmfitgraph(df,acertos,qn1,qn2,fig=None):
    'Só usar com prova 89 ou 49!'

    provid = df['ID_PROVA_CN'].values[0]

    if not fig:
        fig = plt.figure()
        fig.suptitle(u"Curvas Características dos Itens "+str(qn1)+" e "+str(qn2)+" (+ ajuste 2PL TRI)")
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for qn,ax in [(qn1,ax1),(qn2,ax2)]:
        itemstats, teststats = stats(acertos,hs = 'notapadrao',df=df)
        hscale = itemstats['hscale']
        icc = itemstats['icc'][qn-1]
        hbin = icc[:,0]
        nbin = icc[:,1]
        acertos_no_bin = icc[:,2]
        prob = icc[:,3]
        err = icc[:,4]

        ax.errorbar(hbin,prob,yerr=err,fmt='o')

        x = np.linspace(0.9*min(hscale),1.1*max(hscale),200)
        ltmdif,ltmdisc = ltmfitparams(provid)
        a = ltmdisc['value'][qn]
        b = ltmdif['value'][qn]
        p = invlogit(1.0*a*(x-b))
        ax.plot(x,p,'k-')
        ax.set_ylim(0,1)
        ax.set_yticks([0,0.5,1])
        ax.set_xlabel(u"Escore Enem Padronizada")
        ax.set_ylabel(u"Probabilidade")

    ax2.set_ylabel('')
    ax1.text(0.05,0.9,'Q'+str(qn1),transform = ax1.transAxes,fontsize='medium',weight='bold')
    ax2.text(0.05,0.9,'Q'+str(qn2),transform = ax2.transAxes,fontsize='medium',weight='bold')


    return fig

def ltmgrid(df,acertos,ncols=5,nrows=9,fig=None):
    ''
    if not fig:
        fig = plt.figure()
    qn = 1
    provid = df['ID_PROVA_CN'].values[0]
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            itemstats, teststats = stats(acertos,hs = 'notapadrao',df=df)
            hscale = itemstats['hscale']
            icc = itemstats['icc'][qn-1]
            hbin = icc[:,0]
            nbin = icc[:,1]
            acertos_no_bin = icc[:,2]
            prob = icc[:,3]
            err = icc[:,4]

            ax.errorbar(hbin,prob,yerr=err,fmt='o')

            x = np.linspace(0.9*min(hscale),1.1*max(hscale),200)
            ltmdif,ltmdisc = ltmfitparams(provid)
            a = ltmdisc['value'][qn]
            b = ltmdif['value'][qn]
            p = invlogit(1.0*a*(x-b))
            ax.plot(x,p,'k-')
            ax.set_ylim(0,1)
            ax.set_yticks([0,0.5,1])

            qn += 1
        fig.subplots_adjust(left=0.1,right=0.95,bottom=0.05,top=0.9,wspace=0.4,hspace=0.4)
    return fig


def iccgrid(df,acertos,ncols=5,nrows=9,hs='scores',fig=None):
    ''
    if not fig:
        fig = plt.figure()
    qn = 1
    for row in range(nrows):
        for col in range(ncols):
            ax = plt.subplot2grid((nrows,ncols),(row,col))
            fig, ax = iccgraph(df,acertos,qn,hs,fig=fig,ax=ax)
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

    itemstats, teststats = stats(acertos,'nota',df)
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


def csv2df(idprov=89,tipprov='CN',sexo=None, raca=None,hscale=None):
    'Import enem csv to dataframe, clean it up.'
    if idprov in range(49,85):
        csvfile = CSVFILE2009
    elif idprov in range(89,117):
        csvfile = CSVFILE20105percent
    elif idprov in range(121,138):
        csvfile = CSVFILE2011
    else:
        raise Exception("ID da Prova errada!")
    df = pandas.read_table(csvfile)
    print "number of rows in df before filters:", len(df)
    df = df[df['IN_PRESENCA_'+tipprov] == 1]
    print "number of rows after no-show filter (",tipprov,") :",len(df)
    if not idprov:
        df = df[df['NU_NT_'+tipprov] != '         ']
    else:
        # no enem 2009 o pandas não reconhece que NU_ID_PROVA é integer...
        df['ID_PROVA_'+tipprov] = df['ID_PROVA_'+tipprov].apply(int)
        df = df[(df['ID_PROVA_'+tipprov] == idprov) & (df['NU_NT_'+tipprov] != '         ')]

    if sexo is not None:
        df = df[df['TP_SEXO'] == sexo]

    if raca is not None:
        df = df[df['TP_COR_RACA'] == raca]


    print "number of rows in df after filters:", len(df)
    df['nota'] = df['NU_NT_'+tipprov].apply(float)
    df, itemstats, teststats, acertos, acertosn = resvec(df,'TX_RESPOSTAS_'+tipprov,'DS_GABARITO_'+tipprov,hscale=hscale)
    df = resvec2(df,rescol='TX_RESPOSTAS_'+tipprov)
    return df, itemstats, teststats, acertos, acertosn


def iccgriddif(idprov=89,tipprov='CN',ncols=5,nrows=9):
    ''
    df, itemstats, teststats, acertos, acertosn = csv2df(idprov,tipprov)
    dfm, itemstatsm, teststatsm, acertosm, acertosmn = csv2df(idprov,tipprov,sexo=0)
    dff, itemstatsf, teststatsf, acertosf, acertosfn = csv2df(idprov,tipprov,sexo=1)
    #dfm, itemstatsm, teststatsm, acertosm, acertosmn = csv2df(idprov,tipprov,sexo='M')
    #dff, itemstatsf, teststatsf, acertosf, acertosfn = csv2df(idprov,tipprov,sexo='F')
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

def allstats(idprov):
    ''
    df, itemstats, teststats, acertos, acertosn = csv2df(idprov=idprov,tipprov='CN',hscale='scores')
    itemf = itemstats['itemf']
    id25 = itemstats['id25']

    from string import lstrip
    if idprov == 49:
        biscorr = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc49-biscorr.csv',header=None,names=['Q','biscorr'])
        biscorr = biscorr['biscorr'].values
    elif idprov == 89:
        biscorr = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc89-biscorr.csv',header=None,names=['Q','biscorr'])
        biscorr = biscorr['biscorr'].values
    elif idprov == 121:
        biscorr = pandas.read_table('/home/ewout/Dropbox/RIRT/dsc121-biscorr.csv',header=None,names=['Q','biscorr'])
        biscorr = biscorr['biscorr'].values

    else:
        biscorr = "Ainda precisa calcular os coefs bisserial para a prova %i" % idprov



    logfit = np.array(itemstats['iccfitsparam'])[:,2]
    ltmdif,ltmdisc = ltmfitparams(idprov)
    ltmdif = ltmdif['value'].values
    ltmdisc = ltmdisc['value'].values

    statdf = pandas.DataFrame({'itemf':itemf,'id25':id25,'biscorr':biscorr,'logfit':logfit,'ltmdisc':ltmdisc,'ltmdif':ltmdif})

    return statdf

def generate_graphs(graphs = "all"):
    ''
    df, itemstats, teststats, acertos, acertosn = csv2df(idprov=89,tipprov='CN',hscale='scores')
    df9, itemstats9, teststats9, acertos9, acertosn9 = csv2df(idprov=49,tipprov='CN',hscale='scores')

    igr = 2/(1+np.sqrt(5))
    pol = 2.54
    swidth =  8/pol
    sfig = (swidth,igr*swidth)
    dfig = (2*swidth,2*igr*swidth)


    #fig = plt.figure(figsize=dfig)
    fig, ax, ostats = itemfbar2(acertos,acertos9)
    fig.savefig('../figs/itemfbar.png')
    # bar graph da questão 25 de 2010
    #fig = plt.figure(figsize=dfig)
    fig, ax = gradebar(df,25)
    fig.savefig('../figs/gradebar-25-2010-5percent.png')

    fig, ax = gradebar(df9,31)
    fig.savefig('../figs/gradebar-31-2009.png')

    fig, ax, ostats2 = idbar2(acertos,acertos9)
    fig.savefig('../figs/idbar-5percent.png')


    fig, ax, ostats3 = biscorr()
    fig.savefig('../figs/biscorr.png')

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    iccgraph(df,acertos,11,hs='scores',fig=fig,ax=ax1)
    iccgraph(df,acertos,25,hs='scores',fig=fig,ax=ax2)
    ax2.set_ylabel('')
    ax1.text(0.05,0.9,'Q11',transform = ax1.transAxes,fontsize='medium',weight='bold')
    ax2.text(0.05,0.9,'Q25',transform = ax2.transAxes,fontsize='medium',weight='bold')
    fig.savefig('../figs/eicc-5percent.png')

    fig = plt.figure()
    fig = ltmfitgraph(df,acertos,11,25,fig=fig)
    fig.savefig('../figs/ltm89-11-25-5percent.png')

    fig, ax, ostats4 = logfit(itemstats,itemstats9)
    fig.savefig('../figs/logfit-5percent.png')

    ostatstotal = ostats.join(ostats2)
    ostatstotal = ostatstotal.join(ostats3)
    ostatstotal = ostatstotal.join(ostats4)

    ostatstotal.to_excel('../figs/ostats-5percent.xlsx')

    return ostatstotal


def generate_graphs2(idprov,tipprov='CN'):
    ''
    figsize = (6,6)
    df, itemstats, teststats, acertos, acertosn = csv2df(idprov,tipprov,hscale='scores')
    # padrões de resposta
    hs1,hs2,hs3,hs4 = [],[],[],[]
    for qn in range(1,46):
        fn = 'padr-resposta-q'+str(qn)+'.png'
        fig = plt.figure(figsize=figsize)
        print fn
        fig, ax = gradebar(df,qn,fig)
        path = os.path.join("../figs/",tipprov,str(idprov),fn)
        fig.savefig(path)
        path = path[3:]
        hs1.append("<a href=\"%s\"><img src=\"%s\" alt=\"\" /></a>" % (path,path))

        fn = 'cce-escore-total-q'+str(qn)+'.png'
        fig = plt.figure(figsize=figsize)
        print fn
        fig, ax = iccgraph(df,acertos,qn,hs="scores",fig=fig)
        path = os.path.join("../figs/",tipprov,str(idprov),fn)
        fig.savefig(path)
        path = path[3:]
        hs2.append("<a href=\"%s\"><img src=\"%s\" alt=\"\" /></a>" % (path,path))


        fn = 'cce-escore-enem-q'+str(qn)+'.png'
        fig = plt.figure(figsize=figsize)
        print fn
        fig, ax = iccgraph(df,acertos,qn,hs="notapadrao",fig=fig)
        path = os.path.join("../figs/",tipprov,str(idprov),fn)
        fig.savefig(path)
        path = path[3:]
        hs3.append("<a href=\"%s\"><img src=\"%s\" alt=\"\" /></a>" % (path,path))

    path = os.path.join("../figs/",tipprov,str(idprov),"hmtlsnippet1.html")
    open(path,'w').write("\n".join(hs1))


    path = os.path.join("../figs/",tipprov,str(idprov),"hmtlsnippet2.html")
    open(path,'w').write("\n".join(hs2))

    path = os.path.join("../figs/",tipprov,str(idprov),"hmtlsnippet3.html")
    open(path,'w').write("\n".join(hs3))

    for qn in range(1,46):
        fn = 'occ-q'+str(qn)+'.png'
        path = os.path.join("figs/",tipprov,str(idprov),fn)
        hs4.append("<a href=\"%s\"><img src=\"%s\" alt=\"\" /></a>" % (path,path))

    path = os.path.join("../figs/",tipprov,str(idprov),"hmtlsnippet4.html")
    open(path,'w').write("\n".join(hs4))

if __name__ == '__main__':
    generate_graphs()
