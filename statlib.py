#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache 

from __future__ import division
import sys, os, random
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pandas
from collections import Counter

import statsmodels.api as sm

def isinteger(x):
    ''
    ok = True
    try:
        int(x)
    except ValueError:
        return False
    return ok

def sasinput(infile,filtercols=[]):
    ''
    infile = os.path.expanduser(infile)
    colnames = []
    quebras = []
    for line in open(infile):
        if line.startswith('@'):
            cut,name,lenspec,comment = line.split(None,3)
            quebras.append(int(cut[1:]))
            colnames.append(name)
    try:
        lenspec = int(lenspec.split('.')[0])
    except ValueError:
        lenspec = int(lenspec[5:].split('.')[0])
    quebras.append(quebras[-1]+lenspec)
    dic = []
    if not filtercols:
        filtercols = colnames
    for i,name in enumerate(colnames):
        if name in filtercols:
            dic.append((name,(quebras[i],quebras[i+1])))
                
    return dic


def l2u(s):
    'converte string latin1 em utf8'
    return s.decode('latin1').encode('utf8')

def slices(s, dic):
    l = []
    for colname,cut in dic:
        piece = s[cut[0]-1:cut[1]-1]
        l.append(piece)
    return l

def filelen(fn):
    'number of line in file for files with fixed-width lines'
    fn = os.path.expanduser(fn)
    fsize = os.path.getsize(fn)
    linesize = len(open(fn).readline())
    return fsize // linesize


def convert_fff(infile,outfile,dic,max=sys.maxint,sample=1):
    'sample = fraction of lines to sample'
    infile = os.path.expanduser(infile)
    outfile = os.path.expanduser(outfile)
    
    outfile = open(outfile,'w')
    colnames = [col for col,cuts in dic]
    outfile.write('\t'.join(colnames)+'\n')

    flines = filelen(infile)
    if sample < 1:
        print flines, int(round(sample*flines))
        linestoconvert = random.sample(xrange(flines),int(round(sample*flines)))
        print "converting only %i lines." % len(linestoconvert)
    else:
        linestoconvert = range(flines)

    linesdone = 0
    for lnumber, line in enumerate(open(infile)):
        if lnumber <= max:
            if lnumber in linestoconvert:
                print "converting line %i" % lnumber, "(%i lines done)" % linesdone 
                linesdone += 1
                line = slices(line,dic)
                line = '\t'.join(line)
                line = l2u(line)
                outfile.write(line+'\n')
    outfile.close()


def itemdiscrimination(acertos,frac=0.5,scores=None):
    'If scores is given, use it as an ability scale, otherwise use acertos.'
    if not scores:
        scores = acertos.sum(axis=1)
    low = int(round((frac)*len(scores)))
    high = int(round((1-frac)*len(scores)))
    sscores = sorted(scores)
    vhigh,vlow = sscores[high],sscores[low]
    highgroup,lowgroup = acertos[scores > vhigh], acertos[scores < low]
    Nhigh,Nlow = highgroup.shape[0],lowgroup.shape[0]
    phigh,plow =  highgroup.sum(axis=0) / Nhigh, lowgroup.sum(axis=0) / Nlow

    return phigh - plow

def icc(acertos,qn,hscale,bins=10):
    ''
    acertos = acertos[:,qn]
    from collections import Iterable
    if not isinstance(bins,Iterable):
        nbins = int(bins)
        high = hscale.max()
        low = hscale.min()
        bins = np.linspace(low,high,nbins)

    probs = []
    for low, high in zip(bins,bins[1:]):
        acertos_bin = (acertos[(hscale >=low) & (hscale < high)])
        nbin = len(acertos_bin)
        probs.append(((high+low)/2,nbin,acertos_bin.sum(),acertos_bin.mean(),acertos_bin.std()/np.sqrt(nbin)))
    return np.array(probs)

def invlogit(x):
    return 1/(1+np.e**(-x))

def icclogisticfit(acertos,qn,hscale):
    ''
    X = sm.add_constant(hscale,prepend=True)
    q = acertos[:,qn]
    result = sm.GLM(q,X,family=sm.families.Binomial()).fit()
    const = result.params[0]
    nota  = result.params[1]
    snota = result.bse[0]
    sconst = result.bse[1]
    itemd = -1.*const/nota
    sitemd = itemd*np.sqrt((sconst/const)**2+(snota/nota)**2)
    
    return result,const,sconst,nota,snota,itemd,sitemd

def stats(acertos,hs = None,df = None):
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
    itemstats['id50'] = itemdiscrimination(acertos,frac=0.5)
    itemstats['id25'] = itemdiscrimination(acertos,frac=0.25)

    if hs == "nota":
        hsscale = df['nota']
    elif hs == "scores":
        hscale = df['ressum']
    elif hs == "notapadrao":
        notamean = df['nota'].mean()
        notastd = df['nota'].std() 
        hscale = (df['nota'] - notamean)/notastd
    else:
        return itemstats, teststats 


    probv = []
    iccfits = []
    iccfitsparam =[]
    for qn in range(k):
        probv.append(icc(acertos,qn,hscale,bins=20))
        iccfitresult,const,sconst,nota,snota,itemd,sitemd = icclogisticfit(acertos,qn,hscale)
        iccfits.append(iccfitresult)
        iccfitsparam.append((const,sconst,nota,snota,itemd,sitemd))
    itemstats['icc'] = probv           
    itemstats['iccfit'] = iccfits
    itemstats['iccfitsparam'] = iccfitsparam
    itemstats['hscale'] = hscale

    return itemstats, teststats



if __name__ == '__main__':
    dados = '~/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.txt'
    #dados = '~/enem/Microdados ENEM 2009/Dados Enem 2009/DADOS_ENEM_2009.txt'
    dicfile = '~/enem/Microdados ENEM 2010/Input_SAS/INPUT_SAS_ENEM_2010.SAS'
    #dicfile = '~/enem/Microdados ENEM 2009/Input_SAS/INPUT_SAS_ENEM_2009.sas'
    filtercols = ['NU_INSCRICAO','IDADE','TP_SEXO','TP_COR_RACA','COD_MUNIC_INSC','UF_INSC','IN_TP_ENSINO','IN_PRESENCA_CN','IN_PRESENCA_CH','IN_PRESENCA_LC','IN_PRESENCA_MT','ID_PROVA_CN','NU_NT_CN','TX_RESPOSTAS_CN','DS_GABARITO_CN','ID_PROVA_CH','NU_NT_CH','TX_RESPOSTAS_CH','DS_GABARITO_CH','ID_PROVA_LC','NU_NT_LC','TX_RESPOSTAS_LC','DS_GABARITO_LC','ID_PROVA_MT','NU_NT_MT','TX_RESPOSTAS_MT','DS_GABARITO_MT']
    dic = sasinput(dicfile,filtercols=filtercols)
    out = dados[:-4] + '.csv'
    convert_fff(dados,out,dic,sample = 0.01)

