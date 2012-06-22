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


def stats(acertos,hscale = None):
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
    itemstats['id27'] = itemdiscrimination(acertos,frac=0.27)
    probv = []
    if hasattr(hscale,'max'):
        for qn in range(k):
            probs = icc(acertos,qn,hscale,bins=20)
            probv.append(probs)
            itemstats['icc'] = probv
    return itemstats, teststats 



if __name__ == '__main__':
    dados = '~/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.txt'
    dicfile = '~/enem/Microdados ENEM 2010/Input_SAS/INPUT_SAS_ENEM_2010.SAS'
    filtercols = ['NU_INSCRICAO','IDADE','TP_SEXO','TP_COR_RACA','ID_PROVA_CN','NU_NT_CN','TX_RESPOSTAS_CN','DS_GABARITO_CN','ID_PROVA_CH','NU_NT_CH','TX_RESPOSTAS_CH','DS_GABARITO_CH']
    dic = sasinput(dicfile,filtercols=filtercols)
    #chdic = sasinput(dicfile,filtercols=['NU_INSCRICAO','ID_PROVA_CH','NU_NT_CH','TX_RESPOSTAS_CH','DS_GABARITO_CH'])
    #outcn = dados[:-4] + '-CN.csv'
    #outch = dados[:-4] + '-CH.csv'
    out = dados[:-4] + '.csv'
    convert_fff(dados,out,dic,sample = 0.01)
    #convert_fff(dados,outch,chdic,sample = 0.01)
