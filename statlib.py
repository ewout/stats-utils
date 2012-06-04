#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ewout ter Haar <ewout@usp.br>
# License: Apache 

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
    return fsize / linesize


def convert_fff(infile,outfile,dic,max=sys.maxint,sample=1):
    'sample = fraction of lines to sample'
    infile = os.path.expanduser(infile)
    outfile = os.path.expanduser(outfile)
    
    outfile = open(outfile,'w')
    colnames = [col for col,cuts in dic]
    outfile.write('\t'.join(colnames)+'\n')

    flines = filelen(infile)
    if sample < 1:
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



if __name__ == '__main__':
    dados = '~/enem/Microdados ENEM 2010/Dados Enem 2010/DADOS_ENEM_2010.txt'
    dicfile = '~/enem/Microdados ENEM 2010/Input_SAS/INPUT_SAS_ENEM_2010.SAS'
    dic = sasinput(dicfile,filtercols=['NU_INSCRICAO','TX_RESPOSTAS_CN','DS_GABARITO_CN'])
    #dic = sasinput(dicfile)
    print dic
    out = dados[:-3] + 'csv'
    convert_fff(dados,out,dic,sample = 0.001)
