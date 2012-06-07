import pandas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter


def resvec(df,rescol,gabcol):
    ''
    res = df[rescol]
    gab = df[gabcol]
    l = []
    for rvec,gvec in zip(res,gab):
        l.append([1 if x==y else 0 for x,y in zip(rvec,gvec)])

    a = np.array(l)
    return a
        
def resvec2(df,rescol,gabcol):
    ''
    res = df[rescol]
    gab = df[gabcol]
    l = []
    for rvec in res:
        l.append(list(rvec))
    a = np.array(l)
    for qn,x in enumerate(a.T):
        df['Q'+str(qn+1)] = x
    return df


def gradebar(grade,correct):
    ''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = Counter(grade)
    labels = sorted(c.keys())
    values = [c[val] if not val == correct for val in labels]
    N = len(labels)
    x = np.arange(N)
    width = 0.6
    ax.bar(x,values,width,color='r')
    ax.set_xticks(x+0.5*width)
    ax.set_xticklabels(labels)
    return fig


if __name__ == '__main__':

    df = pandas.read_table('DADOS_ENEM_2010.csv')
    df = df[(df['ID_PROVA_CN'] == 89) & (df['NU_NT_CN'] != '         ')]
    df['nota'] = df['NU_NT_CN'].apply(float)

