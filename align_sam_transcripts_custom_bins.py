import numpy as np
import os,sys,argparse
import pandas as pd
import time
import glob
import matplotlib as plt

def alignSamGenic(SamFile, utr3Folder, utr5Folder, cdsFolder, output, sora):
    df100 = pd.read_csv(SamFile, delim_whitespace=True, names=['Qname', 'flag', 'transcript_id', 'start', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'lol', 'a7a'], usecols=['Qname', 'transcript_id', 'start', 'SEQ'])
    df100['length'] = df100['SEQ'].str.len()
    df100['end'] = df100['start'] + df100['length']
    dfs=[]
    all_files = glob.glob(utr3Folder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    utr3 = frame[['transcript_id', 'total_length']]
    utr3['transcript_id'] = utr3['transcript_id'].str.replace('_3utr', '')
    utr3.columns = ['transcript_id', 'utr3_length']
    dfs=[]
    all_files = glob.glob(utr5Folder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    utr5 = frame[['transcript_id', 'total_length']]
    utr5['transcript_id'] = utr5['transcript_id'].str.replace('_5utr', '')
    utr5.columns = ['transcript_id', 'utr5_length']
    dfs=[]
    all_files = glob.glob(cdsFolder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0)
    frame1 = frame.sort_values(by='chr')
    frame = frame1.reset_index()
    cds = frame[['transcript_id', 'total_length']]
    cds['transcript_id'] = cds['transcript_id'].str.replace('_cds', '')
    cds.columns = ['transcript_id', 'cds_length']
    df1 = df100.merge(utr5, on=['transcript_id'], how='left')
    df2 = df1.merge(cds, on=['transcript_id'], how='left')
    df3 = df2.merge(utr3, on=['transcript_id'], how='left')
    df3 = df3[df3['cds_length'].notna()]
    print("number of reads: " + str(len(df3)) + SamFile)
    n = 20
    for i in range(1, 21):
        df3[i] = (df3['utr5_length']/20) * i
    for i in range(1,101):
        df3[i + 20] = ((df3['cds_length']/100) * i) + df3['utr5_length']
    for i in range(1,71):
        df3[i + 120] = ((df3['utr3_length']/70) * i) + df3['utr5_length'] + df3['cds_length']
    x = 190
    names = list(range(1, x + 1))
    names += ['transcript_id']
    names += ['start']
    names += ['end']
    names += ['length']
    df4 = pd.DataFrame(df3, columns=names)    
    lol = pd.DataFrame(columns=names)
    df4['1b'] = np.where((df4['start'] < df4[1]) & (df4[1] > df4['end']), 1, 0)
    df4['1c'] = np.where((df4['start'] < df4[1]) & (df4[1] < df4['end']), ((df4[1] - df4['start'])/df4['length']), 0)
    df4['1_final'] = df4['1b'] + df4['1c']
    for i in range(2,x+1):
        df4[str(i) + 'b'] = np.where((df4['start'] > df4[i - 1]) & (df4['start'] < df4[i]) & (df4[i] > df4['end']), 1, 0)
        df4[str(i) + 'c'] = np.where((df4['start'] > df4[i - 1]) & (df4['start'] < df4[i]) & (df4[i] < df4['end']), ((df4[i] - df4['start'])/df4['length']), 0)
        df4[str(i) + 'd'] = np.where((df4['start'] < df4[i - 1]) & (df4[i] > df4['end']) & (df4['end'] > df4[i - 1]), ((df4['end'] - df4[i - 1])/df4['length']), 0)
        df4[str(i) + 'e'] = np.where((df4['start'] < df4[i - 1]) & (df4['end'] > df4[i]), (df4[i] - df4[i - 1])/df4['length'], 0)
        df4[str(i) + '_final'] = df4[str(i) + 'b'] + df4[str(i) + 'c'] + df4[str(i) + 'e'] + df4[str(i) + 'd']

    names1 = list(range(1, x + 1))
    mystring = '_final'
    names = ["{}{}".format(i,mystring) for i in names1]
    names += ['transcript_id']
    names += ['start']
    names += ['end']
    names += ['length']
    lol = pd.DataFrame(df4, columns=names)
    lol1 = lol.fillna(0)
    del lol1['start']
    del lol1['end']
    del lol1['length']
    lol1 = lol1.set_index('transcript_id')
    lol2 = lol1.T
    lol2['sum'] = lol2.sum(axis=1)
    lol2['normalized'] = lol2['sum']/lol2['sum'].sum()
    lol3 = lol2.reset_index()    
    lol3.to_csv(output, index=False, sep='\t')
    ax = lol3['normalized'].head(20).plot(style='.-', figsize=(20,10), color='red')
    plt.pyplot.xticks(fontsize=22)
    plt.pyplot.yticks(fontsize=22)
    plt.pyplot.xlabel('\nbins', fontsize=24)
    plt.pyplot.ylabel('reads', fontsize = 24)
    fig = ax.get_figure()
    ax = lol3['normalized'].iloc[20:120].plot(style='.-', figsize=(20,10), color='orange')
    plt.pyplot.xticks(fontsize=22)
    plt.pyplot.yticks(fontsize=22)
    plt.pyplot.xlabel('\nbins', fontsize=24)
    plt.pyplot.ylabel('reads', fontsize = 24)
    fig = ax.get_figure()
    ax = lol3['normalized'].tail(70).plot(style='.-', figsize=(20,10), color='blue')
    plt.pyplot.xticks(fontsize=22)
    plt.pyplot.yticks(fontsize=22)
    plt.pyplot.xlabel('\nBinned exonic length', fontsize=24)
    plt.pyplot.ylabel('Reads density', fontsize = 24)
    ax.set_ylim(0,0.02)
    fig = ax.get_figure()

    fig.savefig(sora, dpi=100)
def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-SamFile', required=True)
    parser.add_argument('-utr3Folder', required=True)
    parser.add_argument('-utr5Folder', required=True)
    parser.add_argument('-cdsFolder', required=True)
    parser.add_argument('-output', required=True)
    parser.add_argument('-sora', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    filename = alignSamGenic(args.SamFile, args.utr3Folder, args.utr5Folder, args.cdsFolder, args.output, args.sora)
    end = time.time()
    print ('time elapsed:' + str(end - start))
