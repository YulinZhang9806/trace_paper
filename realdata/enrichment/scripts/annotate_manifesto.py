import numpy as np
import sys
import pandas as pd
import os

def parse_snpinfo_anc(snpinfo):
    infile=open(str(snpinfo))
    lines=infile.readlines()
    infile.close()
    outdict = dict({})
    for i in range(1, len(lines)):
        s=lines[i].strip('\n').strip('\t').split('\t')
        pos = int(s[1])
        ref = str(s[2])
        alt = str(s[3])
        instrict = True if ((int(s[-2]) == 1) & (int(s[-1]) == 1)) else False
        outdict["_".join([str(pos), ref, alt])] = instrict
    return outdict

def merge_df_count_branch_mut(countmerge, snpinfo, unrec_label = "Ghost"):
    data = pd.read_csv(countmerge, sep="\t")
    data = data.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    cur_chr = data.loc[0, 'chromosome']
    print("Processing chromosome:", cur_chr)
    snpinfo_file = str(snpinfo) + '.' + cur_chr + '.txt'
    snpinfo_dict = parse_snpinfo_anc(snpinfo_file)
    nd00_aff = []
    nd10_aff = []
    nd01_aff = []
    nd11_aff = []
    yri_aff = []
    tot_aff = []
    for i in data.index:
        chr = data.loc[i, 'chromosome']
        if chr != cur_chr:
            cur_chr = chr
            print("Processing chromosome:", cur_chr)
            snpinfo_file = str(snpinfo) + '.' + cur_chr + '.txt'
            snpinfo_dict = parse_snpinfo_anc(snpinfo_file)
        nd00 = 0
        nd10 = 0
        nd01 = 0
        nd11 = 0
        tot = 0
        nyri = 0
        snps = data.loc[i, 'dsnps'].split(',')
        mks = data.loc[i, 'dsnps_marks'].split(',')
        ons = data.loc[i, 'branch_mark'].split(',')
        yri = str(data.loc[i, 'DAF_YRI']).split(',')
        new_mks = ""
        for j in range(len(mks)):
            if ons[j] == 'on' and mks[j] in ['ND00_strict', 'ND10_strict', 'ND01_strict', 'ND11_strict'] and snpinfo_dict[snps[j]]:
                tot += 1
                if mks[j] == 'ND00_strict':
                    nd00 += 1
                elif mks[j] == 'ND10_strict':
                    nd10 += 1
                elif mks[j] == 'ND01_strict':
                    nd01 += 1
                elif mks[j] == 'ND11_strict':
                    nd11 += 1
                if float(yri[j]) > 0:
                    nyri += 1
            if mks[j] in ['ND00_strict', 'ND10_strict', 'ND01_strict', 'ND11_strict'] and not snpinfo_dict[snps[j]]:
                if mks[j] == 'ND00_strict':
                    new_mks += 'ND00'
                elif mks[j] == 'ND10_strict':
                    new_mks += 'ND10'
                elif mks[j] == 'ND01_strict':
                    new_mks += 'ND01'
                elif mks[j] == 'ND11_strict':
                    new_mks += 'ND11'
            else:
                new_mks += mks[j]
            new_mks += ','
        data.loc[i, 'dsnps_marks'] = new_mks.strip(',')
        if tot == 0:
            nd00_aff.append(np.nan)
            nd10_aff.append(np.nan)
            nd01_aff.append(np.nan)
            nd11_aff.append(np.nan)
            tot_aff.append(np.nan)
            yri_aff.append(np.nan)
        else:
            nd00_aff.append(nd00)
            nd10_aff.append(nd10)
            nd01_aff.append(nd01)
            nd11_aff.append(nd11)
            tot_aff.append(tot)
            yri_aff.append(nyri)
    nd00_aff = np.array(nd00_aff)
    nd10_aff = np.array(nd10_aff)
    nd01_aff = np.array(nd01_aff)
    nd11_aff = np.array(nd11_aff)
    tot_aff = np.array(tot_aff)
    yri_aff = np.array(yri_aff)
    data['nd00_b'] = nd00_aff
    data['nd10_b'] = nd10_aff
    data['nd01_b'] = nd01_aff
    data['nd11_b'] = nd11_aff
    data['tot_b'] = tot_aff
    data['yri_b'] = yri_aff
    data['nea_prop'] = (nd10_aff + nd11_aff) / tot_aff
    data['den_prop'] = (nd01_aff + nd11_aff) / tot_aff
    data['nd00_prop'] = nd00_aff / tot_aff
    data['nd11_prop'] = nd11_aff / tot_aff
    data['nd10_prop'] = nd10_aff / (nd00_aff + nd10_aff + nd01_aff)
    data['nd01_prop'] = nd01_aff / (nd00_aff + nd10_aff + nd01_aff)
    data['nYRI'] = yri_aff / tot_aff
    data['assign_label'] = "NONE"
    data.loc[(data['nd10_prop'] > 0.2) & (data['nd10_b'] > data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "NEA"
    data.loc[(data['nd01_prop'] > 0.2) & (data['nd10_b'] < data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "DEN"
    data.loc[(data['nd00_prop'] > 0.8) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = unrec_label
    return data

countmerge = sys.argv[1]
snpinfo = sys.argv[2]
output = sys.argv[3]
unrec_label = sys.argv[4] if len(sys.argv) > 4 else "Ghost"

data = merge_df_count_branch_mut(countmerge, snpinfo, unrec_label)
data.to_csv(output, sep="\t", index=False)