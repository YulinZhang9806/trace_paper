import numpy as np
import sys
from utils import SUMMARIZE
import pandas as pd
import os
import argparse

def merge_df_count_branch_mut(shared, count, count_t1t2, output, classify_group = 1):
    data = pd.read_csv(shared, sep="\t")
    len_cols = len(data.columns)
    data = pd.read_csv(count_t1t2, sep="\s+")
    dff = pd.read_csv(count, sep="\s+")
    data = data.merge(dff, on=dff.columns[0:len_cols].tolist(), how='left')
    nd00_aff = []
    nd10_aff = []
    nd01_aff = []
    nd11_aff = []
    yri_aff = []
    tot_aff = []
    for i in data.index:
        nd00 = 0
        nd10 = 0
        nd01 = 0
        nd11 = 0
        tot = 0
        nyri = 0
        mks = data.loc[i, 'dsnps_marks'].split(',')
        ons = data.loc[i, 'branch_mark'].split(',')
        yri = str(data.loc[i, 'DAF_YRI']).split(',')
        assert len(mks) == len(ons) == len(yri)
        for j in range(len(mks)):
            if ons[j] == 'on' and mks[j] in ['ND00_strict', 'ND10_strict', 'ND01_strict', 'ND11_strict']:
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
    data['nd00_prop'] = nd00_aff / tot_aff
    data['nd10_prop'] = nd10_aff / (nd00_aff + nd10_aff + nd01_aff)
    data['nd01_prop'] = nd01_aff / (nd00_aff + nd10_aff + nd01_aff)
    data['nYRI'] = yri_aff / tot_aff
    data['assign_label'] = "NONE"
    if classify_group == 1:
        data.loc[(data['nd10_prop'] > 0.2) & (data['nd10_b'] > data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "NEA"
        data.loc[(data['nd01_prop'] > 0.2) & (data['nd10_b'] < data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "DEN"
        data.loc[(data['nd00_prop'] > 0.8) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "Ghost"
    elif classify_group == 2:
        data.loc[((data['nd00_prop'] <= 0.8) | (data['nYRI'] <= 0.1)) & (data['nd10_b'] > data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "NEA"
        data.loc[((data['nd00_prop'] <= 0.8) | (data['nYRI'] <= 0.1)) & (data['nd10_b'] < data['nd01_b']) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "DEN"
        data.loc[(data['nd00_prop'] > 0.8) & (data['nYRI'] > 0.1) & (data['nderived_strict'] >= 30) & (data[["ND00_strict", "ND10_strict", "ND01_strict", "ND11_strict"]].sum(axis=1) >= 10), 'assign_label'] = "Ghost"
    else:
        raise ValueError("Invalid classify_group value. Must be 1 or 2.")
    data.to_csv(output, sep="\t", index=False)

def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description="Annotate TRACE detected archaic segments.")
    
    # Add arguments with prefixes
    parser.add_argument("--ind", required=True, help="individual ID")
    parser.add_argument("--inpref", required=True, help="input summary file prefix (parts before .summary.txt)")
    parser.add_argument("--chrom", required=True, help="chromosome, must match the chromosome field in input summary file")
    parser.add_argument("--sample-info", required=True, help="sample information file")
    parser.add_argument("--xss-pref", required=True, help="prefix for xss files (parts before .[chrom].xss.npz)")
    parser.add_argument("--hmmix-path", help='path to hmmix decoded outputs, assume file names are [samplename].hap[1/2].txt;' 
    + 'if not provided, will skip appending hmmix info and only use ibdmix for annotation',
    default=None
    )
    parser.add_argument("--ibdmix-path", required=True, help="path to ibdmix outputs, assume file names are [chrom].txt")
    parser.add_argument("--snpinfo-pref", required=True, help="prefix to snpinfo file with mutation age info appended (parts before .[chrom].txt)")
    parser.add_argument("--bcf-pref", required=True, help="prefix to bcf files for extracting allele frequency info, assume file names are [prefix][chrom].bcf")
    parser.add_argument("--classify-group", type=int, choices=[1, 2], default=1, help="classification rule group for final annotation, 1 for more stringent classification and 2 for more relaxed classification")
    parser.add_argument("--outpref", required=True, help="output prefix for the final annotated tables, file names will be [outpref].count_merge.txt")

    # add parse args for these
    args = parser.parse_args()
    ind = int(args.ind)
    inpref = args.inpref
    chrom = args.chrom
    sample_info = args.sample_info
    hmmixpath = args.hmmix_path
    ibdmixpath = args.ibdmix_path
    xsspref = args.xss_pref
    outpref = args.outpref
    snpinfo = args.snpinfo_pref
    bcfpref = args.bcf_pref
    classify_group = args.classify_group

    df = pd.read_csv(sample_info, sep="\t")
    inds = [ind]
    hmmixpath = hmmixpath if hmmixpath is not None and hmmixpath.endswith("/") else (hmmixpath + "/" if hmmixpath is not None else None)
    hmmixfiles = [f"{hmmixpath}{df[df['ID']==inds[i]]['Name'].values[0]}.txt" for i in range(len(inds))]
    ibdmixfile = f"{ibdmixpath}{chrom}.txt" if ibdmixpath.endswith("/") else f"{ibdmixpath}/{chrom}.txt"
    individualID = [df[df['ID']==inds[i]]['SampleID'].values[0] for i in range(len(inds))]
    xssfile = f"{xsspref}.{chrom}.xss.npz"

    for i in range(len(inds)):
        prefix = inpref + ".summary"
        if hmmixpath is not None and os.path.exists(hmmixfiles[i]):
            prefixn = outpref + ".hmmix"
            SUMMARIZE().append_hmmix_info(hmmixfiles[i], prefix + ".txt", prefixn, inference = "hmmix", individualID = individualID[i])
            prefix = prefixn
            prefixn = prefix + ".ibdmix"
            SUMMARIZE().append_hmmix_info(ibdmixfile, prefix + ".txt", prefixn, inference = "ibdmix", individualID = individualID[i])
            os.remove(prefix + ".txt")
        else:
            print(f"HMMix file for individual {inds[i]} not found at {hmmixfiles[i]}. Skipping HMMix annotation.")
            prefixn = outpref + ".ibdmix"
            SUMMARIZE().append_hmmix_info(ibdmixfile, prefix + ".txt", prefixn, inference = "ibdmix", individualID = individualID[i])

    for i in range(len(inds)):
        hap = "left" if inds[i]%2 == 0 else "right"
        if hmmixpath is not None and os.path.exists(hmmixfiles[i]):
            SUMMARIZE().final_ind_count(
                samplename = individualID[i], 
                summary = outpref + ".hmmix.ibdmix.txt", 
                hap = hap, 
                snpinfo = snpinfo, 
                bcfpref = bcfpref, 
                outpref = outpref + ".count"
            )
            SUMMARIZE().append_t1_t2(
                npzpref = xssfile,
                summary = outpref + ".count_DAF.txt", 
                snpinfo = snpinfo + f".{chrom}.txt", 
                func = "median",
                outpref = outpref + ".count_t1t2",
            )
        else:
            SUMMARIZE().final_ind_count(
                samplename = individualID[i], 
                summary = outpref + ".ibdmix.txt", 
                hap = hap, 
                snpinfo = snpinfo, 
                bcfpref = bcfpref, 
                outpref = outpref + ".count"
            )
            SUMMARIZE().append_t1_t2(
                npzpref = xssfile,
                summary = outpref + ".count_DAF.txt", 
                snpinfo = snpinfo + f".{chrom}.txt",
                func = "median",
                outpref = outpref + ".count_t1t2",
            )
        os.remove(outpref + ".count_DAF.txt")
        merge_df_count_branch_mut(
            shared = outpref + ".ibdmix.txt" if hmmixpath is None or not os.path.exists(hmmixfiles[i]) else outpref + ".hmmix.ibdmix.txt",
            count = outpref + ".count.txt", 
            count_t1t2 = outpref + ".count_t1t2.txt", 
            output = outpref + ".count_merge.txt",
            classify_group = classify_group
        )
        os.remove(outpref + ".count.txt")
        os.remove(outpref + ".count_t1t2.txt")
        with np.load(xssfile) as d:
            data = {k: d[k] for k in d.files} 
        states = data["states"]
        treespan_phy = data["treespan_phy"]
        window_size = int(treespan_phy[0, 1] - treespan_phy[0, 0])
        nea_states = np.zeros(shape=states.shape, dtype=states.dtype)
        den_states = np.zeros(shape=states.shape, dtype=states.dtype)
        ghost_states = np.zeros(shape=states.shape, dtype=states.dtype)
        df = pd.read_csv(outpref + ".count_merge.txt", sep="\t")
        for j in range(len(df)):
            if df.loc[j, 'assign_label'] == "NEA":
                nea_states[int(df.loc[j, 'start'] / window_size):int(df.loc[j, 'end'] / window_size)] = 1
            elif df.loc[j, 'assign_label'] == "DEN":
                den_states[int(df.loc[j, 'start'] / window_size):int(df.loc[j, 'end'] / window_size)] = 1
            elif df.loc[j, 'assign_label'] == "Ghost":
                ghost_states[int(df.loc[j, 'start'] / window_size):int(df.loc[j, 'end'] / window_size)] = 1
        data["nea_states"] = nea_states
        data["den_states"] = den_states
        data["ghost_states"] = ghost_states
        np.savez_compressed(
            xssfile, **data
        )

if __name__ == "__main__":
    main()
