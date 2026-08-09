"""Utility functions for TRACE benchmarking and data processing."""
import os
import numpy as np
import pandas as pd
import tskit
import sys
import pybedtools

class SUMMARIZE:
    """A class of functions to summarize the results."""

    def append_hmmix_info(self, hmmixfile, summaryfile, outpref, inference = "hmmix", individualID = None):
        """
        Append the HMMIX info to the summary file.
        """
        try:
            hmmix = pd.read_csv(hmmixfile, sep="\s+")
            summary = pd.read_csv(summaryfile, sep="\s+")
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)
        if inference == "hmmix":
            hmmix['end'] = hmmix['end'] + 1000
            hmmix = hmmix[hmmix["state"] == "Archaic"]
            hmmix['archaic'] = "Ambiguous"
            hmmix.loc[(hmmix[['AltaiNeandertal', 'Vindija33.19', 'Chagyrskaya-Phalanx']].max(axis=1) > hmmix['Denisova']) & (hmmix['state'] == 'Archaic'), 'archaic'] = 'Neanderthal'
            hmmix.loc[(hmmix[['AltaiNeandertal', 'Vindija33.19', 'Chagyrskaya-Phalanx']].max(axis=1) < hmmix['Denisova']) & (hmmix['state'] == 'Archaic'), 'archaic'] = 'Denisova'
            hmmix['archaic'] = hmmix["mean_prob"].astype(str) + "_" + hmmix['archaic']
            out = ("\t").join(summary.columns) + "\thmmix_start\thmmix_end\thmmix_overlap_length(bp)\tmean_pp\thmmix_assign\n"
        elif inference == "ibdmix":
            if individualID is None:
                print("Please provide the individual ID for ibdmix inference.")
                sys.exit(1)
            hmmix = hmmix[(hmmix["ID"] == individualID) & (hmmix["end"] > hmmix["start"])]
            hmmix['archaic'] = (hmmix['end'] - hmmix['start']).astype(str) + "_" + hmmix['archaic'].astype(str)
            out = ("\t").join(summary.columns) + "\tibdmix_start\tibdmix_end\tibdmix_overlap_length(bp)\tslod\tibdmix_assign\n"
        elif inference == "archie":
            hmmix = pd.read_csv(hmmixfile, sep="\s+", header=None, names=["chrom", "start", "end", "pp", "ID", "overlap"])
            if individualID is None:
                print("Please provide the individual ID for ibdmix inference.")
                sys.exit(1)
            hmmix = hmmix[hmmix["ID"] == individualID]
            hmmix = hmmix[(hmmix['pp'] > 0.9) & (hmmix['overlap'] > 0.7)]
            hmmix['archaic'] = hmmix['pp'].astype(str) + "_ghost"
            out = ("\t").join(summary.columns) + "\tarchie_start\tarchie_end\tarchie_overlap_length(bp)\tpp\tarchie_assign\n"
        s_segs = pybedtools.BedTool(summary[["chromosome", "start", "end"]].to_csv(sep="\t", index=False, header=False), from_string=True)
        h_segs = pybedtools.BedTool(hmmix[["chrom", "start", "end", "archaic"]].to_csv(sep="\t", index=False, header=False), from_string=True)
        overlap = s_segs.intersect(h_segs, wao=True)
        ovl = {}
        for i in range(len(overlap)):
            ol = int(overlap[i][-1])
            if ol > 0:
                key = f"{overlap[i][0]}_{overlap[i][1]}_{overlap[i][2]}"
                if key not in ovl:
                    ovl[key] = []
                ovl[key].append([f"{overlap[i][4]}_{overlap[i][5]}", ol, overlap[i][6]])
        for s in range(len(summary)):
            chrom = summary['chromosome'][s]
            start = summary['start'][s]
            end = summary['end'][s]
            out += ("\t").join(summary.iloc[s].astype(str).tolist())
            if f"{chrom}_{start}_{end}" in ovl:
                ss = ovl[f"{chrom}_{start}_{end}"]
                if len(ss) > 1:
                    hmm_start = ""
                    hmm_end = ""
                    hmm_ol = 0
                    hmm_pp = []
                    hmm_assign = ""
                    for i in range(len(ss)):
                        ss_s1 = ss[i][0].split("_")
                        hmm_start += f"{ss_s1[0]},"
                        hmm_end += f"{ss_s1[1]},"
                        hmm_ol += ss[i][1]
                        ss_s2 = ss[i][2].split("_")
                        hmm_pp.append(float(ss_s2[0]))
                        hmm_assign += f"{ss_s2[1]},"
                    out += f"\t{hmm_start[:-1]}\t{hmm_end[:-1]}\t{hmm_ol}\t{np.max(hmm_pp)}\t{hmm_assign[:-1]}\n"
                else:
                    ss = ss[0]
                    ss_s1 = ss[0].split("_")
                    ss_s2 = ss[2].split("_")
                    out += f"\t{ss_s1[0]}\t{ss_s1[1]}\t{ss[1]}\t{ss_s2[0]}\t{ss_s2[1]}\n"
            else:
                out += "\t-1\t-1\t-1\t-1\t-1\n"
        with open(f"{outpref}.txt", "w") as f:
            f.write(out)
        return

    def parse_snpinfo_anc(self, snpinfo):
        """Read in snpinfo file and return a dictionary with ancestral information.

        snpinfo: a str of full file name

        return: a dictionary {pos:[following information in the file]}
        """
        infile=open(str(snpinfo))
        lines=infile.readlines()
        infile.close()
        outdict = dict({})
        col_dict = {}
        cols = lines[0].strip('\n').strip('\t').split('\t')
        for j in range(len(cols)):
            col_dict[cols[j]] = j
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            pos = int(s[col_dict['pos']])
            ref = str(s[col_dict['ref']])
            alt = str(s[col_dict['alt']])
            anc = str(s[col_dict['ancestral']])
            cpg = True if int(s[col_dict['CpG']]) == 1 else False
            dnea = ((s[col_dict['Chagyrskaya-Phalanx']] in ['1','2']) | (s[col_dict['AltaiNeandertal']] in ['1','2']) | (s[col_dict['Vindija33.19']] in ['1','2']))
            # dnea = (s[col_dict['Vindija33.19']] in ['1','2'])
            nea_missing = (s[col_dict['Chagyrskaya-Phalanx']] == '9') & (s[col_dict['AltaiNeandertal']] == '9') & (s[col_dict['Vindija33.19']] == '9')
            # nea_missing = (s[col_dict['Vindija33.19']] == '9')
            dden = (s[col_dict['Denisova']] in ['1','2'])
            den_missing = (s[col_dict['Denisova']] == '9')
            freqs = float(s[col_dict['AltAF_YRI']])
            aout = False if int(s[col_dict['in_outgroup']]) == 1 else True
            instrict = True if int(s[col_dict['in_strictmask']]) == 1 else False
            if ref == anc:
                outdict["_".join([str(pos), ref, alt])] = ["keep", dnea, nea_missing, dden, den_missing, aout, instrict, freqs, cpg]
            elif alt == anc:
                outdict["_".join([str(pos), ref, alt])] = ["switch", dnea, nea_missing, dden, den_missing, aout, instrict, freqs, cpg]
            elif anc in ['A', 'T', 'C', 'G']:
                outdict["_".join([str(pos), ref, alt])] = ["other", dnea, nea_missing, dden, den_missing, aout, instrict, freqs, cpg]
            else:
                outdict["_".join([str(pos), ref, alt])] = ["N", dnea, nea_missing, dden, den_missing, aout, instrict, freqs, cpg]
        return outdict

    def count_snps(self, seg_snps, snp_dict, hap="left"):
        """Count number of derived and archaic snps in the region file.

        seg_snps: str, full name of the region file
        snp_dict: dictionary including snp information

        return: nsnps, nout, ND00, ND10, ND01, ND11, A05, AN01, AD01
        """
        snpcount = {"nsnps":0, "nout":0, "ND00":0, "ND10":0, "ND01":0, "ND11":0}
        snpcount_strict = {"nsnps":0, "nout":0, "ND00":0, "ND10":0, "ND01":0, "ND11":0}
        dsnps = []
        dsnps_freqs = []
        dsnps_marks = []
        hapdict = {"left":0, "right":1}
        infile=open(seg_snps)
        lines=infile.readlines()
        infile.close()
        for i in range(len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            pos = "_".join([s[0], s[1], s[2]])
            if not pos in snp_dict:
                print(pos)
                print("Error: SNP not found in snp_dict.")
                sys.exit(1)
            if str(s[3]) in ['./.', './1', './0', '0/.', '1/.']: # skip missing sites
                continue
            genos = str(s[3]).split('|')
            geno = int(genos[hapdict[hap]])
            ancinfo, dnea, nea_missing, dden, den_missing, aout, instrict, freqs, cpg = snp_dict[pos]
            if ancinfo == 'keep':
                if geno > 0:
                    snpcount['nsnps'] += 1
                    snpcount['nout'] += (aout)
                    mark = "other"
                    if not dnea and not nea_missing and not dden and not den_missing:
                        snpcount['ND00'] += 1
                        mark = "ND00"
                    elif dnea and not dden and not den_missing:
                        snpcount['ND10'] += 1
                        mark = "ND10"
                    elif not dnea and not nea_missing and dden:
                        snpcount['ND01'] += 1
                        mark = "ND01"
                    elif dnea and dden:
                        snpcount['ND11'] += 1
                        mark = "ND11"
                    if not cpg:
                        dsnps.append(pos)
                        dsnps_freqs.append(freqs)
                    if instrict:
                        snpcount_strict['nsnps'] += 1
                        snpcount_strict['nout'] += (aout)
                        if not dnea and not nea_missing and not dden and not den_missing:
                            snpcount_strict['ND00'] += 1
                            mark = "ND00_strict"
                        elif dnea and not dden and not den_missing:
                            snpcount_strict['ND10'] += 1
                            mark = "ND10_strict"
                        elif not dnea and not nea_missing and dden:
                            snpcount_strict['ND01'] += 1
                            mark = "ND01_strict"
                        elif dnea and dden:
                            snpcount_strict['ND11'] += 1
                            mark = "ND11_strict"
                    if not cpg:
                        dsnps_marks.append(mark)
            elif ancinfo == 'switch':
                if (1 - geno) > 0:
                    snpcount['nsnps'] += 1
                    snpcount['nout'] += (aout)
                    mark = "other"
                    if dnea and dden:
                        snpcount['ND00'] += 1
                        mark = "ND00"
                    elif not dnea and not nea_missing and dden:
                        snpcount['ND10'] += 1
                        mark = "ND10"
                    elif dnea and not dden and not den_missing:
                        snpcount['ND01'] += 1
                        mark = "ND01"
                    elif not dnea and not nea_missing and not dden and not den_missing:
                        snpcount['ND11'] += 1
                        mark = "ND11"
                    if not cpg:
                        dsnps.append(pos)
                        dsnps_freqs.append(1 - freqs)
                    if instrict:
                        snpcount_strict['nsnps'] += 1
                        snpcount_strict['nout'] += (aout)
                        if dnea and dden:
                            snpcount_strict['ND00'] += 1
                            mark = "ND00_strict"
                        elif not dnea and not nea_missing and dden:
                            snpcount_strict['ND10'] += 1
                            mark = "ND10_strict"
                        elif dnea and not dden and not den_missing:
                            snpcount_strict['ND01'] += 1
                            mark = "ND01_strict"
                        elif not dnea and not nea_missing and not dden and not den_missing:
                            snpcount_strict['ND11'] += 1
                            mark = "ND11_strict"
                    if not cpg:
                        dsnps_marks.append(mark)
            elif ancinfo == 'other':
                snpcount['nsnps'] += 1
                snpcount['nout'] += (aout)
                if not cpg:
                    dsnps.append(pos)
                    freqs = 0*freqs
                    dsnps_freqs.append(freqs)
                    dsnps_marks.append("other")
                if geno > 0:
                    snpcount['ND00'] += (not dnea and not nea_missing and not dden and not den_missing)
                    snpcount['ND10'] += (dnea and not dden and not den_missing)
                    snpcount['ND01'] += (not dnea and not nea_missing and dden)
                    snpcount['ND11'] += (dnea and dden)
                    if instrict:
                        snpcount_strict['nsnps'] += 1
                        snpcount_strict['nout'] += (aout)
                        snpcount_strict['ND00'] += (not dnea and not nea_missing and not dden and not den_missing)
                        snpcount_strict['ND10'] += (dnea and not dden and not den_missing)
                        snpcount_strict['ND01'] += (not dnea and not nea_missing and dden)
                        snpcount_strict['ND11'] += (dnea and dden)
                else:
                    snpcount['ND00'] += (dnea and dden)
                    snpcount['ND10'] += (not dnea and not nea_missing and dden)
                    snpcount['ND01'] += (dnea and not dden and not den_missing)
                    snpcount['ND11'] += (not dnea and not nea_missing and not dden and not den_missing)
                    if instrict:
                        snpcount_strict['nsnps'] += 1
                        snpcount_strict['nout'] += (aout)
                        snpcount_strict['ND00'] += (dnea and dden)
                        snpcount_strict['ND10'] += (not dnea and not nea_missing and dden)
                        snpcount_strict['ND01'] += (dnea and not dden and not den_missing)
                        snpcount_strict['ND11'] += (not dnea and not nea_missing and not dden and not den_missing)
        return snpcount, snpcount_strict, np.array(dsnps), np.array(dsnps_freqs), np.array(dsnps_marks)

    def get_regions_bed(self, file):
        """Helper function to read a bed file, concatenate records to bcftools regions string."""
        out = ""
        infile=open(file)
        lines=infile.readlines()
        infile.close() 
        for i in range(len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            out += s[0] + ':' + s[1] + '-' + s[2]
            if not i == len(lines) - 1:
                out += ','
        return out

    def final_ind_count(self, samplename, summary, hap, snpinfo, bcfpref, outpref):
        """Create summary table for number of mutations per archaic haplotype for each individual.

        samplename: str
        bed: str, bed file prefix of pure archaic regions for the haplotype.
        hap: left or right
        snpinfo: str, snpinfo file prefix, should not include chr identifier.
        bcfpref: str, prefix of the original bcffile, should not include chr identifier.

        output: a txt file with archaic segment information and counts of mutations.
        """         
        if os.path.exists(str(outpref) + ".txt"):
            os.remove(str(outpref) + ".txt")
        infile=open(summary)
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tnderived\tnoutgroup\tND00\tND10\tND01\tND11\tnderived_strict\tnoutgroup_strict\t'
        out += 'ND00_strict\tND10_strict\tND01_strict\tND11_strict\n'
        out1 = lines[0].strip('\n') + "\tdsnps\tdsnps_marks\tDAF_YRI\n"
        cur_chr = lines[1].strip('\n').split('\t')[0]
        snpfile = str(snpinfo) + '.' + cur_chr + '.txt'
        bcffile = str(bcfpref) + cur_chr + '.bcf'
        snp_dict = self.parse_snpinfo_anc(snpfile)
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            if not cur_chr == str(s[0]):
                if not os.path.exists(str(bcfpref) + str(s[0]) + '.bcf'):
                    print(f"File {bcfpref}{s[0]}.bcf does not exist, exit...")
                    sys.exit(1)
                else:
                    cur_chr = str(s[0])
                    snpfile = str(snpinfo) + '.' + cur_chr + '.txt'
                    bcffile = str(bcfpref) + cur_chr + '.bcf'
                    snp_dict = self.parse_snpinfo_anc(snpfile)
            os.system('echo "' + str(s[0]) + "\t" + str(s[1]) + "\t" + str(s[2]) + '" > ' + str(outpref) + str(i) + 'seg.bed')
            reg = self.get_regions_bed(str(outpref) + str(i) + "seg.bed")
            os.system(
                "bcftools view -s " + str(samplename) + " -v snps -r " + str(reg) + " " + str(bcffile) + " | bcftools query -f'[%POS\t%REF\t%ALT\t%GT\n]' > "+ str(outpref) + str(i) + "seg_snps"
            )
            snpcount, snpcount_strict, dsnps, dsnps_freqs, dsnps_marks = self.count_snps(str(outpref) + str(i) + "seg_snps", snp_dict, hap)
            out += lines[i].strip('\n') + '\t' + str(snpcount['nsnps']) + '\t' + str(snpcount['nout']) + '\t' + str(snpcount['ND00']) + '\t' 
            out += str(snpcount['ND10']) + '\t' + str(snpcount['ND01']) + '\t' + str(snpcount['ND11']) + '\t' 
            out += str(snpcount_strict['nsnps']) + '\t' + str(snpcount_strict['nout']) + '\t' 
            out += str(snpcount_strict['ND00']) + '\t' + str(snpcount_strict['ND10']) + '\t' + str(snpcount_strict['ND01']) + '\t' + str(snpcount_strict['ND11']) + "\n"
            dsnps_freqs = np.round(dsnps_freqs, 3)
            if len(dsnps) > 0:
                out1 += lines[i].strip('\n') + '\t' + ','.join(dsnps) + '\t'
                out1 += ','.join(dsnps_marks) + '\t'
                out1 += ','.join([str(i) for i in dsnps_freqs]) + '\n'
            else:
                out1 += lines[i].strip('\n') + '\t' + 'NA\tNA\tNA\n'
            os.remove(str(outpref) + str(i) + 'seg.bed')
            os.remove(str(outpref) + str(i) + "seg_snps")
        outfile=open(outpref + '.txt','a')
        outfile.write(out)
        outfile.close()
        outfile=open(outpref + '_DAF.txt','w')
        outfile.write(out1)
        outfile.close()

    def append_t1_t2(self, npzpref, summary, snpinfo, func, outpref):
        mutage = {}
        infile = open(snpinfo)
        lines = infile.readlines()
        infile.close()
        col_dict = {}
        cols = lines[0].strip('\n').strip('\t').split('\t')
        for j in range(len(cols)):
            col_dict[cols[j]] = j
        for i in range(1, len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            pos = int(s[col_dict["pos"]])
            if not s[-1] in ['nan', 'Not_mapped', 'NA']:
                ag = float(s[col_dict["mutage"]])
                mutage[pos] = ag
        with np.load(f"{npzpref}") as data:
            t1 = data["t1s"]
            t2 = data["t2s"]
            treespan_phy = data["treespan_phy"]
        windowsize = treespan_phy[0][1] - treespan_phy[0][0]
        infile=open(summary)
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tt1s\tt2s\tmutages\tbranch_mark\n'
        col_dict = {}
        cols = lines[0].strip('\n').strip('\t').split('\t')
        for j in range(len(cols)):
            col_dict[cols[j]] = j
        for i in range(1, len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            st = int(int(s[1]) / windowsize)
            ed = int(int(s[2]) / windowsize)
            t1_val = ",".join(t1[st:ed].astype('str'))
            t2_val = ",".join(t2[st:ed].astype('str'))
            out += lines[i].strip('\n') + '\t' + str(t1_val) + '\t' + str(t2_val) + '\t'
            muts = s[col_dict['dsnps']].split(',')
            mks = []
            ags = []
            for m in muts:
                pos = int(m.split('_')[0])
                if pos in mutage:
                    ags.append(mutage[pos])
                    if mutage[pos] > t1[int(pos / windowsize)] and mutage[pos] < t2[int(pos / windowsize)]:
                        mks.append("on")
                    elif mutage[pos] <= t1[int(pos / windowsize)]:
                        mks.append("below")
                    elif mutage[pos] >= t2[int(pos / windowsize)]:
                        mks.append("above")
                else:
                    ags.append("NA")
                    mks.append("NA")
            out += "\t" + ",".join([str(i) for i in ags]) + "\t" + ",".join(mks) + '\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()

    def append_seg_freq(self, npzpref, indlist, state_arc, outpref, chrom):
        """Append segment frequency in defined population to the summary file."""
        if state_arc == "NEA":
            sn = "nea_states"
        elif state_arc == "DEN":
            sn = "den_states"
        elif state_arc == "Ghost":
            sn = "ghost_states"
        else:
            sys.exit(f"Unrecognized state_arc {state_arc}")
        with np.load(f"{npzpref}{indlist[0]}.{chrom}.xss.npz") as data:
            treespan_phy = data["treespan_phy"]
        states = np.zeros((len(indlist), treespan_phy.shape[0]))
        for idx in range(len(indlist)):
            with np.load(f"{npzpref}{indlist[idx]}.{chrom}.xss.npz") as data:
                try:
                    s = data[sn]
                    states[idx] = s
                except:
                    print(f"Error: {indlist[idx]} does not have {sn} in {npzpref}{indlist[idx]}.{chrom}.xss.npz")
                    continue
        np.savez_compressed(f"{outpref}.npz", states=states, treespan_phy=treespan_phy)

    def append_bscore_recombrate(self, popcountfile, bscore, recombrate):
        """Append Bscore and recombination rate for the region."""
        infile = open(popcountfile, 'r')
        lines = infile.readlines()
        infile.close()
        b2 = pybedtools.BedTool(bscore)
        b3 = pybedtools.BedTool(recombrate)
        out = lines[0].strip('\n') + '\tbscore\trecombrate\n'
        for i in range(1, len(lines)):
            s = lines[i].strip('\n').split('\t')
            bedstring = str(s[0]) + '\t' + str(s[1]) + '\t' + str(s[2])
            b1 = pybedtools.BedTool(bedstring, from_string=True)
            ib1 = b2.intersect(b1)
            bsc = np.sum([int(ib1[i][3]) * (int(ib1[i][2]) - int(ib1[i][1])) for i in range(len(ib1))]) / (int(s[2]) - int(s[1]))
            out += lines[i].strip('\n') + '\t' + str(bsc) + '\t'
            ib1 = b3.intersect(b1)
            rec = np.sum([float(ib1[i][3]) * (int(ib1[i][2]) - int(ib1[i][1])) for i in range(len(ib1))]) / (int(s[2]) - int(s[1]))
            out += str(float(rec)) + '\n'
        outfile = open(popcountfile, 'w')
        outfile.write(out)
        outfile.close()

    def bscore_recombrate_windows(self, npzfile, bscore, recombrate, chrom):
        """Get per 1000bp Bscore and recombination rate."""
        d = np.load(npzfile)
        data = {k: d[k] for k in d.files} 
        treespan_phy = data["treespan_phy"]
        treespan_cM = TRACE().add_recombination_map(treespan_phy, recombrate)
        recrate = (treespan_cM[:, 1] - treespan_cM[:, 0]) / (1000 * 1e-6)
        bstring = ""
        for i in range(treespan_phy.shape[0]):
            bstring += f"{chrom}\t{treespan_phy[i, 0]}\t{treespan_phy[i, 1]}\n"
        b1 = pybedtools.BedTool(bstring, from_string=True).sort()
        b2 = pybedtools.BedTool(bscore).sort()
        maped = b1.map(b2, c=4, o='mean', null=0)
        bsc = np.array([float(field[3]) for field in maped])
        data["bscore"] = bsc
        data["recombrate"] = recrate
        np.savez_compressed(npzfile, **data)