"""Utility functions for SNP info annotation."""
import os
import numpy as np
from pyfaidx import Fasta

class SNPINFO:
    """A class of functions to get SNP ancestral and archaic information."""
    
    ## checked
    def get_snp_info(self, bcffile, ancesterfa, outpref):
        """Get SNP info file with polarization information.
    
        bcffile: a file end with ".bcf" and could be read by bcftools.
        ancesterfa: a fasta file of ancestral genome, end with ".fa"

        return: a txt file with snp information
        """
        if os.path.exists(str(outpref) + ".txt"):
            os.remove(str(outpref) + ".txt")
        # ignore indels
        os.system(
            "bcftools view -v snps " + str(bcffile) + " | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > " + str(outpref) + ".bcftools_temp.txt"
        )
        type_dict = {'AT':'TA', 'TA':'TA', 'TC':'TC', 'AG':'TC', 'TG':'TG', 'AC':'TG', 'CT':'CT', 'GA':'CT', 'CA':'CA', 'GT':'CA', 'CG':'CG', 'GC':'CG', 'NN':'NN'}
        ancestral = Fasta(ancesterfa)
        c = list(ancestral.keys())[0]
        infile=open(str(outpref) + ".bcftools_temp.txt")
        lines=infile.readlines()
        infile.close()
        out = "chr\tpos\tref\talt\tancestral\tderived\tupstream\tdownstream\ttype\ttype_fold\tCpG\n"
        for i in range(len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            pos = int(s[1])
            ref = str(s[2])
            if len(ref) > 1: # skip indels
                continue
            alts = str(s[3]).strip('\t').strip(',').split(',')
            for alt in alts:
                out += "\t".join(lines[i].strip('\n').split('\t')[:-1]) + '\t' + alt + '\t'
                triplit = ancestral[c][pos-2:pos+1].seq
                up = list(triplit)[0]
                anc = list(triplit)[1]
                down = list(triplit)[2]
                out += anc.upper() + '\t'
                dr = 'N'
                if anc.upper() == ref.upper():
                    dr = alt.upper()
                    out += dr + '\t'
                elif anc.upper() == alt.upper():
                    dr = ref.upper()
                    out += dr + '\t'
                elif anc.upper() in ['A', 'T', 'C', 'G']:
                    dr = ref.upper() + ',' + alt.upper()
                    out += dr + '\t'
                else:
                    anc = 'N'
                    out += 'N\t'
                out += up.upper() + '\t' + down.upper() + '\t'
                if len(dr) == 1:
                    ty = anc.upper() + dr.upper()
                else:
                    ty = anc.upper() + dr.split(',')[0] + ',' + anc.upper() + dr.split(',')[1]
                out += ty + '\t'
                try:
                    out += type_dict[ty] + '\t'
                except:
                    out += type_dict[ty.split(',')[0]] + ',' + type_dict[ty.split(',')[1]] + '\t'
                if anc.upper() == 'C' and down.upper() == 'G':
                    out += str(1)
                elif anc.upper() == 'G' and up.upper() == 'C':
                    out += str(1)
                else:
                    out += str(0)
                out += '\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()
        os.remove(str(outpref) + ".bcftools_temp.txt")

    ## checked
    def parse_archaic_geno(self, archaic_geno, arc_return_dict):
        """Helper function to parse archaic_geno file.

        archaic_geno: a file of "bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]'" output
        arc_return_dict: a dictionary indicating the output genotype order

        return: a set and a dictionary (pos:[genotypes(order by arc_return_dict), ref, alt])
        """
        geno_dict = {'0/0':0, '0/1':1, '1/0':1, '1/1':2, './.':9, './0':9, '0/.':9, './1':9, '1/.':9}
        infile=open(archaic_geno)
        lines=infile.readlines()
        infile.close()
        existing_sites = set()
        geno_info = dict()
        for i in range(len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            if len(str(s[3])) > 1: # archaics are multiallelic at the site, ignore
                continue
            if not "_".join([s[1], s[2], s[3]]) in existing_sites:
                existing_sites.add("_".join([s[1], s[2], s[3]]))
                geno_info["_".join([s[1], s[2], s[3]])] = [0] * (len(arc_return_dict) + 2)
                geno_info["_".join([s[1], s[2], s[3]])][-2] = str(s[2])
                geno_info["_".join([s[1], s[2], s[3]])][-1] = str(s[3])
            geno_info["_".join([s[1], s[2], s[3]])][arc_return_dict[str(s[4])]] = geno_dict[str(s[5])]
        return existing_sites, geno_info

    ## checked
    def append_archaic_snpinfo(self, snpinfo, archaicbcf, outpref):
        """Append archaic snp info to the existing snpinfo file.
        
        snpinfo: str, no file extension
        archaicbcf: str, no file extension, need to have multiallelic sites split to multiple records

        return: new snpinfo txt file.
        """
        # chrom, pos from snpinfo file, remove header
        os.system(
            "cut -f 1,2 " + str(snpinfo) + ".txt | tail -n +2 > " + snpinfo + "allsnp"
        )
        # used split multiallelic site versions here
        os.system(
            "bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' -R " + snpinfo + "allsnp " + str(archaicbcf) + ".bcf" + " > " + snpinfo + "archaic_geno"
        )
        os.remove(snpinfo + "allsnp")
        os.system(
            "bcftools query -l " + str(archaicbcf) + ".bcf" + " > " + snpinfo + "archaic_samples"
        )
        with open(snpinfo + "archaic_samples") as f:
            archaic_samples = f.read().strip().split('\n')
        arc_return_dict = {x: idx for idx, x in enumerate(archaic_samples)}
        os.remove(snpinfo + "archaic_samples")
        existing_sites, geno_info = self.parse_archaic_geno(snpinfo + "archaic_geno", arc_return_dict)
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\t' + '\t'.join(archaic_samples) + '\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if "_".join([s[1], s[2], s[3]]) in existing_sites:
                ginfo = geno_info["_".join([s[1], s[2], s[3]])]
            elif "_".join([s[1], s[2], "."]) in existing_sites:
                ginfo = geno_info["_".join([s[1], s[2], "."])]
            else:
                out += "\t".join([str(9)] * len(archaic_samples)) + '\n'
                continue # missing site
            if ginfo[-2].upper() == str(s[2]) and (ginfo[-1].upper() == str(s[3]) or ginfo[-1] == '.'):
                for ars in archaic_samples:
                    out += str(ginfo[arc_return_dict[ars]]) + '\t'
                out = out.strip('\t') + '\n'
            else:
                out += "\t".join(["REF_DONT_MATCH"] * len(archaic_samples)) + '\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()
        os.remove(snpinfo + "archaic_geno")
        return existing_sites, geno_info
    
    ## checked
    def parse_AFR_info(self, bcffile, tempfile, poplabelfile):
        """Helper function to parse AFR SNP info."""
        os.system(
            "bcftools view -v snps -Ov " + str(bcffile) + "| vcftools --vcf - --keep "+str(poplabelfile)+" --freq --stdout > " + str(tempfile)
        )
        infile = open(str(tempfile))
        lines = infile.readlines()
        infile.close()
        afr_dict = dict()
        afr_set = set()
        for i in range(1, len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            nalleles = int(s[2])
            afr_set.add(int(s[1]))
            if not int(s[1]) in afr_dict:
                afr_dict[int(s[1])] = dict()
            for j in range(len(s) - nalleles, len(s)):
                allele, freq = s[j].split(':')
                afr_dict[int(s[1])][allele] = float(freq)
        os.remove(str(tempfile))
        return afr_dict, afr_set

    ## checked
    def append_AFR_info(self, snpinfo, bcffile, poplabelfile, linelabel, outpref):
        """Append AFR info to the existing snpinfo file."""
        afr_dict, afr_set = self.parse_AFR_info(bcffile, outpref + "afr_frq.txt", poplabelfile)
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\t' + str(linelabel) + '\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if int(s[1]) in afr_set:
                out += str(afr_dict[int(s[1])][str(s[3])]) + '\n'
            else:
                out += "missing_site\n"
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()

    def parse_outgroup_hmmix(self, outgroupfile):
        """Parse outgroup info from hmmix outgroup file."""
        infile=open(outgroupfile)
        lines=infile.readlines()
        infile.close()
        outgroupset = set()
        for i in range(1, len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            outgroupset.add("_".join([str(s[0]), str(s[1])]))
        return outgroupset

    def append_outgroup_info(self, snpinfo, outgroupfile, outpref):
        """Append outgroup info to the existing snpinfo file."""
        outgroupset = self.parse_outgroup_hmmix(outgroupfile)
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tin_outgroup\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if "_".join([s[0], s[1]]) in outgroupset:
                out += '1\n'
            else:
                out += '0\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()

    def parse_strictmask(self, strictmask, bcffile, tempfile):
        """Parse mutations in strictmask."""
        os.system(
            "bcftools view -v snps -R " + str(strictmask) + " " + str(bcffile) + " | bcftools query -f '%CHROM\t%POS\n' > " + str(tempfile)
        )
        infile=open(tempfile)
        lines=infile.readlines()
        infile.close()
        strictmaskset = set()
        for i in range(len(lines)):
            s = lines[i].strip('\n').strip('\t').split('\t')
            strictmaskset.add("_".join([str(s[0]), str(s[1])]))
        os.remove(tempfile)
        return strictmaskset

    def append_strictmask_info(self, snpinfo, strictmask, bcffile, outpref):
        """Append strictmask info to the existing snpinfo file."""
        strictmaskset = self.parse_strictmask(strictmask, bcffile, outpref + "strictmask.txt")
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tin_strictmask\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if "_".join([s[0], s[1]]) in strictmaskset:
                out += '1\n'
            else:
                out += '0\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()

    def append_manifesto_info(self, snpinfo, manifesto, bcffile, outpref):
        """Append strictmask info to the existing snpinfo file."""
        strictmaskset = self.parse_strictmask(manifesto, bcffile, outpref + "manifesto.txt")
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tin_manifesto\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if "_".join([s[0], s[1]]) in strictmaskset:
                out += '1\n'
            else:
                out += '0\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()

    def append_mutage_info(self, snpinfo, mutage_pref, mutage_range, outpref):
        """Append mutage info to the existing snpinfo file."""
        mutage = {}
        for i in mutage_range:
            mutage_sub = {}
            infile = open(str(mutage_pref) + str(i) + ".txt")
            lines = infile.readlines()
            infile.close()
            for j in range(1, len(lines)):
                s = lines[j].strip('\n').strip('\t').split('\t')
                ss = s[2].split('_')
                l = float(ss[0])
                h = float(ss[1])
                if int(s[1]) not in mutage_sub: # get uniquely mapped mutations
                    if l == h:
                        mutage_sub[int(s[1])] = np.nan # remove map to root node mutations
                    else:
                        mutage_sub[int(s[1])] = (l + h) / 2
                else:
                    mutage_sub[int(s[1])] = np.nan
            for k, v in mutage_sub.items():
                if not k in mutage:
                    mutage[k] = []
                mutage[k].append(v)
        infile=open(snpinfo + '.txt')
        lines=infile.readlines()
        infile.close()
        out = lines[0].strip('\n') + '\tmutage\n'
        for i in range(1, len(lines)):
            s=lines[i].strip('\n').strip('\t').split('\t')
            out += lines[i].strip('\n') + '\t'
            if int(s[1]) in mutage:
                if len(mutage[int(s[1])]) > 0:
                    out += str(np.nanmean(mutage[int(s[1])])) + '\n'
                else:
                    out += 'NA\n'
            else:
                out += 'Not_mapped\n'
        outfile=open(outpref + '.txt','w')
        outfile.write(out)
        outfile.close()