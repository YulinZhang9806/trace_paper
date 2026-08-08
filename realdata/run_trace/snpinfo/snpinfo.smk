#!python3

from utils import SNPINFO
import numpy as np


rule snpinfo:
    """Get SNP info for all biallelic snps in 1000 Genomes."""
    input:
        tkg_bcf=paths["1000g_bcf"] + "1000g_hg38_chr{chrom}.bcf",
        ancestral=paths["human_ancestor"] + "chr{chrom}.fa",
    params:
        outpref="1000g_biall_snpinfo.human_ancestor.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.chr{chrom}.txt",
    run:
        SNPINFO().get_snp_info(input.tkg_bcf, input.ancestral, params.outpref)


rule snpinfo_archaic:
    """Get SNP info for archaics given SNPs."""
    input:
        txtfile="1000g_biall_snpinfo.human_ancestor.chr{chrom}.txt",
        archaic=paths["archaic_bcf_pref"] + ".chr{chrom}.bcf",
    params:
        snpinfo="1000g_biall_snpinfo.human_ancestor.chr{chrom}",
        archaicbcf=paths["archaic_bcf_pref"] + ".{chrom}",
        outpref="1000g_biall_snpinfo.human_ancestor.archaic.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.chr{chrom}.txt",
    run:
        SNPINFO().append_archaic_snpinfo(
            params.snpinfo, params.archaicbcf, params.outpref
        )

rule snpinfo_afr:
    """Get SNP info for African DAF given SNPs."""
    input:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.chr{chrom}.txt",
        tkg_bcf = paths["1000g_bcf"] + "1000g_hg38_chr{chrom}.bcf",
        yri_poplabel = paths["YRI_sample"],
    params:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.chr{chrom}",
        outpref = "1000g_biall_snpinfo.human_ancestor.archaic.afr.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.afr.chr{chrom}.txt",
    run:
        SNPINFO().append_AFR_info(params.snpinfo, input.tkg_bcf, input.yri_poplabel, "AltAF_YRI", params.outpref)

rule snpinfo_outgroup:
    """Get SNP info for whether they present in hmmix outgroup SNP list."""
    input:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.chr{chrom}.txt",
        outgroupfile = paths["hmmix_outgroup"],
    params:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.chr{chrom}",
        outpref = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.chr{chrom}.txt",
    run:
        SNPINFO().append_outgroup_info(params.snpinfo, input.outgroupfile, params.outpref)


rule snpinfo_strictmask:
    """Get SNP info for if in strictmask."""
    input:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.chr{chrom}.txt",
        strictmask = paths["strictmask"],
        bcffile = paths["1000g_bcf"] + "1000g_hg38_chr{chrom}.bcf",
    params:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.chr{chrom}",
        outpref = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}.txt",
    run:
        SNPINFO().append_strictmask_info(
            params.snpinfo, input.strictmask, input.bcffile, params.outpref
        )

rule snpinfo_manifesto:
    """Get SNP info for if in manifesto filter."""
    input:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}.txt",
        manifesto = paths["manifesto_bed"] + "chr{chrom}.bed",
        bcffile = paths["1000g_bcf"] + "1000g_hg38_chr{chrom}.bcf",
    params:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}",
        outpref = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.manifesto.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.manifesto.chr{chrom}.txt",
    run:
        SNPINFO().append_manifesto_info(params.snpinfo, input.manifesto, input.bcffile, params.outpref)

rule extract_mutage:
    """Extract mutation ages for mutations from SINGER trees"""
    input:
        trees = expand(
            paths["singer_trees"] + "singer_chr{chrom}_{PPID}.tsz", PPID = range(250, 300), allow_missing=True
        ),
    output:
        mutage = expand(
            paths["singer_trees"] + "mutage_chr{chrom}_{PPID}.txt", PPID = range(250, 300), allow_missing=True
        )
    params:
        chrom = lambda wildcards: wildcards.chrom,
    run:
        for i in range(len(input.trees)):
            try:
                ts = tszip.decompress(input.trees[i])
            except Exception as e:
                print(f"Error decompressing {input.trees[i]}: {e}")
                sys.exit(1)
            chrom = params.chrom
            out = "chromosome\tposition\tmutation_age\n"
            for tree in ts.trees():
                for mut in tree.mutations():
                    if tree.parent(mut.node) != tskit.NULL:
                        out += f"{chrom}\t{int(ts.site(mut.site).position)}\t{tree.time(mut.node)}_{tree.time(tree.parent(mut.node))}\n"
                    else:
                        out += f"{chrom}\t{int(ts.site(mut.site).position)}\t{tree.time(mut.node)}_{tree.time(mut.node)}\n"
            with open(output.mutage[i], 'w') as f:
                f.write(out)

rule snpinfo_mutage:
    """Get SNP info for mutage."""
    input:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}.txt",
        mutage = expand(
            paths["singer_trees"] + "mutage_chr{chrom}_{PPID}.txt", PPID = range(250, 300), allow_missing=True
        ),
    params:
        snpinfo = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.chr{chrom}",
        mutage_pref = paths["singer_trees"] + "mutage_chr{chrom}_",
        mutage_range = range(250, 300),
        outpref = "1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.mutage.chr{chrom}",
    output:
        txtfile="1000g_biall_snpinfo.human_ancestor.archaic.afr.outgroup.strictmask.mutage.chr{chrom}.txt",
    run:
        SNPINFO().append_mutage_info(
            params.snpinfo, params.mutage_pref, params.mutage_range, params.outpref
        )
