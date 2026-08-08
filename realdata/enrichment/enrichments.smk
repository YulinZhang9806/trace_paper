#!python3

import numpy as np
import pandas as pd


configfile: "config.yaml"


rule all:
    input:
        expand(
            "results/enrichments/{dataset}/great/{p}.hg38.cell_component.csv",
            p=["NEA", "DEN", "Ghost"],
            dataset=config["enrichments"]["datasets"],
        ),


rule annotate_manifesto:
    """Annotate the input tables with SNPs in archaic manifesto regions."""
    input:
        filtered_segments=lambda wildcards: config["enrichments"]["datasets"][
            wildcards.population
        ],
        manifesto_snpinfo=expand(
            config["paths"]["manifesto_snpinfo_pref"] + ".chr{CHROM}.txt",
            CHROM=range(1, 23),
        ),
    output:
        manifesto_table="results/enrichments/{population}/manifesto_filtered_segments.tsv",
    params:
        snpinfo_pref=config["paths"]["manifesto_snpinfo_pref"],
    shell:
        """
python scripts/annotate_manifesto.py {input.filtered_segments} {params.snpinfo_pref} {output.manifesto_table}

"""


rule determine_peaks_bed:
    """Create a bed file for peaks of archaic ancestry."""
    input:
        filtered_segments="results/enrichments/{population}/manifesto_filtered_segments.tsv",
        chrom_sizes=config["paths"]["chrom_sizes"],
        bed_mask=config["paths"]["strict_mask"],
    output:
        bed=temp("results/enrichments/{population}/bed/{target_pop}.hg38.peaks.bed"),
    wildcard_constraints:
        population="AFR|EUR|OCN|EAS|SAS",
        target_pop="NEA|DEN|Ghost",
    params:
        target_pop=lambda wildcards: f"{wildcards.target_pop}",
        n_derived_strict=30,
        nd_filter=10,
        sigma=2,
    script:
        "scripts/peak_identification.py"


rule detect_enrichment_depletion:
    """Detect regions with enrichment or depletion of specific ancestries."""
    input:
        filtered_segments=lambda wildcards: config["enrichments"]["datasets"][
            wildcards.population
        ],
        chrom_sizes=config["paths"]["chrom_sizes"],
        bed_mask=config["paths"]["strict_mask"],
    output:
        raw_enrichments="results/enrichments/{population}/bed/{target_pop}.hg38.enrichments.raw.tsv",
        raw_depletions="results/enrichments/{population}/bed/{target_pop}.hg38.depletion.raw.tsv",
        enrichments="results/enrichments/{population}/bed/{target_pop}.hg38.enrichments.tsv",
        depletions="results/enrichments/{population}/bed/{target_pop}.hg38.depletion.tsv",
    wildcard_constraints:
        population="AFR|EUR|OCN|EAS|SAS",
        target_pop="NEA|DEN|Ghost",
    params:
        target_pop=lambda wildcards: f"{wildcards.target_pop}",
        n_derived_strict=30,
        nd_filter=10,
        quantile=0.01,
    script:
        "scripts/region_enrichments.py"


rule detect_deserts:
    input:
        filtered_segments=lambda wildcards: config["enrichments"]["datasets"][
            wildcards.population
        ],
        chrom_sizes=config["paths"]["chrom_sizes"],
        bed_mask=config["paths"]["strict_mask"],
    output:
        deserts="results/deserts/{population}/bed/{target_pop}.hg38.deserts.bed",
    wildcard_constraints:
        population="AFR|EUR|OCN|EAS|SAS",
        target_pop="NEA|DEN|Ghost",
    params:
        target_pop=lambda wildcards: f"{wildcards.target_pop}",
        n_derived_strict=30,
        nd_filter=10,
        min_length=10e6,
        max_freq=0.001,
    script:
        "scripts/detect_deserts.py"


rule filter_deserts_strict_mask:
    """Ensure that detected deserts are sufficiently within the strict mask."""
    input:
        bed=rules.detect_deserts.output.deserts,
        bed_mask=config["paths"]["strict_mask"],
        centromere=config["paths"]["centromeres"],
    output:
        raw_bed="results/deserts/{population}/bed/{target_pop}.hg38.deserts.overlaps.bed",
        bed="results/deserts/{population}/bed/{target_pop}.hg38.deserts.filt.bed",
    params:
        centromere_overlap=0.001,
        strict_overlap=0.8,
    shell:
        """
bedtools intersect -a {input.bed} -b {input.centromere} {input.bed_mask} -names centromere strict -wao | awk \'{{a[$1\"-\"$2\"-\"$3\"-\"$4\"-\"$5\"-\"$6\"-\"$7\"-\"$8] += $NF}} END {{for (i in a){{print i, a[i]}}}}\' | sed \'s/-/\t/g\' | awk \'{{OFS=\"\t\"; print $0,$NF/($3-$2)}}\' > {output.raw_bed}
awk \'($8 == \"strict\") && ($9 > {params.strict_overlap}) || (($8 == \"centromere\") && ($9 < {params.centromere_overlap}))\' {output.raw_bed} | awk \'{{OFS=\"\t\";print $1,$2,$3,$4,$5,$6,$7}}\' | sort | uniq | bedtools sort > {output.bed}

"""


rule sort_peaks_file:
    input:
        bed=rules.determine_peaks_bed.output.bed,
    output:
        bed="results/enrichments/{population}/bed/{target_pop}.hg38.peaks.sorted.bed",
    shell:
        "bedtools sort -i {input.bed} > {output.bed}"


rule run_enrichment_broad_greatR:
    """Run enrichment testing using greatR for the dataset."""
    input:
        bed=rules.sort_peaks_file.output.bed,
        enrichments=rules.detect_enrichment_depletion.output.enrichments,
        depletions=rules.detect_enrichment_depletion.output.depletions,
        deserts=rules.filter_deserts_strict_mask.output.bed,
    output:
        go_molecular_tsv="results/enrichments/{population}/great/{target_pop}.hg38.go_molecular.csv",
        go_bioprocess_tsv="results/enrichments/{population}/great/{target_pop}.hg38.bioprocess.csv",
        go_cellcomponent_tsv="results/enrichments/{population}/great/{target_pop}.hg38.cell_component.csv",
    params:
        outfix=lambda wildcards: f"results/enrichments/{wildcards.population}/great/{wildcards.target_pop}",
    shell:
        "Rscript scripts/enrichment_great.R {input.bed} {params.outfix}"
