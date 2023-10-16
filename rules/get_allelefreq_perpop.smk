
rule get_mafs_POP_randomSNP:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list',
    sites = 'list/mafs/PRE_long_LC_vs_ALL_gal_vs_PEL_long_LC_1000random_SNPs_sorted.tsv',
    chrom = 'list/mafs/PRE_long_LC_vs_ALL_gal_vs_PEL_long_LC_1000random_SNPs_sorted.chr'
  output:
    touch('mafs/{POPS}_1000random_SNPs.mafs.done')
  log:
   'log/mafs/{POPS}_1000random_SNPs.mafs.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    angsd sites index {input.sites} &&
    angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out mafs/{wildcards.POPS}_1000random_SNPs -doSaf 1 -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -sites {input.sites} -rf {input.chrom} 2> {log}
    """

rule get_mafs_POP_randomSite:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list',
    sites = 'list/mafs/PRE_long_LC_vs_ALL_gal_vs_PEL_long_LC_1000random_sites_sorted.tsv',
    chrom = 'list/mafs/PRE_long_LC_vs_ALL_gal_vs_PEL_long_LC_1000random_sites_sorted.chr'
  output:
    touch('mafs/{POPS}_1000random_sites.mafs.done')
  log:
   'log/mafs/{POPS}_1000random_sites.mafs.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    angsd sites index {input.sites} &&
    angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out mafs/{wildcards.POPS}_1000random_sites -doSaf 1 -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -sites {input.sites} -rf {input.chrom} 2> {log}
i    """


rule get_mafs_POP_randomoutliers:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list',
    sites = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_1000random_fst_outliers_sorted.tsv',
    chrom = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_1000random_fst_outliers_sorted.chr'
  output:
    touch('mafs/{POPS}_1000random_fst_outliers.mafs.done')
  log:
   'log/mafs/{POPS}_1000random_fst_outliers.mafs.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    angsd sites index {input.sites} &&
    angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out mafs/{wildcards.POPS}_1000random_fst_outliers -doSaf 1 -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -sites {input.sites} -rf {input.chrom} 2> {log}
    """

rule get_mafs_POP_top10NOTfixed:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list',
    sites = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_top10_NOTfixed_fst_outliers_sorted.tsv',
    chrom = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_top10_NOTfixed_fst_outliers_sorted.chr'
  output:
    touch('mafs/{POPS}_top10_NOTfixed_fst_outliers.mafs.done')
  log:
   'log/mafs/{POPS}_top10_NOTfixed_fst_outliers.mafs.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    angsd sites index {input.sites} &&
    angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out mafs/{wildcards.POPS}_top10_NOTfixed_fst_outliers -doSaf 1 -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -sites {input.sites} -rf {input.chrom} 2> {log}
    """


rule get_mafs_POP_random10fixed:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list',
    #sites = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_10random_fixed_fst_outliers_sorted.tsv',
    #chrom = 'list/mafs/PRE_long_LC_vs_PEL_long_LC_10random_fixed_fst_outliers_sorted.chr'
    sites = '../Daphnia_RestEggs_snakemake_pbs_2.0_HiC/ngsLD/LC/LDpruned_snps_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_0.5Minweight.list',
    chrom = '../Daphnia_RestEggs_snakemake_pbs_2.0_HiC/ngsLD/LC/LDpruned_snps_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_0.5Minweight.chr'
  output:
    touch('mafs/{POPS}_10random_fixed_fst_outliers.mafs.done')
  log:
   'log/mafs/{POPS}_10random_fixed_fst_outliers.mafs.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    angsd -b {input.bamlist} -ref {input.ref} -out mafs/{wildcards.POPS}_10random_fixed_fst_outliers -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doBcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 8 -sites {input.sites} -rf {input.chrom} 2> {log}
    """
