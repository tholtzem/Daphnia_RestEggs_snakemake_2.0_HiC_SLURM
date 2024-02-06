
#rule newline2csv:
#  input:
#    pops = 'list/pop_list/{pops}.txt'
#  output:
#    pops = 'list/pop_list/{pops}.csv'
#  log: 'log/{pops}_newline2csv.log'
#  threads: 2
#  message: """ Convert newline character separated list to comma separated list """
#  shell:
#    """
#    cat {input.pops} | cut -f1 -d'.' | cut -f2 -d'/' | sed -z 's/\n/,/g;s/,$/\n/' > {output.pops} 2> {log}
#    """

rule get_fixedSites:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20.vcf',
    pop1 = 'ancestry/LC/list/{pop1}.csv',
    pop2 = 'ancestry/LC/list/{pop2}.csv',
    hybrids = 'ancestry/LC/list/hybrids2.csv'
    #hybrids = 'ancestry/LC/list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_withoutRefs.csv'
  output:
    fixed_sites = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites.txt',
    report = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites.report.tsv'
  log:
    'log/get_fixedSites_{pop1}_{pop2}.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of at least 100% in parental species) --- """
  threads: 1 
  shell:
    """
    pop1=$(cat {input.pop1})
    pop2=$(cat {input.pop2})
    hy=$(cat {input.hybrids})
    ruby scripts/get_fixed_site_gts.rb {input.vcf} {output.fixed_sites} $pop1 $hy $pop2 1.0 > {output.report} 2> {log}
    """ 


rule plot_fixedSites:
  input:
    fixed_sites = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites.txt'
  output:
    svg = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites_100percent_thinned1000bp.svg',
    report = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites_100percent_thinned1000bp.report.tsv'
  log:
    'log/plot_fixedSites_{pop1}_{pop2}.fixed_sites_100percent_thinned1000bp.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of 80-100% in all individuals) --- """
  threads: 12
  shell:
    """
    ruby scripts/plot_fixed_site_gts.rb {input.fixed_sites} {output.svg} 1.0 1000 > {output.report} 2> {log}
    """


rule list_fixedSites:
  input:
    fixed_sites = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites.txt'
  output:
    fixed_sites = 'ancestry/LC/all/{pop1}_{pop2}.fixed_sites_SNP.tsv'
  log: 'log/{pop1}_{pop2}.fixed_sites_SNP.log'
  message: """ Get chrom and pos info of fixed sites """
  shell:
    """
    cat {input} | cut -f1,2 | tail -n+2 > {output}
    """


rule vcf_ancestry_sites:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20.vcf',
    sites = 'ancestry/LC/all/ancestry/LC/all/LONGoverall_GALoverall.fixed_sites_SNP.tsv'
  output:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_ancestrySites.vcf.gz'
  log: 'log/ancestrySites.vcf.log'
  message: """ Remove sites with more than 20 % missingness, output uncompressed vcf for downstream analysis """
  shell:
    """
    bcftools view -T {input.sites} -Oz -o {output.vcf} {input.vcf} 2> {log}
    """


rule vcf_select_samples:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_ancestrySites.vcf.gz',
    samples = 'ancestry/LC/list/DW_{DW}_ID.list'
  output:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites.vcf.gz'
  log: 'log/DW_{DW}_ancestrySites.vcf.log'
  message: """ Subsamples vcf for selected samples only   """
  shell:
    """
    bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf} && bcftools index -t {output.vcf} 2> {log}
    """


rule vcf_update_tags:
  input:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites.vcf.gz'
  output:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites_new.vcf.gz'
  log: 'log/DW_{DW}_ancestrySites_new.vcf.log'
  message: """ update tags """
  shell:
    """
    bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t AN,AC,AF,MAF &&
    bcftools index -t {output.vcf} 2> {log}
    """


rule extract_AF:
  input:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites_new.vcf.gz'
  output:
    AF = 'ancestry/LC/data/imissRM/DW_{DW}_AF.tsv'
  log: 'log/DW_{DW}_AF.log'
  message: """ Extract allele frequencies """
  shell:
    """
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' {input.vcf} > {output.AF} 2> {log}
    """




rule index_fixedSites:
  input:
    'ancestry/fixed_sites.list'
  output:
    touch('ancestry/index_fixed_sites.done')
  log: 'log/index_fixed_sites.log'
  threads: 12
  message:
    """ Index list of fixed sites """
  shell:
    """
    module load angsd/0.938
    /apps/uibk/bin/sysconfcpus -n {threads} angsd sites index {input} 2> {log}
   i """



#rule saf_fixedSites:
#  input:
#    ref = config["ref_rapid"],
#    bam = 'realigned/{sample}.realigned.bam',
#    fixed_sites = 'ancestry/fixed_sites.list',
#    touched = 'ancestry/index_fixed_sites.done'
#  output:
#    touch('ancestry/{sample}.GL2.saf.idx.done')
#  log:
#    'log/saf_fixedSites_{sample}.GL2.saf.idx.log'
#  threads: 12
#  message:
#    """ Compute site allele frequency likelihood (.saf.idx) for single samples using angsd (required for the estimation of heterozygosity from single samples and fixed sites)"""
#  shell:
#    """
#    module load angsd/0.938
#    /apps/uibk/bin/sysconfcpus -n 12 angsd -i {input.bam} -ref {input.ref} -anc {input.ref} -out ancestry/{wildcards.sample}.GL2 -doSaf 1 -GL 2 -sites {input.fixed_sites} 2> {log}
#    i"""
#
#
#rule realSFS_samples_fixedSites:
#  input:
#    touched1 = 'ancestry/{sample}.GL2.saf.idx.done',
#    touched2 = 'ancestry/index_fixed_sites.done'
#  output:
#    'ancestry/{sample}.GL2.est.ml'
#  log:
#    'log/saf_fixedSites_{sample}.GL2.est.ml.log'
#  threads: 12
#  message:
#    """ Optimize .saf.idx and estimate the site frequency spectrum (SFS) using realSFS and fold (required for the global heterozygosity estimate for single samples) """
#  shell:
#    """
#    module load singularity/2.x
#    module load angsd/0.938
#    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS ancestry/{wildcards.sample}.GL2.saf.idx -fold 1 > {output} 2> {log}
#    """



