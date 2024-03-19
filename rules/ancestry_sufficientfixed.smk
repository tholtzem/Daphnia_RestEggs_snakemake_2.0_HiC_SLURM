
rule get_sufficient_fixedSites:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20.vcf',
    pop1 = 'ancestry/LC/list/{pop1}.csv',
    pop2 = 'ancestry/LC/list/{pop2}.csv',
    hybrids = 'ancestry/LC/list/hybrids2.csv'
    #hybrids = 'ancestry/LC/list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_withoutRefs.csv'
  output:
    fixed_sites = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}.txt',
    report = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}.report.tsv'
  log:
    'log/new_get_fixedSites_{pct}percent_sufficient{pct}_{pop1}_{pop2}.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of at least 100% in parental species) --- """
  threads: 1 
  shell:
    """
    prop=$(echo "scale=1; {wildcards.pct} / 100" | bc)
    pop1=$(cat {input.pop1})
    pop2=$(cat {input.pop2})
    hy=$(cat {input.hybrids})
    ruby scripts/get_fixed_site_gts.rb {input.vcf} {output.fixed_sites} $pop1 $hy $pop2 $prop $prop > {output.report} 2> {log}
    """ 



rule plot_fixedSites:
  input:
    fixed_sites = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}.txt'
  output:
    svg = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_thinned1000bp.svg',
    report = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_thinned1000bp.report.tsv'
  log:
    'log/plot_fixedSites_{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_thinned1000bp.log'
  message: """ --- Run ruby script (https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/get_fixed_site_gts.rhttps://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_fixed_site_gts.rbb) to get the alleles at sites that are fixed differently in two parental species (complete fixation of 80-100% in all individuals) --- """
  threads: 1
  shell:
    """
    prop=$(echo "scale=1; {wildcards.pct} / 100" | bc)
    ruby scripts/plot_fixed_site_gts.rb {input.fixed_sites} {output.svg} $prop 1000 > {output.report} 2> {log}
    """


rule list_fixedSites:
  input:
    fixed_sites = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}.txt'
  output:
    fixed_sites = 'ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_SNP.tsv'
  log: 'log/new_{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_SNP.log'
  message: """ Get chrom and pos info of fixed sites """
  shell:
    """
    cat {input} | cut -f1,2 | tail -n+2 > {output}
    """


rule vcf_ancestry_sites:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20.vcf',
    sites = 'ancestry/LC/overall/LONGoverall0.99_GALoverall0.99.fixed_sites_{pct}percent_sufficient{pct}_SNP.tsv'
  output:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites{pct}.vcf.gz'
  log: 'log/ancestrySites{pct}_MajorMinor4.vcf.log'
  message: """ Subset vcf """
  shell:
    """
    bcftools view -T {input.sites} -Oz -o {output.vcf} {input.vcf} 2> {log}
    """

rule vcf_select_samples:
  input:
    vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites{pct}.vcf.gz',
    samples = 'ancestry/LC/list/DW_{DW}_ID.list'
  output:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4.vcf.gz'
  log: 'log/DW_{DW}_ancestrySites{pct}_MajorMinor4.vcf.log'
  message: """ Subsamples vcf for selected samples only   """
  shell:
    """
    bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf} && bcftools index -t {output.vcf} 2> {log}
    """


rule vcf_update_tags:
  input:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4.vcf.gz'
  output:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4_new.vcf.gz'
  log: 'log/DW_{DW}_ancestrySites{pct}_MajorMinor4_new.vcf.log'
  message: """ update tags """
  shell:
    """
    bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t AN,AC,AF,MAF &&
    bcftools index -t {output.vcf} 2> {log}
    """


rule extract_AF:
  input:
    vcf = 'ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4_new.vcf.gz'
  output:
    AF = 'ancestry/LC/data/imissRM/DW_{DW}_AF_ancestrySites{pct}_MajorMinor4_new.tsv'
  log: 'log/DW_{DW}_AF_ancestrySites{pct}_MajorMinor4_new.log'
  message: """ Extract allele frequencies """
  shell:
    """
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' {input.vcf} > {output.AF} 2> {log}
    """


