rule faidx_ref:
  input:
    ref = config["ref_HiC"],
  output:
    'ref/ref_HiC_chromsize.tsv'
  log: 'log/new_{pop1}_{pop2}.fixed_sites_{pct}percent_SNP.log'
  message: """ Get chrom and pos info of fixed sites """
  shell:
    """
    cat {ref} | cut -f1,2 > {output}
    """
