rule angsd_ngsRelate_POP:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list'
  output:
    touch('ngsRelate/{POPS}.glf3.done')
  log:
   'log/ngsRelate/{POPS}.glf3.log'
  threads: 2
  resources: mem_mb=100000, walltime="24:00:00"
  message:
    """ Estimate allele frequencies and calculate GLs while SNP calling for each population using angsd """
  shell:
    """
    angsd -b {input.bamlist} -ref {input.ref} -out ngsRelate/{wildcards.POPS}.glf3 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -skipTriallelic 1 -GL 2 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -doCounts 1 2> {log}
    """


rule extract_freq:
  input:
    'ngsRelate/{POPS}.glf3.done'
  output:
    'ngsRelate/{POPS}.glf3.freq'
  log:
   'log/ngsRelate/{POPS}.glf3.freq.log'
  threads: 1
  resources: mem_mb=100, walltime="00:05:00"
  message:
    """ Extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)  """
  shell:
    """
    zcat ngsRelate/{wildcards.POPS}.glf3.mafs.gz | cut -f6 | sed 1d > {output} 2> {log}
    """


#rule ngsRelate:
#  input:
#    touched = 'ngsRelate/{POPS}.glf3.done',
#    freq = 'ngsRelate/{POPS}.glf3.freq',
#    bamlist = 'list/saf/{POPS}.list'
#  output:
#    'ngsRelate/{POPS}.glf3.stats.txt'
#  log: 'ngsRelate/{POPS}.glf3.stats.log'
#  threads: 2
#  resources: mem_mb=80000, walltime="24:00:00"
#  message:
#    """--- Calculating relatedness (etc.) with ngsRelate2 from genotype likelihoods (glf3) ---"""
#  shell:
#    """
#    N=`cat {input.bamlist} | wc -l`
#    ngsrelate -g angsdput.glf.gz -p {threads} -n $N -f {input.freq} -O {output} 2> {log}
#    """



