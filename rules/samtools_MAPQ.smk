rule samtools_MAPQ:
  input:
    '/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{sample}.realigned.bam'
  output:
    '/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{sample}.realigned.bam.MAPQ.gz'
  log: 'log/{sample}.realigned.bam.MAPQ.log'
  threads: 12
  message: """ Extract MAPQ per sample using samtools view """
  shell:
    """
    samtools view {input} | cut -f3,4,5 | gzip > {output} 2> {log}
    """
