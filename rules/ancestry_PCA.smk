localrules: getChrom_from_ancestry_sites

rule index_ancestry_sites:
  input:
    'ancestry/LC/all/{prefix}_SNP.tsv'
  output:
    touch('ancestry/LC/all/index_{prefix}_SNP.done')
  log: 'log/ancestry/LC/all/index_{prefix}_SNP.log'
  threads: 1
  resources: mem_mb=500, walltime="00:30:00"
  message:
    """ Index ancestry sites file """
  shell:
    """
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd/angsd sites index {input} 2> {log}
    """


rule getChrom_from_ancestry_sites:
  input:
    'ancestry/LC/all/{prefix}_SNP.tsv'
  output:
     'ancestry/LC/all/{prefix}_SNP.chr'
  log:
    'log/ancestry/LC/all/get_chrom_{prefix}.log'
  message:
    """ Get a list of chromosomes/scaffolds """
  shell:
    """
     cut -f1 {input} | uniq > {output}
    """



rule get_ancestry_GLs:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/LC_realignedBAM_df1.list',
    touched = 'ancestry/LC/all/index_{prefix}_SNP.done',
    sites = 'ancestry/LC/all/{prefix}_SNP.tsv',
    chroms = 'ancestry/LC/all/{prefix}_SNP.chr'
  output:
    touch('ancestry/LC/all/data/get_ancestry_GLs_{prefix}.done')
  log: 'log/ancestry/LC/all/get_ancestry_GLs_{prefix}log'
  threads: 1
  resources: mem_mb=240000, walltime="24:00:00"
  message:
    """ Call genotype likelihoods from realigned bam files using a list of ancestry informative sites """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -b {input.bamlist} -ref {input.ref} -out ancestry/LC/all/data/{wildcards.prefix} -GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 -sites {input.sites} -rf {input.chroms} 2> {log}
    """



rule PCAngsd_ancestry:
  input: 
    touched = 'ancestry/LC/all/data/get_ancestry_GLs_{prefix}.done' 
  output:
    touch('pcangsd/ancestry/PCAngsd_ancestry_GLs_{prefix}_cov_admix.done')
  log: 'log/pcangsd/ancestry/PCAngsd_ancestry_GLs_{prefix}_cov_admix.log'
  threads: 1
  resources: mem_mb=500, walltime="00:15:00"
  message:
    """ Estimate covariance matrix and admixture proportions from ancestry informative GLs using PCAngsd """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    for i in {{2..10}}; do
        singularity exec --home $PWD:$HOME /scratch/c7701178/bio/ngs+tools.sif pcangsd --beagle ancestry/LC/all/data/{wildcards.prefix}.beagle.gz --out pcangsd/ancestry/PCAngsd_ancestry_GLs_{wildcards.prefix} --maf 0.05 --admix --admix_K $i --n_eig 10 --threads {threads} 2> {log}
    done
    """

