rule get_chrom2:
	input:
		sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.list',
	output:
		chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.chr'
	message:
		""" Get a list of chromosomes/scaffolds """
	shell:
		"""
		cut -f1 {input} | uniq > {output}
		"""


rule index_sites_file:
	input:
		sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.list',
	output:
		touch('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP_index.done')
	log: 'log/ancestry/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP_index.log'
	message:
		""" Index sites file """
	shell:
		"""
		module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
		singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd sites index {input} 2> {log}
		"""


rule get_mafs_bams_sites:
  input:
    ref = config["ref_HiC"],
    #bamlist = 'synthesis_daphnia/list/df{df}_realignedbam.list',
    bams = '/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{bams}.realigned.bam',
    touched = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP_index.done',
    sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.list',
    chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.chr'
  output:
    touch('mafs/{bams}.mafs.done')
  log:
   'log/ancestry/{bams}.mafs.log'
  threads: 2
  resources: mem_mb=5000, walltime="12:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -i {input.bams} -ref {input.ref} -anc {input.ref} -out ancestry/LC/hybrids_LC/mafs/{wildcards.bams}.mafs -GL 2 -doMaf 1 -doMajorMinor 3 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd 1 -sites {input.sites} -rf {input.chrom} 2> {log}
    """
