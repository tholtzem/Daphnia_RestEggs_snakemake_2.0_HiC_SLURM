rule get_ancestry_sites:
  input:
	  vcf = 'ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.vcf.gz'
  output:
	  sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.tsv'
  message: "Get information of chrom/pos/ref/alt from vcf"
  shell:
  	"""
	bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} | tail -n+2 > {output.sites} 2> {log}
	"""


rule get_chrom2:
	input:
		sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.tsv'
	output:
		chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.chr'
	message:
		""" Get a list of chromosomes/scaffolds """
	shell:
		"""
		cut -f1 {input.sites} | uniq > {output.chrom}
		"""


rule index_sites_file:
	input:
		sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.tsv'
	output:
		touch('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.index.done'),
	log: 'log/ancestry/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.index.log'
	message:
		""" Index sites file """
	shell:
		"""
		module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
		singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd sites index {input.sites} 2> {log}
		"""


rule get_mafs_bams_sites:
  input:
    ref = config["ref_HiC"],
    #bamlist = 'synthesis_daphnia/list/df{df}_realignedbam.list',
    bams = '/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{bams}.realigned.bam',
    touched = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.index.done',
    sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.tsv',
    chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.chr'
  output:
    touch('mafsNEW/{bams}.mafs.done')
  log:
   'log/ancestry/{bams}.mafs.log'
  threads: 2
  resources: mem_mb=5000, walltime="12:00:00"
  message:
    """ Compute 'allele frequencies' for each individual in angsd using the reference allele as major """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -i {input.bams} -ref {input.ref} -out mafsNEW/{wildcards.bams} -GL 2 -doMaf 1 -doMajorMinor 4 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -sites {input.sites} -rf {input.chrom} 2> {log}
    """


rule get_mafs_ancestry_DW:
  input:
    ref = config["ref_HiC"],
    bamlist = 'mafsNEW/list/DW_{N}_bam.list',
    touched = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.index.done',
    sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.tsv',
    chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.chr'
  output:
    touch('mafsNEW/DW/DW_{N}.mafs.done')
  log:
   'log/mafsNEW/DW/DW_{N}.mafs.log'
  threads: 1
  resources: mem_mb=10000, walltime="24:00:00"
  message:
    """ Compute allele frequencies for discrete windows in angsd using the reference allele as major """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -b {input.bamlist} -ref {input.ref} -out mafsNEW/DW/DW_{wildcards.N} -GL 2 -doMaf 1 -doMajorMinor 4 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -sites {input.sites} -rf {input.chrom} 2> {log}
    """


rule get_mafs_DW_MajorMinor3:
  input:
    ref = config["ref_HiC"],
    bamlist = 'mafsNEW/list/DW_{N}_bam.list',
    touched = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_ancestrySites_chrom_pos.index.done',
    sites = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_ancestrySites_chrom_pos_ref_alt.tsv',
    chrom = 'list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_ancestrySites_chrom_pos.chr'
  output:
    touch('mafsNEW/DW_sites3/DW_{N}.mafs.done')
  log:
   'log/mafsNEW/sites3/DW_{N}.mafs.log'
  threads: 1
  resources: mem_mb=10000, walltime="24:00:00"
  message:
    """ Compute allele frequencies for discrete windows in angsd using the reference allele as major """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -b {input.bamlist} -ref {input.ref} -out mafsNEW/DW_sites2/DW_{wildcards.N} -GL 2 -doMaf 1 -doMajorMinor 3 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -sites {input.sites} -rf {input.chrom} 2> {log}
    """
