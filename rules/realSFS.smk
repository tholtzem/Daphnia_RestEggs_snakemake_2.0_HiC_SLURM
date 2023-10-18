# Note: The estimation of the SFS is a 2-step procedure: first, we generate a ".saf.idx" file (site allele frequency likelihood) using angsd (-anc -doSaf 1), followed by an optimization of the .saf.idx file which will estimate the Site frequency spectrum (SFS) ".sfs" using realSFS.  Folding should now be done in realSFS and not in the saf file generation (-fold 1).

#rule get_SAF_POP_list:
#  input:
#    'depth/stats/{POPS}_realignedBAM_df1.list'
#  output:
#    'list/saf/{POPS}.list'
#  log:
#   'log/saf/get_SAF_{POPS}_list.log'
#  threads: 12
#  message:
#    """ copy population files for saf """
#  shell:
#    """
#    cp {input} {output} 2> {log}
#    """

rule angsd_SAF_POP:
  input:
    ref = config["ref_HiC"],
    bamlist = 'list/saf/{POPS}.list'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done')
  log:
   'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.log'
  threads: 2
  resources: mem_mb=80000, walltime="24:00:00"
  message:
    """ Compute site allele frequency likelihood (.saf.idx) for each population using angsd """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/angsd -b {input.bamlist} -ref {input.ref} -anc {input.ref} -out saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -minInd {wildcards.IND} -setMinDepth {wildcards.MinDepth} -setMaxDepth {wildcards.MaxDepth} -GL 2 -doCounts 1 -doSaf 1 2> {log}
    """


rule onedsfs_folded:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done'
  output:
    sfs = 'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs'
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs.log'
  threads: 2
  resources: mem_mb=10000, walltime="01:00:00"
  message:
    """ Optimize .saf.idx and calculate 1dsfs's using realSFS and fold for theta estimates """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/misc/realSFS saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.saf.idx -fold 1 > {output.sfs} 2> {log}
    """


rule sfs2theta:
  input:
    touched = 'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done',
    sfs = 'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done')
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.log'
  threads: 2
  resources: mem_mb=800, walltime="1:00:00"
  message:
    """ calculate theta estimates for each site """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/misc/realSFS saf2theta saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.saf.idx -sfs {input.sfs} -outname saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2 2> {log} 2> {log}
    """


rule theta_stat_chrom:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done')
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.log'
  threads: 2
  resources: mem_mb=100, walltime="00:10:00"
  message:
    """ Estimate Tajimas D and other statistics """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/misc/thetaStat do_stat saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.thetas.idx 2> {log}
    """


rule theta_stat_SW:
  input:
    'saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done'
  output:
    touch('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window{bp}bp.done')
  log:
    'log/saf_POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window{bp}bp.log'
  threads: 2
  resources: mem_mb=600, walltime="06:00:00"
  message:
    """ Estimate Tajimas D and other statistics using a sliding window analysis """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/angsd+ngsrelate.sif /opt/angsd-0.939/misc/thetaStat do_stat saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.thetas.idx -win {wildcards.bp} -step {wildcards.bp} -outnames saf/POPS/{wildcards.POPS}_{wildcards.IND}_{wildcards.MinDepth}_{wildcards.MaxDepth}.GL2.theta_stats_Window{wildcards.bp}bp 2> {log}
    """


rule twodsfs_folded:
  input:
    'saf/POPS/PEL_cuc_LC_8_8_215.GL2.1D.sfs'
    #pop2 = 'saf/POPS/df_{POP2}.GL2.saf.idx.done'
    #pop1 = 'saf/POPS/{POP1}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done',
    #pop2 = 'saf/POPS/{POP2}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done'
  output:
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}_2D.sfs'
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}_2D.sfs.log'
  threads: 2
  resources: mem_mb=10000, walltime="48:00:00"
  message:
    """ Optimize .saf.idx and calculate all pairwise 2dsfs's (joint site frequency spectra) using realSFS and fold for FST """
  shell:
    """ 
    realSFS saf/POPS/{wildcards.POP1}*.GL2.saf.idx saf/POPS/{wildcards.POP2}*.GL2.saf.idx -fold 1 > {output.pop_pair} 2> {log.pop_pair}
    """


rule Fst_index_2pops:
  input:
    #pop1 = 'saf/POPS/{POP1}.GL2.saf.idx.done',
    #pop2 = 'saf/POPS/{POP2}.GL2.saf.idx.done',
    pop_pair = 'saf/POPS/{POP1}_vs_{POP2}_2D.sfs'
  output:
    pop_pair = touch('saf/POPS/2pops/{POP1}_vs_{POP2}.fstIndex.done')
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}.fstIndex.log'
  threads: 4
  resources: mem_mb=1000, walltime="6:00:00"
  message:
    """ Index fst for easy window analysis (2 pops) """
  shell:
    """
    realSFS fst index saf/POPS/{wildcards.POP1}*.GL2.saf.idx saf/POPS/{wildcards.POP2}*.GL2.saf.idx -sfs {input.pop_pair} -fstout saf/POPS/2pops/{wildcards.POP1}_vs_{wildcards.POP2} -whichFst 1 2> {log.pop_pair}
    """


rule Fst_index_3pops:
  input:
    pop_pair1 = 'saf/POPS/{POP1}_vs_{POP2}_2D.sfs',
    pop_pair2 = 'saf/POPS/{POP1}_vs_{POP3}_2D.sfs',
    pop_pair3 = 'saf/POPS/{POP2}_vs_{POP3}_2D.sfs'
  output:
    touch('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done')
  log:
    pop_pair = 'log/saf_POPS/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.log'
  threads: 4
  resources: mem_mb=1000, walltime="6:00:00"
  message:
    """ Index fst for easy window analysis """
  shell:
    """
    realSFS fst index saf/POPS/{wildcards.POP1}*.GL2.saf.idx saf/POPS/{wildcards.POP2}*.GL2.saf.idx saf/POPS/{wildcards.POP3}*.GL2.saf.idx -sfs {input.pop_pair1} -sfs {input.pop_pair2} -sfs {input.pop_pair3} -fstout saf/POPS/3pops/{wildcards.POP1}_vs_{wildcards.POP2}_vs_{wildcards.POP3} -whichFst 1 2> {log.pop_pair}
    """


rule Fst_global:
  input:
    'saf/POPS/2pops/{POP1}_vs_{POP2}.fstIndex.done'
  output:
    'saf/POPS/2pops/{POP1}_vs_{POP2}.fst.global.tsv'
  log:
    'log/saf_POPS/{POP1}_vs_{POP2}.fst_Global.log'
  threads: 2
  #resources: mem_mb=, walltime="6:00:00"
  message:
    """ Get the global Fst estimate """
  shell:
    """
    realSFS fst stats saf/POPS/2pops/{wildcards.POP1}_vs_{wildcards.POP2}.fst.idx > {output} 2> {log}
    """


rule Fst_global_3pops:
  input:
    'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done'
  output:
    'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.global.tsv'
  log:
    'log/saf_POPS/{POP1}_vs_{POP2}_vs_{POP3}.fst_Global.log'
  threads: 2
  #resources: mem_mb=, walltime="6:00:00"
  message:
    """ Get the global Fst estimate """
  shell:
    """
    realSFS fst stats saf/POPS/3pops/{wildcards.POP1}_vs_{wildcards.POP2}_vs_{wildcards.POP3}.fst.idx > {output} 2> {log}
    """




rule Fst_windows_3pops:
  input:
    'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done'
  output:
    'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.Window10kb.type{types}.tsv'
  log:
    'log/saf_POPS/{POP1}_vs_{POP2}_vs_{POP3}.fst.Window10kb.type{types}.log'
  threads: 2
  resources: mem_mb=100, walltime="00:30:00"
  message:
    """ Get the Fst estimates per window """
  shell:
    """
      realSFS fst stats2 saf/POPS/3pops/{wildcards.POP1}_vs_{wildcards.POP2}_vs_{wildcards.POP3}.fst.idx -win 10000 -step 1000 -type {wildcards.types} > {output} 2> {log}
    """



rule Fst_print:
	input:
		#'saf/POPS/2pops/{POP1}_vs_{POP2}.fstIndex.done'
		'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done'
	output:
		#'saf/POPS/2pops/{POP1}_vs_{POP2}.A_B.tsv.gz'
		'saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.A_B.tsv.gz'
	log:
		#'log/saf_POPS/{POP1}_vs_{POP2}.A_B.log'
		'log/saf_POPS/{POP1}_vs_{POP2}_vs_{POP3}.A_B.log'
	threads: 2
	resources: mem_mb=1000, walltime="00:10:00"
	message: """ Get the numerator and denominator of Fst per site """
	shell:
		"""
		realSFS fst print saf/POPS/3pops/{wildcards.POP1}_vs_{wildcards.POP2}_vs_{wildcards.POP3}.fst.idx | gzip > {output} 2> {log}
		"""

#realSFS fst print saf/POPS/2pops/{wildcards.POP1}_vs_{wildcards.POP2}.fst.idx | gzip > {output} 2> {log}











#rule plotHeterozygosity_samples:
#  input:
#    'saf/{sample}.est.ml'
#  output:
#    'saf/{sample}.heterozygosity.txt'
#  log:
#    'log/{sample}.heterozygosity.log'
#  message:
#    """ Plot heterozygosity of single samples """
#  shell:
#    """
#   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
#    """



#
#
#rule plotHeterozygosity_pops:
#  input:
#    'saf/{pops}.saf.est.ml'
#  output:
#    'saf/{pops}.heterozygosity.txt'
#  log:
#    'log/{pops}.heterozygosity.log'
#  message:
#    """ Plot heterozygosity """
#  shell:
#    """
#   Rscript scripts/plotHeterozygosity.R {input} {output} 2> {log}
#    """
#
#
#rule angsd_saf_samples:
#  input:
#    ref = config["ref_rapid"],
#    bam = 'realigned/{sample}.realigned.bam'
#  output:
#    touch('saf/saf_samples/{sample}.GL2.saf.idx.done')
#  log:
#    'log/saf_samples/{sample}.GL2.saf.idx.log'
#  threads: 12
#  message:
#    """ Compute site allele frequency likelihood (.saf.idx) for single samples using angsd (required for the global heterozygosity estimate for single samples)"""
#  shell:
#    """
#    module load angsd/0.938
#    /apps/uibk/bin/sysconfcpus -n 12 angsd -i {input.bam} -out saf/saf_samples/{wildcards.sample}.GL2 -anc {input.ref} -ref {input.ref} -doSaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minQ 20 -minMapQ 30 -GL 2 2> {log}
#    """
#
#
#rule real_SFS_samples:
#  input:
#    'saf/saf_samples/{sample}.GL2.saf.idx.done'
#  output:
#    'saf/saf_samples/{sample}.GL2.est.ml'
#  log:
#    'log/saf_samples/{sample}.GL2.est.ml.log'
#  threads: 12
#  message:
#    """ Optimize .saf.idx and estimate the site frequency spectrum (SFS) using realSFS and fold (required for the global heterozygosity estimate for single samples) """
#  shell:
#    """
#    module load singularity/2.x
#    module load angsd/0.938
#    singularity exec /apps/uibk/angsd/0.938/angsd.sandbox /opt/angsd/misc/realSFS saf/saf_samples/{wildcards.sample}.GL2.saf.idx -fold 1 > {output} 2> {log}
#    """
