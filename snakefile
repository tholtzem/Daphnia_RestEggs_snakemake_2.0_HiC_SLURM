include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		#expand('raw/qc/fastqc/{sample}_{pair}_001_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
		#expand('raw/qc/multiqc/'),
		## !!! works, but some intermediate files from trimming were deleted due to space limitation on mach2 !!! ##
		#expand('trm/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmd.fq.gz']),
		#expand('trm/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmdfilt.fq.gz']),
		#expand('trm/qc/fastqc/{sample}_{pair}.trmdfilt_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
		#expand('trm/qc/multiqc/')#,
		## !!! works, but some intermediate files from Kraken were deleted due to space limitation on mach2 !!! ##
		#expand('KRAKEN2_RESULTS/{sample}_kraken2_classified.done', sample=sample_names),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['trmdfilt.classified_retain.fq']),
		#expand('KRAKEN2_RESULTS/{sample}_{pair}.{ext}', sample=sample_names, pair=['R1','R2'], ext=['trmdfilt.keep.fq.gz']),
		#expand('KRAKEN2_RESULTS/qc/fastqc/{sample}_{pair}.trmdfilt.keep_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['zip', 'html']),
		#expand('KRAKEN2_RESULTS/qc/multiqc/'),
		#expand('ref/bb_indexRef.done'),
		#expand('bbmap/{sample}.{ext}', sample=sample_names, ext=['bam']),
		#expand('bbmap/{sample}.{ext}', sample=sample_names, ext=['bam.bai']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['dedup.bam', 'dedup.metrics.txt']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['overlapclipped.bam']),
		#expand('deDup/{sample}.{ext}', sample=sample_names, ext=['overlapclipped.bam.bai']),
		#exp#and('ref/{ref_name}.{ext}', ref_name=['Dgaleata_M5_PBasm.FINAL'] ,ext=['fasta.fai', 'dict']),
		#expand('deDup/{sample}.symlink.done', sample=outgroup_names),
		#expand('deDup/{sample}.added2ClippedList.done', sample=outgroup_names),
		#expand('deDup/{sample}.REFsadded2ClippedList.done', sample=refclone_names),
		##Run until here
		#expand('list/indels_{ref_name}.list', ref_name=['Dgaleata_M5_PBasm.FINAL']),
		#expand('list/indels_{ref_name}_curvi.list', ref_name=['Dgaleata_M5_PBasm.FINAL']),
		#expand('realigned/{sample}.{ext}', sample=outgroup_names, ext=['realigned.bam']),
                #expand('depth/{sample}.{ext}', sample=outgroup_names, ext=['realigned.bam.coverage.hist']),
		#expand('depth/{sample}.{ext}', sample=outgroup_names, ext=['realigned.bam.depth.gz']),
		#expand('/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{sample}.{ext}', sample=refclone_names, ext=['realigned.bam.MAPQ.gz']),
		#expand('depth/{sample}.added2DepthList.list', sample=outgroup_names),
		#expand('depth/stats/depth_statistics_curvi.txt'),
		#expand('depth/plots/plot_summary_{sets}.done', sets=['LC', 'LZ', 'LCwithoutREF', 'LZwithoutREF', 'longREF7', 'longREF9', 'cucREF7', 'galREF9', 'ALL_gal_LC', 'hatched_gal_LC', 'EGGS_gal_LC', 'PEL_cuc_LC', 'PEL_long_LC', 'PRE_long_LC', 'curvirostris']),
                ## run until here
		#expand('depth/stats/rbind_depthFilter.done'),
		#expand('bedtools/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.bed']),
		#expand('bedtools/plots/{sample}.realigned.{ext}', sample=ALL_names, ext=['genomecov.pdf']),
		#expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.list', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('angsd/{sets}/index_globalSNP_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.chr', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_cov_admix.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_geno.beagle.gz', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub_pos.gz', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('ngsLD/{sets}/run_ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure),
		#expand('ngsLD/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.1', '0.1', '0.1', '0.1']),
		#expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.list', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.1', '0.1', '0.1', '0.1']),
		#expand('ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.chr', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.1', '0.1', '0.1', '0.1']),
		#expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_geno.beagle.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
		#expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_pos.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
		#expand('LD_decay/{sets}/ngsLD_windows_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.ld.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
                #expand('saf/saf_samples/{sample}.GL2.saf.idx.done', sample=ALL_names),
                #expand('saf/saf_samples/{sample}.GL2.est.ml', sample=ALL_names),
		#expand('list/saf/{POPS}.list', POPS=sets_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window{bp}bp.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf, bp=['10000', '10000', '10000', '10000', '10000', '10000', '10000', '10000', '10000', '10000'])
		# , '20000', '20000', '20000', '20000', '20000', '20000', '20000', '20000', '30000', '30000', '30000', '30000', '30000', '30000', '30000', '30000', '50000', '50000', '50000', '50000', '50000', '50000', '50000', '50000'
		#expand('saf/POPS/{POP1}_vs_{POP2}_2D.sfs', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fstIndex.done', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done', zip, POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fst.global.tsv', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fst.Window10kb.type{types}.tsv', zip, POP1=POP1, POP2=POP2, types=['1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.global.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.Window10kb.type{types}.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC', types=['1', '2', '3']),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.A_B.tsv.gz', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.A_B.tsv.gz', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('ngsRelate/{POPS}.glf3.freq', POPS=['PEL_cuc_LC', 'PEL_long_LC', 'PRE_long_LC', 'PEL_F1_LC']),
		#expand('ngsRelate/{POPS}.glf3.freq', POPS=sets_saf),
		#expand('ngsRelate/{POPS}.glf3.stats.txt', POPS=sets_saf),
		#expand('mafs/{POPS}_1000random_SNPs.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_1000random_sites.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_1000random_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_top10_NOTfixed_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_10random_fixed_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('list/abbababa/4pop_combinations.bam.list.done'),
		#expand('list/abbababa/4pop_combinations.sizefile.done'),
		#expand('abbababa/{prefix}/{combi}.abbababa.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi[4:12], blocksize=['2500000']),
		#expand('list/abbababa/errorFile.txt'),
		#expand('list/abbababa/4pop_combinations.Popnames.list.done'),
		#expand('abbababa/{prefix}/{combi}.Dstat_results.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi[4:12], blocksize=['2500000']),
		##expand('realigned/{sample}.{ext}', sample=samples_ALL_names, ext=['realigned.bai']),
		#expand('saf/POPS/sites/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sites.saf.idx.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/sites/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sites.1D.sfs', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/sites/{POP1}_vs_{POP2}.sites.2D.sfs', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/2pops/{POP1}_vs_{POP2}.sites.fstIndex.done', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/2pops/{POP1}_vs_{POP2}.sites.fst.global.tsv', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fstIndex.done', zip, POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fst.global.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fst.Window{bp}bp.type{types}.tsv', POP1=['PRE_long_LC'], POP2=['ALL_gal_LC'], POP3=['PEL_long_LC'], types=['2'], bp=['10000', '20000', '30000', '50000', '100000']),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.A_B.tsv.gz', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('synthesis_daphnia/list/PRE_long_LC_ALL_gal_LC_PEL_long_LC.SNPs.chr'),
		#expand('synthesis_daphnia/list/PRE_long_LC_ALL_gal_LC_PEL_long_LC.SNPs_index.done'),
		#expand('synthesis_daphnia/mafs/{bams}.mafs.done', bams=bams),
		#expand('synthesis_daphnia/mafs/df{df}.mafs.done', df=range(1,19)),
		#expand('ancestry/LC/data/get_bcf_{prefix}_MajorMinor4.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['vcf.gz']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['INFODP.txt', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['INFODP.pdf', 'INFODP.stats.txt']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_hf.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['vcf.gz', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_hf.DP{genoDP}.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf.gz', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_hf.DP{genoDP}.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['imiss']),
		#expand('ancestry/LC/data/{prefix}_MajorMinor4_hf.DP{genoDP}.imiss.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['pdf', 'tsv']),
		#expand('ancestry/LC/data/imissRM/{prefix}_MajorMinor4_hf.DP{genoDP}.imissRM.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf.gz']),
		#expand('ancestry/LC/data/imissRM/{prefix}_MajorMinor4_hf.DP{genoDP}.imissRM.vmiss20.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf', 'nbrSites.txt']),
		#expand('ancestry/LC/data/imissRM/{prefix}_MajorMinor4_hf.DP{genoDP}.imissRM.vmiss20.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['bed'])
		#expand('ancestry/LC/new/{pop1}_{pop2}.{ext}', pop1=['LONGoverall_0.99'], pop2=['GALoverall'], ext=['fixed_sites.txt', 'fixed_sites.report.tsv']),
		#expand('ancestry/LC/new/{pop1}_{pop2}.{ext}', pop1=['LONGoverall_0.99'], pop2=['GALoverall'], ext=['fixed_sites_100percent_thinned1000bp.svg', 'fixed_sites_100percent_thinned1000bp.report.tsv']),
		#expand('ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent.{ext}', pop1=['LONGoverall0.99'], pop2=['GALoverall0.99'], pct=['100'], ext=['txt', 'report.tsv']),
		#expand('ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_thinned1000bp.{ext}', pop1=['LONGoverall0.99'], pop2=['GALoverall0.99'], pct=['100'], ext=['svg', 'report.tsv']),
		expand('ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}.{ext}', pop1=['LONGoverall0.99'], pop2=['GALoverall0.99'], pct=['99'], ext=['txt', 'report.tsv']),
		expand('ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_thinned1000bp.{ext}', pop1=['LONGoverall0.99'], pop2=['GALoverall0.99'], pct=['99'], ext=['svg', 'report.tsv']),
		expand('ancestry/LC/overall/{pop1}_{pop2}.fixed_sites_{pct}percent_sufficient{pct}_{ext}', pop1=['LONGoverall0.99'], pop2=['GALoverall0.99'], pct=['99'], ext=['SNP.tsv']),
		expand('ancestry/LC/data/imissRM/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites{pct}.{ext}', pct=['99'], ext=['vcf.gz']),
		expand('ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4.{ext}', DW=range(1,19), pct=['99'], ext=['vcf.gz']),
		expand('ancestry/LC/data/imissRM/DW_{DW}_ancestrySites{pct}_MajorMinor4_new.{ext}', DW=range(1,19), pct=['99'], ext=['vcf.gz']),
		expand('ancestry/LC/data/imissRM/DW_{DW}_AF_ancestrySites{pct}_MajorMinor4_new.{ext}', DW=range(1,19), pct=['99'], ext=['tsv']),
		#expand('ref/ref_HiC_chromsize.tsv')
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.{ext}', ext=['tsv']),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.{ext}', ext=['chr']),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_MajorMinor4_hf.DP10.imissRM.vmiss20_ancestrySites.index.done'),
		#expand('mafsNEW/DW/DW_{N}.mafs.done', N=range(1,19)),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window10kb.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('list/abbababa/4pop_combinations.bam.list.done'),
		#expand('list/abbababa/4pop_combinations.sizefile.done'),
		#expand('abbababa/{prefix}/{combi}.abbababa.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi2, blocksize=['2500000']),
		#expand('list/abbababa/errorFile.txt'),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.chr'),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP_index.done'),
		#expand('ancestry/LC/unfiltered/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites.txt', 'fixed_sites.report.tsv']),
		#expand('ancestry/LC/notfiltered/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites_100percent_thinned1000bp.svg', 'fixed_sites_100percent_thinned1000bp.report.tsv']),
		#expand('list/abbababa/4pop_combinations.Popnames.list.done'),
		#expand('abbababa/{prefix}/{combi}.Dstat_results.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi2, blocksize=['2500000']),
		#expand('ancestry/LC/all/index_{prefix}_SNP.done', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('ancestry/LC/all/{prefix}_SNP.chr', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('ancestry/LC/all/data/get_ancestry_GLs_{prefix}.done', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('pcangsd/ancestry/PCAngsd_ancestry_GLs_{prefix}_cov_admix.done', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('ancestry/LC/data/imissRM/{prefix}/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_{prefix}.vcf', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('ancestry/LC/data/imissRM/{prefix}/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_hf.DP10.imissRM.vmiss20_{prefix}_selected_samples.vcf', prefix=['longispina7_galeata9.fixed_sites', 'longispina7_galeata9.fixed_sites_100completeness', 'longispina7_galeata9.fixed_sites_100completeness_1000bp']),
		#expand('TAFT/taft/test.done'),
		#expand('newhybrids/sim/two_BX/{sites}sites/N_IND{N_IND}/test_{sites}sites_N_IND{N_IND}.done', sites=['1000'], N_IND=['25']),
		#expand('mafsNEW/{bams}.mafs.done', bams=bams),
		#expand('mafsNEW/DW_sites2/DW_{N}.mafs.done', N=range(1,19)),
		#expand('mafsNEW/DW_sites3/DW_{N}.mafs.done', N=range(1,19))
		#
		


# -----------------------------------------------


#include: "rules/eggs2.0.smk"
#include: "rules/samtools_MAPQ.smk"
#include: "rules/get_beagle_allSNPs.smk"
#include: "rules/PCAngsd_allSNPs.smk"
#include: "rules/ngsLD.smk"
#include: "rules/get_beagle_LDprunedSNPs.smk"
#include: "rules/PCAngsd_LDprunedSNPs.smk"
#include: "rules/estimate_ngsLD.smk"
#include: "rules/realSFS.smk"
#include: "rules/get_allelefreq_perpop.smk"
#include: "rules/ngsRelate.smk"
#include: "rules/realSFS_Fst_sites.smk"
#include: "rules/abbababa2.smk"
#include: "rules/get_vcf.smk"
#include: "rules/ancestry.smk"
include: "rules/ancestry_sufficientfixed.smk"
#include: "rules/faidx.smk"
#include: "rules/get_allelefreq_perpop_sites.smk"
#include: "rules/ancestry_PCA.smk"
#include: "rules/taft.smk"
#include: "rules/newhybrids.smk"
