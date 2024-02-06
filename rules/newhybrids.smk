rule simulate_2GENs:
  input: 
    codom = 'newhybrids/sim/codom_{sites}.in',
    GtypFreq = 'newhybrids/sim/TwoGensGtypFreq.txt'
  output:
     touch('newhybrids/sim/two_BX/{sites}sites/N_IND{N_IND}/test_{sites}sites_N_IND{N_IND}.done')
  log: 'log/test_{sites}sites_N_IND{N_IND}.log'
  threads: 4
  resources: mem_mb=500, walltime="00:15:00"
  message:
	  """ Simulate pure individuals, hybrids and backcrosses """
  shell:
    """
    mkdir newhybrids/sim/two_BX
    cd newhybrids/sim/two_BX/
    mkdir {wildcards.sites}sites
    cd {wildcards.sites}sites/
    /scratch/c7701178/bio/newhybrids/simdata_nh -g ../../../TwoGensGtypFreq.txt -c ../../../codom_{wildcards.sites}.in -i 0 z {wildcards.N_IND} -i 1 z {wildcards.N_IND} -i 2 {wildcards.N_IND} -i 3 {wildcards.N_IND} -i 4 {wildcards.N_IND} -i 5 {wildcards.N_IND} 
    """

