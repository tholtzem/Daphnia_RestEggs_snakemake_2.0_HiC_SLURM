rule taft:
	input:
		'TAFT/taft/info2.dat'
	output:
		touch('TAFT/taft/test.done')
	log:
		'log/taft/test.log'
	threads: 4
	resources: mem_mb=500, walltime="00:15:00"
	message:
		""" TAFT """
	shell:
		"""
		cd TAFT/taft/
		./taftx.exe <<<info2
		"""


