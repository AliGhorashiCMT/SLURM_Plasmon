# SLURM_Plasmon

A few notes-

This has been written for a particular DFT calculation. The volume (area) of whatever 2d material have to be altered for different calculations. Also the filenames are by JJDFTX.jl conventions. Note that this code is from JJDFTX.jl but specialized to just plasmonics calculations through SLURM. The basic idea is to create a job array with different wavevectors being solved by different arrays. 
