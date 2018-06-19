import os

script_name = "/mnt/home/anyaejo/Mutualism/mutualism_model/Symbulation2/mixed"

treatments = ["${PBS_ARRAYID} 0.001 5 0.01","${PBS_ARRAYID} 0.001 5 0.02","${PBS_ARRAYID} 0.001 5 0.03","${PBS_ARRAYID} 0.001 5 0.04","${PBS_ARRAYID} 0.001 5 0.05","${PBS_ARRAYID} 0.001 5 0.06","${PBS_ARRAYID} 0.001 5 0.07","${PBS_ARRAYID} 0.001 5 0.08","${PBS_ARRAYID} 0.001 5 0.09","${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.12", "${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.13", "${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.14", "${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.15","${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.16","${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.17","${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.18","${PBS_ARRAYID} 0.001 5 0.11","${PBS_ARRAYID} 0.001 5 0.19"]





header = "#!/bin/bash -login\n\n#PBS -l walltime=4:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=2gb\n#PBS -N mut_tests\n#PBS -t 1001-1020\n\ncd /mnt/home/anyaejo/Mutualism/vert_trans_paper/onetwenty_verttrans\nmodule load GNU/5.2\n"


for t in range(len(treatments)):
    tempfilename = 'temp_'+str(t)+'.qsub'
    tempfile = open(tempfilename, 'w')
    tempfile.write(header)
    tempfile.write(script_name)
    tempfile.write(" ")
    tempfile.write(treatments[t])
    tempfile.close()
    os.system("qsub {0}".format(tempfilename))
    print "submitted", t
    os.system("rm {0}".format(tempfilename))
