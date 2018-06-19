import os

script_name = "/mnt/home/anyaejo/Mutualism/mutualism_model/Symbulation2/denovo"

treatments = ["${PBS_ARRAYID} 0.001 5 0.0","${PBS_ARRAYID} 0.001 5 0.1","${PBS_ARRAYID} 0.001 5 0.2","${PBS_ARRAYID} 0.001 5 0.3","${PBS_ARRAYID} 0.001 5 0.4","${PBS_ARRAYID} 0.001 5 0.5","${PBS_ARRAYID} 0.001 5 0.6","${PBS_ARRAYID} 0.001 5 0.7","${PBS_ARRAYID} 0.001 5 0.8","${PBS_ARRAYID} 0.001 5 0.9","${PBS_ARRAYID} 0.001 5 1.0", "${PBS_ARRAYID} 0.01 5 0.0","${PBS_ARRAYID} 0.01 5 0.1","${PBS_ARRAYID} 0.01 5 0.2","${PBS_ARRAYID} 0.01 5 0.3","${PBS_ARRAYID} 0.01 5 0.4","${PBS_ARRAYID} 0.01 5 0.5","${PBS_ARRAYID} 0.01 5 0.6","${PBS_ARRAYID} 0.01 5 0.7","${PBS_ARRAYID} 0.01 5 0.8","${PBS_ARRAYID} 0.01 5 0.9","${PBS_ARRAYID} 0.01 5 1.0"]





header = "#!/bin/bash -login\n\n#PBS -l walltime=4:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=2gb\n#PBS -N mut_tests\n#PBS -t 1001-1020\n\ncd /mnt/home/anyaejo/Mutualism/vert_trans_paper/denovo\nmodule load GNU/5.2\n"


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
