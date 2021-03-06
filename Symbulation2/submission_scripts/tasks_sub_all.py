import os

script_name = "/mnt/home/anyaejo/Mutualism/mutualism_model/Symbulation2/tasks"

treatments = ["${PBS_ARRAYID} 0.001 0.0 0","${PBS_ARRAYID} 0.001 0.1 0","${PBS_ARRAYID} 0.001 0.2 0","${PBS_ARRAYID} 0.001 0.3 0","${PBS_ARRAYID} 0.001 0.4 0","${PBS_ARRAYID} 0.001 0.5 0","${PBS_ARRAYID} 0.001 0.6 0","${PBS_ARRAYID} 0.001 0.7 0","${PBS_ARRAYID} 0.001 0.8 0","${PBS_ARRAYID} 0.001 0.9 0","${PBS_ARRAYID} 0.001 1.0 0"]

#treatments = ["${PBS_ARRAYID} 0.001 0.4 1"]



header = "#!/bin/bash -login\n\n#PBS -l walltime=4:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=2gb\n#PBS -N mixed_tasks\n#PBS -t 1001-1020\n\ncd /mnt/home/anyaejo/Mutualism/tasks\nmodule load GNU/5.2\n"


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
