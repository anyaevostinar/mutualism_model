import os

script_name = "/mnt/home/anyaejo/Mutualism/mutualism_model/Symbulation2/tasks_generalist"

treatments = ["${PBS_ARRAYID} 0.001 0.0","${PBS_ARRAYID} 0.001 0.1","${PBS_ARRAYID} 0.001 0.2","${PBS_ARRAYID} 0.001 0.3","${PBS_ARRAYID} 0.001 0.4","${PBS_ARRAYID} 0.001 0.5","${PBS_ARRAYID} 0.001 0.6","${PBS_ARRAYID} 0.001 0.7","${PBS_ARRAYID} 0.001 0.8","${PBS_ARRAYID} 0.001 0.9","${PBS_ARRAYID} 0.001 1.0"]

#treatments = ["${PBS_ARRAYID} 0.001 0.6"]




header = "#!/bin/bash -login\n\n#PBS -l walltime=4:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=2gb\n#PBS -N mut_tests\n#PBS -t 1001-1020\n\ncd /mnt/home/anyaejo/Mutualism/host_general\nmodule load GNU/5.2\n"

#header = "#!/bin/bash -login\n\n#PBS -l walltime=4:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=2gb\n#PBS -N mut_tests\n#PBS -t 1005\n\ncd /mnt/home/anyaejo/Mutualism/host_general\nmodule load GNU/5.2\n"


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
