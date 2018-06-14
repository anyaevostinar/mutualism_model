import os.path
import gzip

folder = ''

#treatment_postfixes = ["mut0.01_mult10_vert0._start0.5","mut0.01_mult10_vert0.1_start0.5", "mut0.01_mult10_vert0.2_start0.5", "mut0.01_mult10_vert0.3_start0.5", "mut0.01_mult10_vert0.4_start0.5", "mut0.01_mult10_vert0.5_start0.5",  "mut0.01_mult10_vert0.6_start0.5", "mut0.01_mult10_vert0.7_start0.5", "mut0.01_mult10_vert0.8_start0.5", "mut0.01_mult10_vert0.9_start0.5",  "mut0.01_mult10_vert1._start0.5",  "mut0.1_mult10_vert0.5_start0.5", "mut0.01_mult10_vert0.5_start0.","mut0.01_mult10_vert0.5_start1.", "mut0.01_mult1_vert0.5_start0.5", "mut0.01_mult5_vert0.5_start0.5", "mut0.01_mult5_vert0._start0.5","mut0.01_mult5_vert0.1_start0.5", "mut0.01_mult5_vert0.2_start0.5", "mut0.01_mult5_vert0.3_start0.5", "mut0.01_mult5_vert0.4_start0.5", "mut0.01_mult5_vert0.5_start0.5",  "mut0.01_mult5_vert0.6_start0.5", "mut0.01_mult5_vert0.7_start0.5", "mut0.01_mult5_vert0.8_start0.5", "mut0.01_mult5_vert0.9_start0.5",  "mut0.01_mult5_vert1._start0.5"]

treatment_postfixes = ["mut0.001_mult5_vert0._start0.","mut0.001_mult5_vert0.1_start0.","mut0.001_mult5_vert0.2_start0.","mut0.001_mult5_vert0.3_start0.","mut0.001_mult5_vert0.4_start0.","mut0.001_mult5_vert0.5_start0.","mut0.001_mult5_vert0.6_start0.","mut0.001_mult5_vert0.7_start0.","mut0.001_mult5_vert0.8_start0.","mut0.001_mult5_vert0.9_start0.","mut0.001_mult5_vert1._start0."]

reps = range(1001, 1021)

header = "uid treatment rep update donate count Partner\n"
#header = "uid treatment rep update host_donate sym_donate host_count sym_count\n"

outputFileName = "munged_alltogether.dat"

outFile = open(outputFileName, 'w')
outFile.write(header)



for t in treatment_postfixes:
    for r in reps:
#        try:
#        fname = folder+"/"+t + "_" + str(r) + ".csv"
        fname = "avg_donation_"+str(r)+"_"+t + ".csv"
        uid = folder+"_"+t + "_" + str(r)
        curFile = open(fname, 'r')
        for line in curFile:
            if (line[0]!="U"):
                
                splitline = line.split(', ')
                if int(splitline[0]) % 1000 == 0:
                    outstring1 = "{} {} {} {} {} {} {}\n".format(uid,t,r,splitline[0], splitline[1], splitline[3], "Host" )
                    outstring2 = "{} {} {} {} {} {} {} \n".format(uid, t, r, splitline[0], splitline[2], splitline[4].strip("\n"), "Symbiont")
#                    outstring1 = "{} {} {} {} {} {} {} {}\n".format(uid,t,r,splitline[0], splitline[1], splitline[2], splitline[3], splitline[4].strip("\n" ))
                    outFile.write(outstring1)
                    outFile.write(outstring2)
                        
        curFile.close()
            
#except:
 #           print fname, ' error'
outFile.close()

