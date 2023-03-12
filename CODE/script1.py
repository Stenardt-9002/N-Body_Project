

import os 
import re

# try:
#     os.system("rm answer.csv")
# except :
#     pass
f1 = open("finalstuff.txt","w")

therad_list = [1,2,4,8,16,32]
# therad_list = [ 64,128,256,512,1024]

# therad_list = [512,1024]

# numer_of_bodies = [1024,2000,3000,4000,5000]
numer_of_bodies = [1024,1536]


# therad_list = [1,2,4,8,10,12,16,32]

for eachbody in numer_of_bodies:
    for i in therad_list:
        for j in range(4):
        # for j in range(2):

            # MAIN_CMMAND =  "OMP_NUM_THREADS="+str(i)+" ./run1 > sample2.txt"
            # MAIN_CMMAND =  "OMP_NUM_THREADS="+str(i)+" ./ran > sample2.txt"
            MAIN_CMMAND = "./test1 "+str(i) + " " +str(eachbody)+" > captureout.txt"
            os.system(MAIN_CMMAND)

            f = open("captureout.txt","r")

            daata = f.read()
            pos = daata.find("kernel time :")
            pos2 = daata.find("\nTotal")
            strin_g = daata[pos+13:pos2+1]

            strin_g.strip()
            num = float(strin_g)
            
            f1.write("thread list = ")
            f1.write(str(i))
            f1.write(" bodies list = ")
            f1.write(str(eachbody))
            f1.write(" kernel ")
            f1.write(str(num))
            f1.write(" ")
            pos = daata.find("Total time:")
            strin_g = daata[pos+13:-1]
            strin_g.strip()
            num = float(strin_g)
            f1.write(" total time ")
            f1.write(str(num))
            f1.write(" ")
            f1.write("\n")

            f.close()


          
f1.close()
       

    # f = open("answer.txt","a+")
    # f.write("Hi there")
    # f.close()
