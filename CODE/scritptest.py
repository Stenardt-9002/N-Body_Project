

import os 
import re

# try:
#     os.system("rm answer.csv")


# except :
#     pass

# # therad_list = [6,8,10,12,14,16,20,24]
# # therad_list = [1,2,4,6,8,10,12,14,16,20,24]
# therad_list = [1,2,4,8,16,32,64,128,256,512,1024]
# numer_of_bodies = [2000,3000,4000,5000]

# therad_list = [1,2,4,8,10,12,16,32]


# for eachbody in numer_of_bodies:
#     for i in therad_list:


        # for j in range(4):
            # MAIN_CMMAND =  "OMP_NUM_THREADS="+str(i)+" ./run1 > sample2.txt"
            
            # MAIN_CMMAND =  "OMP_NUM_THREADS="+str(i)+" ./ran > sample2.txt"



# MAIN_CMMAND = "./test1 "+str(1024) + " " +str(1024)+" > captureout.txt"


# os.system(MAIN_CMMAND)

f = open("captureout.txt","r")
daata = f.read()
pos = daata.find("kernel time :")
pos2 = daata.find("\nTotal")

# strin_g = daata[pos+4:-1]
strin_g = daata[pos+13:pos2+1]

strin_g.strip()
num = float(strin_g)
print(num)


pos = daata.find("Total time:")
pos2 = daata.find("\nTotal")
strin_g = daata[pos+13:-1]
strin_g.strip()
num = float(strin_g)
print(num)




f.close()

# os.system("rm captureout.txt")

# f = open("answer.txt","a+")

    # f = open("answer.csv","a+")

    # f.write(str(num))        
    # f.write(",")
    # # f.write("\n")

    # f.close()

    # f = open("answer.txt","a+")
    # # f.write("\n\nThread "+str(i)+"\n")
    # f.write("\n")

    # f.close()

# f = open("answer.txt","a+")
# f.write("Hi there")
# f.close()
