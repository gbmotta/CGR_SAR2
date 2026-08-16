# import OS module
from ntpath import join
import os
import zipfile
import gzip

# Get the list of all files and directories
path = 'C:\BaseSpace\Seq_01_110123-378631891\FASTQ_Generation_2023-01-13_21_28_28Z-646697053'
path_final = "D:\OPAS\JAN_23"
dir_list = os.listdir(path)

for x in dir_list:
    for item in os.listdir(path + "\\" + "".join(x)): 
         if item.endswith('.gz'): 
            file_name = path + "\\" + "".join(x) + "\\" + item 
            out_name = open((path_final +'\\' + item.replace('.gz','')), 'w')
            with gzip.open(file_name,"rb") as ip_byte:
                out_name.write(ip_byte.read().decode("utf-8"))
                out_name.close()

                
