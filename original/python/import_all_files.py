from ntpath import join
import os
import shutil
from pathlib import Path

#Aqui vamos abrir e importar os arquivos.
path = "Z:\\MAR23\\Seq_02_03042023-385496337\\FASTQ_Generation_2023-04-04_19_48_23Z-663216607"
path_final = "Z:\\MAR23"
dir_list = os.listdir(path)

for x in dir_list:
        files=os.listdir(path + "\\" + "".join(x))
        print(files)
        for fname in files:
                shutil.copy2(os.path.join(path,x,fname), path_final)
for fq in os.listdir(path_final):
        print(fq.replace("_001",""))
        os.rename(os.path.join(path_final, fq), os.path.join(path_final, fq.replace("_001","")))
