import json
import csv
from pathlib import Path
import os
import itertools

lista_codigos = []
lista_reads = []
lista_resultante_reads = []
results = 'D:\OPAS\JAN_23'

results_pastas = os.listdir(results)
for x in results_pastas:
     #print(x)   
     if "results" in x:
        pastas_sequencias = os.listdir(results + "\\" + "".join(x))
        for file in (pastas_sequencias):
            if "fastp.json" in file:
                #print(type(file))
                arquivos_fastp = (results + "\\" + x + "\\" + file)
                #print(arquivos_fastp)
                with open(arquivos_fastp) as jsonFile:
                    jsonObject = json.load(jsonFile)
                    jsonFile.close()

                    xx_values = arquivos_fastp.split("D:\OPAS\JAN_23\\")[1].split("_")[0]
                    reads = jsonObject['filtering_result']
                    header =['Código', "Reads"]
                    lista_codigos.append(xx_values)
                    lista_reads.append(reads['passed_filter_reads'])
                    #print(reads['passed_filter_reads'])
                    #print(xx_values)
#print(lista_reads)
#print(lista_codigos)

lista_resultante_reads = [list(sublist) for sublist in zip(lista_codigos, lista_reads)]
print(lista_resultante_reads)

with open('D:\OPAS\JAN_23\RESULTS\csv_reads.csv', 'w', encoding='UTF8', newline='') as reads_csv:
    writer = csv.writer(reads_csv)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(lista_resultante_reads) 