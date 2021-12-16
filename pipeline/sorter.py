import os
import sys
import gzip

sample_path = str(input('Input the path to files: '))
n = 0

for file in os.listdir(sample_path):
    with open(f'{sample_path}/{file}', 'r') as my_file:
        i = 0
        for line in my_file:
            if i <= 400000:
                if line != '\n':
                    i += 1
            else:
                next
                
    if i > 400000:
        with open(f'{sample_path}/tmp.fastq', 'w') as new_file:
            with open(f'{sample_path}/{file}', 'r') as my_file:
                j = 0
                for line in my_file:
                    if j < 400000:
                        new_file.write(f'{line}')
                        j += 1
                    else:
                        next
        
        with open(f'{sample_path}/tmp.fastq', 'r') as my_file:
            with open(f'{sample_path}/{file}', 'w') as new_file:
                for line in my_file:
                    new_file.write(f'{line}')
    
    n += 1
    print(n)
os.remove(f'{sample_path}/tmp.fastq')
print('The job is completed!!!')
