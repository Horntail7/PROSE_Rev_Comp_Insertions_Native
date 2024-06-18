import sys
import gzip
from Bio import SeqIO
import re

# Command-line arguments: input file, output directory, barcode file, native barcode file
input_file = sys.argv[1]
output_dir = sys.argv[2]
barcode_file = sys.argv[3]
chip_number = sys.argv[4]
native_bar = sys.argv[5]

# Read barcode sequence maps
def read_barcode_map(file_path):
    barcode_map = {}
    with open(file_path, 'r') as f:
        for line in f:
            name = line.strip().split(',')[0]
            barcode_map[name] = line.strip().split(',')[1:]
            
    return barcode_map

barcode_map = read_barcode_map(barcode_file)

# Create combined barcode map
#combined_barcode_map = {}
#for barcode_name, barcode_seq in barcode_map.items():
#    for native_name, native_seq in native_barcode_map.items():
#        combined_name = f"{chip_number}_{barcode_name}_{native_name}"
#        combined_seq = (barcode_seq, native_seq)
#        combined_barcode_map[combined_name] = combined_seq

# Modified demultiplex function
def demultiplex(input_file, output_dir, barcode_map):
    barcode_reads = {}
    for name, seqs in barcode_map.items():
        barcode_reads[name] = 0
        
    with gzip.open(input_file, "rt") as handle:
        output_files = {name: gzip.open(f"{output_dir}/{chip_number}_{native_bar}_{name}.fastq.gz", "wt") for name in barcode_map}

        
        for record in SeqIO.parse(handle, "fastq"):
            for name, seqs in barcode_map.items():
                
                flag=1
                for seq in seqs:
                    if seq not in record.seq:
                        flag=0
                        break
                if flag == 1:
                    SeqIO.write(record, output_files[name], "fastq")
                    barcode_reads[name]+=1
                    break

        for file in output_files.values():
            file.close()
            
        outfile=open(f"{output_dir}/{chip_number}_{native_bar}_summary.csv","wt")
        outfile.write("Barcode,Reads\n")
        for name,reads in barcode_reads.items():
            outfile.write("%s,%d\n" %(name,reads))
        outfile.close()
        
# Call the function
demultiplex(input_file, output_dir, barcode_map)


