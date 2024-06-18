import argparse
import os
import csv
import gzip
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_OPEN_SCORE = -1
GAP_EXTEND_SCORE = -1
DEGENERATE_MATCH_SCORE = 0

def create_substitution_matrix():
    matrix = substitution_matrices.load("NUC.4.4")
    for base in "ATCG":
        matrix[base, base] = MATCH_SCORE
        for other_base in "ATCG":
            if base != other_base:
                matrix[base, other_base] = MISMATCH_SCORE
        matrix[base, "N"] = DEGENERATE_MATCH_SCORE
        matrix["N", base] = DEGENERATE_MATCH_SCORE
    matrix["N", "N"] = DEGENERATE_MATCH_SCORE
    return matrix

def align_sequence(seq, template, matrix):
    if not seq or not template:
        return 0
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = GAP_OPEN_SCORE
    aligner.extend_gap_score = GAP_EXTEND_SCORE
    return aligner.score(seq, template)

def find_pattern(seq, pattern, matrix, cut_off):
    if not seq or not pattern:
        return -1
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = GAP_OPEN_SCORE
    aligner.extend_gap_score = GAP_EXTEND_SCORE

    alignments = aligner.align(seq, pattern)
    best_score = 0
    best_start = -1

    for alignment in alignments:
        if alignment.score > best_score:
            best_score = alignment.score
            best_start = alignment.path[0][0]

    max_possible_score = len(pattern) * MATCH_SCORE
    threshold = cut_off * max_possible_score

#    print(seq)
#    print(pattern)
#    print(best_score,best_start)

    return best_start if best_score >= threshold else -1

def check_sequential_patterns(sequence, pattern_ids, templates, matrix):
    current_pos = 0
    for pattern_id in pattern_ids:
        if current_pos >= len(sequence):
            return False
        pattern = templates[pattern_id]
        pattern_pos = find_pattern(sequence[current_pos:], pattern, matrix)
        if pattern_pos != -1:
            current_pos += pattern_pos + len(pattern)
        else:
            return False
    return True

def process_reads(reads, templates, min_length, max_length, score_cut_off):
    substitution_matrix = create_substitution_matrix()
    non_sequential_summary = {key: 0 for key in templates}
    sequential_summary = {key: 0 for key in templates}
    non_sequential_summary["reads_not_meeting_length_criteria"] = 0
    sequential_summary["reads_not_meeting_length_criteria"] = 0
    aligned_reads = []
    bc_err_pair = {}
    
    for read in reads:
        sequence = str(read.seq)
        if not (min_length <= len(sequence) <= max_length):
            non_sequential_summary["reads_not_meeting_length_criteria"] += 1
            sequential_summary["reads_not_meeting_length_criteria"] += 1
            continue

        matched = False

        for template_name, template_seq in templates.items():
            score = align_sequence(sequence, template_seq, substitution_matrix)
            if score >= (score_cut_off * len(template_seq.replace("N", "")) * MATCH_SCORE):
                non_sequential_summary[template_name] += 1
                matched = True

        if not matched:
            non_sequential_summary["reads_not_meeting_length_criteria"] += 1
            sequential_summary["reads_not_meeting_length_criteria"] += 1

        flag = 1
        starts=[0]
        ends=[0]
        for template_name, template_seq in templates.items():
            current_pos = find_pattern(sequence, template_seq, substitution_matrix, score_cut_off)
            if current_pos == -1:
                flag = 0
                break
            else:
                sequential_summary[template_name]+=1
                sequence=sequence[current_pos+len(template_seq):]
                starts.append(ends[-1]+current_pos)
                ends.append(ends[-1]+current_pos+len(template_seq))
        if flag == 1:
            aligned_reads.append(read)
            sequence=str(read.seq)
            snips=[]
            for i in range(1,len(starts)-1):
                snips.append(sequence[ends[i]+1:starts[i+1]])

            temp_bc=''
            for snip in snips[:-1]:
                temp_bc+=snip

            if temp_bc in bc_err_pair:
                bc_err_pair[temp_bc].append(snips[-1])
            else:
                bc_err_pair[temp_bc]=[snips[-1]]

    return non_sequential_summary, sequential_summary, aligned_reads, bc_err_pair

def write_summary_to_csv(filename, summary, total_reads):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Name', 'Count', 'Percentage'])
        for key, count in summary.items():
            percentage = 0 if total_reads == 0 else (count / total_reads) * 100
            writer.writerow([key, count, "{:.2f}%".format(percentage)])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--iq', '--input_fastq', type=str, required=True, help='Full path to the fastq file')
    parser.add_argument('--ia', '--input_fasta', type=str, required=True, help='Full path to the fasta file')
    parser.add_argument('--min_length', type=int, default=25, help='Minimum length of the sequence')
    parser.add_argument('--max_length', type=int, default=1000, help='Maximum length of the sequence')
    parser.add_argument('--score_cut_off', type=float, default=0.8, help='Alignment Score')
    parser.add_argument('--od', '--output_dir', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--chip', type=str, required=False, help='Chip identifier (e.g., ChipX)', default=None)
    parser.add_argument('--native', type=str, required=False, help='Native barcode identifier (e.g., barcodeX)', default=None)
    args = parser.parse_args()

# Determine the chip identifier from the input file name if not provided
    if not args.chip:
        filename_parts = os.path.basename(args.iq).split('_')
        if len(filename_parts) >= 3:
            args.chip = filename_parts[0]  # Assuming the filename starts with the chip identifier

    if not args.chip:
        raise ValueError("Chip identifier could not be determined. Please provide it using --chip.")

    os.makedirs(args.od, exist_ok=True)

#    sequential_patterns = {
#        "ABA": ["AB", "BA"], "ABAA": ["AB", "BA", "AA"], "ABGb": ["AB", "BGb"],
#        "AAGb": ["AA", "AGb"], "ABA_not_Gb": ["AB", "BA"], "AB_not_AGb": ["AB"],
#        "ABAGb": ["AB", "BA", "AGb"], "ABCGb": ["AB", "BC", "CGb"]
#    }

    base_name = os.path.basename(args.iq).split('.')[0]
    # Check if the base name already starts with the chip identifier, to avoid redundancy
    if not base_name.startswith(args.chip):
        non_sequential_csv_filename = os.path.join(args.od, f"{args.chip}_{args.native}_{base_name}_non_sequential_summary.csv")
        sequential_csv_filename = os.path.join(args.od, f"{args.chip}_{args.native}_{base_name}_sequential_summary.csv")
    else:
        non_sequential_csv_filename = os.path.join(args.od, f"{base_name}_non_sequential_summary.csv")
        sequential_csv_filename = os.path.join(args.od, f"{base_name}_sequential_summary.csv")

    templates = {}
    with open(args.ia, "rt") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            templates[record.id] = str(record.seq)

    with gzip.open(args.iq, "rt") as fastq_file:
        reads = list(SeqIO.parse(fastq_file, "fastq"))

    total_reads = len(reads)
    non_sequential_summary, sequential_summary, aligned_reads, bc_err_pair  = process_reads(
        reads, templates, args.min_length, args.max_length, args.score_cut_off)

    maximum = -1
    for key in bc_err_pair:
        if len(bc_err_pair[key]) > maximum:
            maximum = len(bc_err_pair[key])
    print(f"No. of barcodes: {len(bc_err_pair)}")
    print(f"Maximum reads per barcode: {maximum}")
        
    write_summary_to_csv(non_sequential_csv_filename, non_sequential_summary, total_reads)
    print(f"Non-sequential results written to {non_sequential_csv_filename}")

    write_summary_to_csv(sequential_csv_filename, sequential_summary, total_reads)
    print(f"Sequential results written to {sequential_csv_filename}")

    # Always create the aligned FASTQ.GZ file, even if empty
    aligned_filename = os.path.join(args.od, f"{base_name}_aligned.fastq.gz")
    with gzip.open(aligned_filename, "wt") as aligned_file:
        SeqIO.write(aligned_reads, aligned_file, "fastq")
    print(f"Aligned reads (including potentially empty file) written to {aligned_filename}")

    out_filename = os.path.join(args.od, f"{base_name}_errors.out")
    outfile = open(out_filename, "wt")
    for key in bc_err_pair:
        outfile.write("%s\t" % key)
        for err in bc_err_pair[key]:
            outfile.write("%s\t" % err)
        outfile.write("\n")
    outfile.close()

    
if __name__ == "__main__":
    main()

