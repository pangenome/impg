#!/usr/bin/env python3
import argparse

def adjust_coordinates(index_file, paf_file, output_file):
    # Dictionary to store both start and end values for each sequence
    seq_ranges = {}
    
    # Parse the index file to get ranges
    with open(index_file, "r") as idx:
        for line in idx:
            parts = line.strip().split()
            seq_name, range_str = parts[0].split(':')
            start, end = range_str.split('-')
            # Store both start and end as a tuple
            seq_ranges[seq_name] = (int(start), int(end))
    
    # Process the PAF file
    with open(paf_file, "r") as inpaf, open(output_file, "w") as outpaf:
        for line in inpaf:
            parts = line.strip().split("\t")
            if len(parts) >= 12:
                # Get query sequence name and coordinates
                qname = parts[0]
                qstart = int(parts[2])
                qend = int(parts[3])
                
                # Get target sequence name and coordinates
                tname = parts[5]
                tstart = int(parts[7])
                tend = int(parts[8])
                
                # Adjust coordinates if ranges are available
                if qname in seq_ranges:
                    start_pos, end_pos = seq_ranges[qname]
                    adjusted_qstart = qstart - start_pos
                    adjusted_qend = qend - start_pos
                    parts[2] = str(adjusted_qstart)
                    parts[3] = str(adjusted_qend)
                    # Replace query name with old_name:start-end format
                    parts[0] = f"{qname}:{start_pos}-{end_pos}"
                
                if tname in seq_ranges:
                    start_pos, end_pos = seq_ranges[tname]
                    adjusted_tstart = tstart - start_pos
                    adjusted_tend = tend - start_pos
                    parts[7] = str(adjusted_tstart)
                    parts[8] = str(adjusted_tend)
                    # Replace target name with old_name:start-end format
                    parts[5] = f"{tname}:{start_pos}-{end_pos}"
                
                # Write adjusted line to output file
                outpaf.write("\t".join(parts) + "\n")

def main():
    parser = argparse.ArgumentParser(description='Adjust PAF coordinates based on FASTA index')
    parser.add_argument('-i', '--index', required=True, help='Path to the FASTA index file')
    parser.add_argument('-p', '--paf', required=True, help='Path to the input PAF file')
    parser.add_argument('-o', '--output', default='adjusted.paf', help='Path to the output PAF file (default: adjusted.paf)')
    
    args = parser.parse_args()
    
    adjust_coordinates(args.index, args.paf, args.output)
    print(f"Coordinates adjusted and names replaced, saved to {args.output}")

if __name__ == "__main__":
    main()
