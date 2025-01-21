import re
import argparse

def parse_gbwt(file_path, submer_size):
    segments = []
    links = []
    node_id = 0
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if line.startswith("V"):
                node_id += 1
                segments.append({"id": node_id})
            elif line.startswith("E") or line.startswith("e"):
                reverse = int(parts[1]) < 0
                parts[1] = str(abs(int(parts[1])))
                from_node, to_node = (parts[1], node_id) if line.startswith("E") else (node_id, parts[1])
                from_direction = "-" if reverse and line.startswith("E") else "+"
                to_direction = "-" if reverse and line.startswith("e") else "+"
                overlap = int(parts[2])
                links.append({"from": from_node, "to": to_node, "from_direction": from_direction, "to_direction": to_direction, "overlap": submer_size - overlap})
    return segments, links

def parse_paths(file_path):
    paths = []
    current_path = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("P"):
                current_path += 1
                path_id = current_path
                path_sequence = []
                offsets = []
            elif line.startswith("z"):
                parts = line.strip().split()
                path_sequence = [f"{part[1:]}-" if part.startswith("-") else f"{part}+" for part in parts[2:]]
            elif line.startswith("o"):
                parts = line.strip().split()
                offsets = parts[2:]
                paths.append({"id": path_id, "sequence": path_sequence, "offsets": offsets})
    return paths

def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        current_id = None
        for line in f:
            if line.startswith(">"):
                current_id = int(line[1:].strip())
            else:
                sequences[current_id] = line.strip()
    return sequences

def write_gfa(output_path, segments, links, paths, sequences):
    with open(output_path, 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        for seg in segments:
            sequence = sequences.get(seg["id"], "*")
            f.write(f"S\t{seg['id']}\t{sequence}\n")
        for link in links:
            f.write(f"L\t{link['from']}\t{link['from_direction']}\t{link['to']}\t{link['to_direction']}\t{link['overlap']}M\n")
        for path in paths:
            path_nodes = ",".join(path["sequence"])
            path_offsets = ",".join(path["offsets"])
            f.write(f"P\t{path['id']}\t{path_nodes}\t{path_offsets}\n")

def main():
    parser = argparse.ArgumentParser(description="Process GBWT, path, and FASTA files to generate a GFA file.")
    parser.add_argument("gbwt_file", help="Path to the GBWT file")
    parser.add_argument("path_file", help="Path to the path file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_gfa", help="Path to the output GFA file")
    
    args = parser.parse_args()

    sequences = parse_fasta(args.fasta_file)
    print("parsed seq")
    submer_size = len(next(iter(sequences.values())))
    print("submer_size",  submer_size)
    segments, links = parse_gbwt(args.gbwt_file, submer_size)
    print("parsed gbwt")
    paths = parse_paths(args.path_file)
    print("parsed paths")
    write_gfa(args.output_gfa, segments, links, paths, sequences)
    print(f"GFA file written to {args.output_gfa}")

if __name__ == "__main__":
    main()
