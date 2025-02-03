import re
import argparse

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


def write_gfa(output_path, paths, sequences):
    submer_size = len(next(iter(sequences.values())))
    print("submer_size",  submer_size)

    with open(output_path, 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        
        for node_id, sequence in sequences.items():
            f.write(f"S\t{node_id}\t{sequence}\n")
        
        for path in paths:
            path_nodes = path["sequence"]
            path_offsets = path["offsets"]
            
            # Generate L lines for consecutive pairs of nodes
            for i in range(len(path_nodes) - 1):
                from_node = path_nodes[i]
                to_node = path_nodes[i + 1]
                from_offset = path_offsets[i]
                to_offset = path_offsets[i + 1]
                offset = (int(from_offset) + submer_size - int(to_offset))
                f.write(f"L\t{from_node[:-1]}\t{from_node[-1]}\t{to_node[:-1]}\t{to_node[-1]}\t{offset}M\n")

        for path in paths:
            path_nodes = path["sequence"]
            path_offsets = path["offsets"]
            offsets = []
             # Generate L lines for consecutive pairs of nodes
            for i in range(len(path_nodes) - 1):
                from_offset = path_offsets[i]
                to_offset = path_offsets[i + 1]
                offset = (int(from_offset) + submer_size - int(to_offset))
                offsets.append(str(offset)+"M")
            f.write(f"P\t{path['id']}\t{','.join(path_nodes)}\t{','.join(offsets)}\n")


def main():
    parser = argparse.ArgumentParser(description="Process path, and syncmer FASTA files to generate a GFA file. \n Requires .txt files (generated with e.g. ONEview -h .1path > .path)")
    parser.add_argument("path_file", help="Path to the path file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_gfa", help="Path to the output GFA file")
    
    args = parser.parse_args()

    # TODO: parse out directly from ONEfile using Python API (instead of converting to text; 
    # currently need to take in .txt files which are the result of e.g. ONEview graph.gbwt > graph.gbwt.txt)
    
    #####################
    # e.g.
    # import ONEcode
    # schema = ONEcode.ONEschema(open("...").read())
    # Read schema from a file
    # Construct a ONEfile object
    # onefile = ONEcode.ONEfile("...", "r", schema, "...", 1)
    # Various methods
    # onefile.readLine()
    #####################

    sequences = parse_fasta(args.fasta_file)
    print("parsed seq")
    
    paths = parse_paths(args.path_file)
    print("parsed paths", paths)
    write_gfa(args.output_gfa, paths, sequences)
    print(f"GFA file written to {args.output_gfa}")

if __name__ == "__main__":
    main()
