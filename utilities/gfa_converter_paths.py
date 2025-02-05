import re
import argparse

def parse_fasta(file_path, nodes_to_avoid):
    sequences = {}
    with open(file_path, 'r') as f:
        current_id = None
        for line in f:
            if line.startswith(">"):
                current_id = int(line[1:].strip())
            else:
                if current_id not in nodes_to_avoid:
                    sequences[current_id] = line.strip()
    return sequences

def parse_paths(file_path, nodes_to_avoid):
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
                filtered_sequence = []
                filtered_offsets = []
                for i, node in enumerate(path_sequence):
                    node_id = int(node[:-1])
                    if node_id not in nodes_to_avoid:
                        filtered_sequence.append(node)
                        filtered_offsets.append(offsets[i])
                paths.append({"id": path_id, "sequence": filtered_sequence, "offsets": filtered_offsets})
    return paths

def write_gfa(output_path, paths, sequences):
    submer_size = len(next(iter(sequences.values())))
    print("submer_size",  submer_size)

    with open(output_path, 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        
        for node_id, sequence in sequences.items():
            f.write(f"S\t{node_id}\t{sequence}\n")
        
        written_links = set()
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
                link = (from_node, to_node, offset)
                if link not in written_links:
                    f.write(f"L\t{from_node[:-1]}\t{from_node[-1]}\t{to_node[:-1]}\t{to_node[-1]}\t{offset}M\n")
                    written_links.add(link)

        for path in paths:
            path_nodes = path["sequence"]
            path_offsets = path["offsets"]
            offsets = []
            for i in range(len(path_nodes) - 1):
                from_offset = path_offsets[i]
                to_offset = path_offsets[i + 1]
                offset = (int(from_offset) + submer_size - int(to_offset))
                offsets.append(str(offset)+"M")
            f.write(f"P\t{path['id']}\t{','.join(path_nodes)}\t{','.join(offsets)}\n")

def read_nodes_to_avoid(file_path):
    nodes_to_avoid = set()
    with open(file_path, 'r') as f:
        for line in f:
            nodes = line.strip().split()
            for node in nodes:
                nodes_to_avoid.add(int(node))
    return nodes_to_avoid

def main():
    parser = argparse.ArgumentParser(description="Process path, and syncmer FASTA files to generate a GFA file. \n Requires .txt files (generated with e.g. ONEview -h .1path > .path)")
    parser.add_argument("path_file", help="Path to the path file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_gfa", help="Path to the output GFA file")
    parser.add_argument("--nodes_to_avoid", help="Path to the file containing list of nodes to avoid", default=None)
    
    args = parser.parse_args()

    nodes_to_avoid = set()
    if args.nodes_to_avoid:
        nodes_to_avoid = read_nodes_to_avoid(args.nodes_to_avoid)

    sequences = parse_fasta(args.fasta_file, nodes_to_avoid)
    print("parsed seq")
    
    paths = parse_paths(args.path_file, nodes_to_avoid)
    write_gfa(args.output_gfa, paths, sequences)
    print(f"GFA file written to {args.output_gfa}")

if __name__ == "__main__":
    main()
