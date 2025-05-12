import argparse
from Bio.PDB import PDBParser, PDBIO

def trim_pdb(input_pdb, output_pdb, start_residue, end_residue):
    # Load the PDB file
    parser = PDBParser()
    structure = parser.get_structure('structure', input_pdb)

    # Loop through models and chains to keep only the desired residues
    for model in structure:
        for chain in model:
            residues_to_keep = []
            for residue in chain:
                if start_residue <= residue.id[1] <= end_residue:
                    residues_to_keep.append(residue)
            chain.child_list = residues_to_keep  # Replace chain's residues with filtered list

    # Save the trimmed structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)  # Set the structure object
    io.save(output_pdb)  # Save to output file
    print(f"Trimmed structure saved to {output_pdb}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Trim a PDB structure to specific residue range.')
    parser.add_argument('-i', '--input', required=True, help='Input PDB file')
    parser.add_argument('-o', '--output', required=True, help='Output PDB file')
    parser.add_argument('-s', '--start', type=int, required=True, help='Starting residue number')
    parser.add_argument('-e', '--end', type=int, required=True, help='Ending residue number')

    # Parse arguments
    args = parser.parse_args()

    # Run the trimming function with the provided arguments
    trim_pdb(args.input, args.output, args.start, args.end)

