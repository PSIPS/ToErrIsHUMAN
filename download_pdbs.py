from Bio.PDB import PDBList

pdb_ids = ["1LYZ", "5PTI", "1UBQ", "3I3Z"]

pdbl = PDBList()
for pdb_id in pdb_ids:
    pdbl.retrieve_pdb_file(pdb_id, pdir="native", file_format="pdb")
    print(f"Downloaded {pdb_id}")
