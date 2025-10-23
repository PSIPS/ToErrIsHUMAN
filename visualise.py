import py3Dmol

# Protein IDs and paths
proteins = {
    "1lyz": r"D:\CASP_length_project\predictions\1lyz.pdb",
    "5pti": r"D:\CASP_length_project\predictions\5pti.pdb",
    "1ubq": r"D:\CASP_length_project\predictions\1ubq.pdb",
    "3i3z": r"D:\CASP_length_project\predictions\3i3z.pdb",
}

for name, path in proteins.items():
    with open(path) as f:
        pdb_str = f.read()
    view = py3Dmol.view(width=400, height=400)
    view.addModel(pdb_str, 'pdb')
    view.setStyle({'cartoon': {'color':'spectrum'}})
    view.zoomTo()
    print(f"{name}")
    view.show()