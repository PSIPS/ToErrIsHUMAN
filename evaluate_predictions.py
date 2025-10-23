import os
import subprocess
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
import matplotlib.pyplot as plt

# --- SETTINGS ---
PDB_IDS = ["1lyz", "5pti", "1ubq", "3i3z"]

EXP_DIR = "native"          # your .ent files
PRED_DIR = "predictions"    # AlphaFold predictions
OUT_DIR = "results"
os.makedirs(OUT_DIR, exist_ok=True)

parser = PDBParser(QUIET=True)
rows = []

def get_ca_atoms(struct):
    return [a for a in struct.get_atoms() if a.get_id() == "CA"]

def compute_rmsd(exp_path, pred_path):
    ref = parser.get_structure("exp", exp_path)
    mob = parser.get_structure("pred", pred_path)
    ref_ca = get_ca_atoms(ref)
    mob_ca = get_ca_atoms(mob)
    n = min(len(ref_ca), len(mob_ca))
    sup = Superimposer()
    sup.set_atoms(ref_ca[:n], mob_ca[:n])
    return sup.rms, len(ref_ca)

def run_tmalign(exp_path, pred_path):
    try:
        result = subprocess.run(
            ["TMalign", pred_path, exp_path],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            if line.strip().startswith("TM-score="):
                return float(line.split()[1])
    except FileNotFoundError:
        return None
    except Exception:
        return None
    return None

def mean_plddt(pred_path):
    """Compute mean pLDDT from AlphaFold PDB (B-factor column)"""
    pl_values = []
    with open(pred_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    bfactor = float(line[60:66])
                    pl_values.append(bfactor)
                except ValueError:
                    continue
    return sum(pl_values)/len(pl_values) if pl_values else None

# --- MAIN LOOP ---
for pid in PDB_IDS:
    exp_path = os.path.join(EXP_DIR, f"{pid}.ent")   # lowercase .ent
    pred_path = os.path.join(PRED_DIR, f"{pid}.pdb") # lowercase .pdb
    if not os.path.exists(exp_path) or not os.path.exists(pred_path):
        print(f"Skipping {pid} (file missing)")
        continue

    rmsd, length = compute_rmsd(exp_path, pred_path)
    tm = run_tmalign(exp_path, pred_path)
    plddt = mean_plddt(pred_path)

    rows.append({
        "pdb_id": pid,
        "length": length,
        "rmsd": rmsd,
        "tm_score": tm,
        "mean_plddt": plddt
    })

    print(f"{pid}: RMSD={rmsd:.2f} Å  TM={tm}  mean_pLDDT={plddt:.2f}")

# --- SAVE RESULTS ---
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "results.csv"), index=False)
print("\nSaved results.csv")

# --- PLOT RMSD vs LENGTH ---
if not df.empty:
    plt.scatter(df["length"], df["rmsd"], s=60)
    plt.xlabel("Sequence length (aa)")
    plt.ylabel("RMSD (Å)")
    plt.title("Prediction accuracy vs. sequence length")
    plt.grid(True)
    plt.savefig(os.path.join(OUT_DIR, "rmsd_vs_length.png"), dpi=300)
    plt.show()
else:
    print("No data to plot.")



