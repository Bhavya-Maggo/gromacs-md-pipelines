# GROMACS MD Simulation Pipelines
A collection of step-by-step molecular dynamics simulation pipelines built with **GROMACS**, **CHARMM36** force field. Covers three simulation types — protein-only, protein–ligand complex, and protein–protein complex.

---

## Pipelines

### 1. `protein-only/`
MD simulation of a single protein in explicit solvent.

- Force field: CHARMM36
- Water model: SPC/E
- Analysis: RMSD, RMSF, Radius of Gyration, SASA, H-bonds

**MDP files reference:** [GROMACS Lysozyme Tutorial](http://www.mdtutorials.com/gmx/lysozyme/index.html)

---

### 2. `protein-ligand/`
MD simulation of a protein–ligand complex with full ligand parameterization via CGenFF.

- Force field: CHARMM36
- Water model: TIP3P
- Ligand parameterization: CGenFF / `cgenff_charmm2gmx.py`
- Analysis: RMSD (protein, ligand, complex), RMSF, Rg, SASA, H-bonds (number, distance, angle)

**MDP files reference:** [GROMACS Protein–Ligand Complex Tutorial](http://www.mdtutorials.com/gmx/complex/index.html)

---

### 3. `protein-protein/`
MD simulation of a two-chain protein–protein complex with per-chain and interface analysis.

- Force field: CHARMM36
- Water model: SPC/E
- Custom index groups: Chain A, Chain B, combined complex
- Analysis: RMSD (complex, chain A, chain B), RMSF, Rg, SASA, interface H-bonds

**MDP files:** Included in `protein-protein/mdp_files/`

---

## Repository Structure

```
gromacs-md-pipelines/
├── README.md
├── protein-only/
│   ├── md_pipeline.sh
│   └── mdp_files/        ← see lysozyme tutorial link above
├── protein-ligand/
│   ├── md_pipeline.sh
│   └── mdp_files/        ← see protein-ligand tutorial link above
└── protein-protein/
    ├── md_pipeline.sh
    └── mdp_files/
        ├── ions.mdp
        ├── minim.mdp
        ├── nvt.mdp
        ├── npt.mdp
        └── md.mdp
```

---

## Requirements

- [GROMACS](https://www.gromacs.org/) 2025 or later
- [OpenBabel](https://openbabel.org/) 
- [xmgrace](https://plasma-gate.weizmann.ac.il/Grace/) for `.xvg` visualization
- [CGenFF Server](https://cgenff.com/) *(protein–ligand only)*
- Perl + `sort_mol2_bonds.pl` *(protein–ligand only)*
- Python + `cgenff_charmm2gmx.py` *(protein–ligand only)*

---

## General Workflow

Every pipeline follows the same backbone:

```
Input PDB
   └─► Topology (pdb2gmx)
          └─► Solvation (editconf + solvate)
                 └─► Ion Addition (grompp + genion)
                        └─► Energy Minimization
                               └─► NVT Equilibration
                                      └─► NPT Equilibration
                                             └─► Production MD
                                                    └─► Analysis
```

---

## References

- [GROMACS Documentation](https://manual.gromacs.org/)
- [GROMACS Lysozyme Tutorial](http://www.mdtutorials.com/gmx/lysozyme/index.html)
- [GROMACS Protein–Ligand Tutorial](http://www.mdtutorials.com/gmx/complex/index.html)
- [CGenFF Server](https://cgenff.com/)
- Vanommeslaeghe, K., & MacKerell Jr, A. D. (2012). Automation of the CHARMM General Force Field (CGenFF) I: bond perception and atom typing. Journal of chemical information and modeling, 52(12), 3144-3154.

---

