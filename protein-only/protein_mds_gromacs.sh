#!/bin/bash
# =============================================================================
# Protein-Only Molecular Dynamics Simulation Pipeline
# GROMACS + CHARMM36 Force Field + SPC/E Water
# =============================================================================
# Requirements: GROMACS, xmgrace
# Input: protein.pdb (cleaned structure — no ligand, no solvent, no HETATM)
# =============================================================================

# =============================================================================
# STEP 1: Load GROMACS Environment
# =============================================================================

source /usr/local/gromacs/bin/GMXRC

# =============================================================================
# STEP 2: Generate Topology
# =============================================================================

gmx pdb2gmx \
  -f protein.pdb \
  -o protein.gro \
  -water spce \
  -ignh
# Interactive prompt: Select force field — choose 1 (CHARMM36)

# Verify topol.top ends with:
#   [ molecules ]
#   ; Compound        #mols
#   Protein_chain_A     1

# =============================================================================
# STEP 3: Define Simulation Box & Solvate
# =============================================================================

gmx editconf \
  -f protein.gro \
  -o newbox.gro \
  -c \
  -d 1.0 \
  -bt dodecahedron

gmx solvate \
  -cp newbox.gro \
  -cs spc216.gro \
  -o solv.gro \
  -p topol.top

# Verify topol.top now ends with:
#   [ molecules ]
#   ; Compound        #mols
#   Protein_chain_A     1
#   SOL             XXXXX

# =============================================================================
# STEP 4: Add Ions
# =============================================================================

# Download ions.mdp from http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
gmx grompp \
  -f mdp_files/ions.mdp \
  -c solv.gro \
  -p topol.top \
  -o ions.tpr

gmx genion \
-s ions.tpr \
-o solv_ions.gro \
-p topol.top \
-pname NA \
-nname CL \
-conc 0.15
# Interactive prompt: Select group 13 (SOL)

# If system is not neutralizing (NA and CL counts are equal), use this to force neutralization at 0.15 M:
# echo "SOL" | gmx genion \
#   -s ions.tpr \
#   -o solv_ions.gro \
#   -p topol.top \
#   -pname NA \
#   -nname CL \
#   -neutral \
#   -conc 0.15

# Verify topol.top now ends with:
#   [ molecules ]
#   ; Compound        #mols
#   Protein_chain_A     1
#   SOL             12588
#   NA                  8   (or CL, depending on protein charge)

# =============================================================================
# STEP 5: Energy Minimization
# =============================================================================

# Download minim.mdp from http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
gmx grompp \
  -f mdp_files/minim.mdp \
  -c solv_ions.gro \
  -p topol.top \
  -o em.tpr

gmx mdrun -v -deffnm em

# =============================================================================
# STEP 6: Equilibration — NVT (Constant Volume & Temperature)
# =============================================================================

# Download nvt.mdp from http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
gmx grompp \                  
  -f mdp_files/nvt.mdp \
  -c em.gro \
  -r em.gro \
  -p topol.top \
  -o nvt.tpr

gmx mdrun -v -deffnm nvt

# =============================================================================
# STEP 7: Equilibration — NPT (Constant Pressure & Temperature)
# =============================================================================

# Download npt.mdp from http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
gmx grompp \
  -f mdp_files/npt.mdp \
  -c nvt.gro \
  -r nvt.gro \
  -t nvt.cpt \
  -p topol.top \
  -o npt.tpr

gmx mdrun -v -deffnm npt

# =============================================================================
# STEP 8: Production MD Run
# =============================================================================

# Download npt.mdp from http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
# Edit md.mdp: set nsteps to desired simulation length
#   e.g., nsteps = 5000000 → 10 ns at 2 fs/step
gmx grompp \
  -f mdp_files/md.mdp \
  -c npt.gro \
  -t npt.cpt \
  -p topol.top \
  -o md.tpr

gmx mdrun -v -deffnm md

# To resume a stopped or crashed run:
# gmx mdrun -s md.tpr -cpi md.cpt -deffnm md -append -v

# =============================================================================
# STEP 9: Trajectory Post-Processing
# =============================================================================

# Remove periodic boundary condition artifacts and center the protein
gmx trjconv \
  -s md.tpr \
  -f md.xtc \
  -o md_noPBC.xtc \
  -pbc mol \
  -center
# Interactive prompt:
#   Select "Protein" (1) for centering
#   Select "System"  (0) for output

# =============================================================================
# STEP 10: RMSD Analysis
# =============================================================================

gmx rms \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsd.xvg \
  -tu ns
# Interactive prompt:
#   Select "Backbone" (4) for least-squares fitting
#   Select "Backbone" (4) for RMSD calculation

xmgrace rmsd.xvg

# =============================================================================
# STEP 11: RMSF Analysis
# =============================================================================

gmx rmsf \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsf.xvg \
  -res
# Interactive prompt: Select "Backbone" (4)

xmgrace rmsf.xvg

# =============================================================================
# STEP 12: Radius of Gyration (Rg)
# =============================================================================

gmx gyrate \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rg.xvg \
  -sel Protein \
  -tu ns

xmgrace rg.xvg

# =============================================================================
# STEP 13: Solvent Accessible Surface Area (SASA)
# =============================================================================

gmx sasa \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o sasa.xvg
# Interactive prompt: Select "Protein" (1)

xmgrace sasa.xvg

# =============================================================================
# STEP 14: Hydrogen Bond Analysis
# =============================================================================

gmx hbond \
  -s md.tpr \
  -f md_noPBC.xtc \
  -tu ns \
  -num hbnum.xvg
# Interactive prompt:
#   Select "MainChain+H" (7) for donor group
#   Select "MainChain+H" (7) for acceptor group

xmgrace hbnum.xvg

# =============================================================================
# END OF PIPELINE
# =============================================================================
