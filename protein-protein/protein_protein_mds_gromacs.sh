#!/bin/bash
# =============================================================================
# Protein–Protein Complex Molecular Dynamics Simulation Pipeline
# GROMACS + CHARMM36 Force Field + SPC/E Water
# =============================================================================
# Requirements: GROMACS, xmgrace
# Input: complex.pdb — both protein chains in a single PDB file
# =============================================================================
#
# NOTE on atom ranges (Step 8):
#   Before running the index-building step, open em.gro and identify:
#     - Last atom of Chain A  → update variable CHAIN_A_END
#     - Last atom of Chain B  → update variable CHAIN_B_END
#
# =============================================================================

# --- Edit these before running Step 8 ---
CHAIN_A_START=1		  # First atom index of Chain A in em.gro
CHAIN_A_END=7990      # Last atom index of Chain A in em.gro
CHAIN_B_START=7991    # First atom index of Chain B in em.gro
CHAIN_B_END=10642     # Last atom index of Chain B in em.gro
# -----------------------------------------

# =============================================================================
# STEP 1: Load GROMACS Environment
# =============================================================================

source /usr/local/gromacs/bin/GMXRC

# =============================================================================
# STEP 2: Generate Topology
# =============================================================================

gmx pdb2gmx \
  -f complex.pdb \
  -o complex.gro \
  -water spce \
  -ignh
# Interactive prompt: Select force field — choose 1 (CHARMM36)

# =============================================================================
# STEP 3: Define Simulation Box
# =============================================================================

gmx editconf \
  -f complex.gro \
  -o boxed.gro \
  -c \
  -d 1.2 \
  -bt dodecahedron

# =============================================================================
# STEP 4: Solvate
# =============================================================================

gmx solvate \
  -cp boxed.gro \
  -cs spc216.gro \
  -o solv.gro \
  -p topol.top

# Verify topol.top [ molecules ] section now includes SOL

# =============================================================================
# STEP 5: Add Ions
# =============================================================================

# Requires: mdp_files/ions.mdp
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
  -neutral
# Interactive prompt: Select group 13 (SOL)

# Verify topol.top [ molecules ] now lists NA/CL after SOL

# =============================================================================
# STEP 6: Energy Minimization
# =============================================================================

# Requires: mdp_files/minim.mdp
gmx grompp \
  -f mdp_files/minim.mdp \
  -c solv_ions.gro \
  -p topol.top \
  -o em.tpr

gmx mdrun -v -deffnm em

# =============================================================================
# STEP 7: Equilibration — NVT
# =============================================================================

# Requires: mdp_files/nvt.mdp
gmx grompp \
  -f mdp_files/nvt.mdp \
  -c em.gro \
  -r em.gro \
  -p topol.top \
  -o nvt.tpr

gmx mdrun -v -deffnm nvt

# =============================================================================
# STEP 8: Equilibration — NPT
# =============================================================================

# Requires: mdp_files/npt.mdp
gmx grompp \
  -f mdp_files/npt.mdp \
  -c nvt.gro \
  -r nvt.gro \
  -t nvt.cpt \
  -p topol.top \
  -o npt.tpr

gmx mdrun -v -deffnm npt

# =============================================================================
# STEP 9: Production MD Run
# =============================================================================

# Requires: mdp_files/md.mdp
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
# STEP 10: Build Custom Index Groups
# =============================================================================

# ── 10a: Define Chain A and Chain B groups ──
# IMPORTANT: Check em.gro for correct atom ranges and update
#            CHAIN_A_END, CHAIN_B_START, CHAIN_B_END at the top of this script.

gmx make_ndx -f em.gro -o index.ndx << EOF
a ${CHAIN_A_START}-${CHAIN_A_END}
name 17 Protein_A
a ${CHAIN_B_START}-${CHAIN_B_END}
name 18 Protein_B
q
EOF

# ── 10b: Add combined Protein_A_Protein_B group ──
gmx make_ndx -f em.gro -o index.ndx -n index.ndx << EOF
17 | 18
q
EOF
# This creates group 19: Protein_A_Protein_B

# =============================================================================
# STEP 11: Trajectory Post-Processing
# =============================================================================

gmx trjconv \
  -s md.tpr \
  -f md.xtc \
  -o md_noPBC.xtc \
  -pbc mol \
  -center \
  -n index.ndx
# Interactive prompt:
#   Select "Protein_A_Protein_B" (19) for centering
#   Select "System" (0) for output

# =============================================================================
# STEP 12: RMSD Analysis
# =============================================================================

# RMSD — Full Complex
gmx rms \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsd_complex.xvg \
  -tu ns \
  -n index.ndx
# Interactive prompt: Select "Protein_A_Protein_B" (19) for both fit and RMSD
xmgrace rmsd_complex.xvg

# RMSD — Chain A only
gmx rms \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsd_chainA.xvg \
  -tu ns \
  -n index.ndx
# Interactive prompt: Select "Protein_A" (17) for both fit and RMSD
xmgrace rmsd_chainA.xvg

# RMSD — Chain B only
gmx rms \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsd_chainB.xvg \
  -tu ns \
  -n index.ndx
# Interactive prompt: Select "Protein_B" (18) for both fit and RMSD
xmgrace rmsd_chainB.xvg

# =============================================================================
# STEP 13: RMSF Analysis
# =============================================================================

# Add backbone atoms group for complex
gmx make_ndx -f em.gro -o index.ndx -n index.ndx << EOF
19 & a C | a CA | a N
name 22 BB_Complex
q
EOF
# Creates group 20: Protein_A_Protein_B_&_C_CA_N (backbone of complex)

gmx rmsf \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rmsf_complex.xvg \
  -res \
  -n index.ndx
# Interactive prompt: Select group 20 (Protein_A_Protein_B_&_C_CA_N)
xmgrace rmsf_complex.xvg

# =============================================================================
# STEP 14: Radius of Gyration (Rg)
# =============================================================================

gmx gyrate \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o rg_complex.xvg \
  -tu ns \
  -n index.ndx
# Interactive prompt: Select "Protein_A_Protein_B" (19)
xmgrace rg_complex.xvg

# =============================================================================
# STEP 15: Solvent Accessible Surface Area (SASA)
# =============================================================================

gmx sasa \
  -s md.tpr \
  -f md_noPBC.xtc \
  -o sasa_complex.xvg \
  -tu ns \
  -n index.ndx
# Interactive prompt: Select "Protein_A_Protein_B" (19)
xmgrace sasa_complex.xvg

# =============================================================================
# STEP 16: Protein–Protein Hydrogen Bond Analysis
# =============================================================================

gmx hbond \
  -s md.tpr \
  -f md_noPBC.xtc \
  -n index.ndx \
  -tu ns \
  -num hbonds_num.xvg
# Interactive prompt:
#   Select "Protein_A" (17) as donor/acceptor group 'r'
#   Select "Protein_B" (18) as reference group 't'
xmgrace hbonds_num.xvg

# =============================================================================
# END OF PIPELINE
# =============================================================================
