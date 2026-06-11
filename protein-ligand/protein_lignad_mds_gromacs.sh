#!/bin/bash
# =============================================================================
# Protein-Ligand Molecular Dynamics Simulation Pipeline
# GROMACS + CHARMM36 Force Field
# =============================================================================
# Requirements: GROMACS, OpenBabel, Perl (sort_mol2_bonds.pl), Python (cgenff_charmm2gmx.py)
# =============================================================================

# =============================================================================
# STEP 1: Ligand Optimization Before Docking
# =============================================================================

obabel ligand.sdf -O ligand_final.sdf \
  --gen3d \
  --ph 7.4 \
  --minimize \
  --ff MMFF94 \
  --steps 2500

# Convert optimized ligand to PDB format
obabel ligand_final.sdf -O ligand.pdb

# =============================================================================
# STEP 2: Load GROMACS Environment
# =============================================================================

source /usr/local/gromacs/bin/GMXRC

# =============================================================================
# STEP 3: Protein Preparation
# =============================================================================

gmx pdb2gmx \
  -f protein.pdb \
  -o protein.gro \
  -p protein.top \
  -water tip3p \
  -ignh

# =============================================================================
# STEP 4: Ligand Parameterization
# =============================================================================

# -- Step 4.1: Convert SDF to MOL2 --
# Use Discovery Studio:
#   a) Convert ligand .sdf to Sybyl .mol2
#   b) Add hydrogen atoms via Chemistry tools

# -- Step 4.2: Assign Partial Charges (Gasteiger) --
obabel lig.mol2 -O lig_charge.mol2 --partialcharge gasteiger

# -- Step 4.3: Rename Ligand Residue --
# Replace UNL11 with LIG
sed -i 's/UNL11/LIG/g' lig_charge.mol2

# -- Step 4.4: Sort Bond Order --
perl sort_mol2_bonds.pl lig_charge.mol2 lig_fix.mol2
# IMPORTANT: In lig_fix.mol2, update the ligand name on the 2nd line after the @<TRIPOS>MOLECULE header to "LIG"
# Download [sort_mol2_bonds.pl] (http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.txt)

# -- Step 4.5: Generate CGenFF Parameters --
# Upload lig_fix.mol2 to https://cgenff.com/
# Download the resulting .str file

# -- Step 4.6: Convert CGenFF Files for GROMACS --
python cgenff_charmm2gmx.py lig lig_fix.mol2 lig_fix.str charmm36-jul2022.ff
# Generates: lig_ini.pdb, lig.top, lig.itp, lig.prm

# -- Step 4.7: Convert Ligand Structure to GROMACS Format --
gmx editconf -f lig_ini.pdb -o lig.gro

# =============================================================================
# STEP 5: Build Protein–Ligand Complex
# =============================================================================

# Manually merge protein.gro and lig.gro into complex.gro:
#   - Append ligand coordinates below protein coordinates
#   - Update the total atom count at the top of complex.gro

# =============================================================================
# STEP 6: Modify Topology File (topol.top)
# =============================================================================

# Key additions:
#   - #include "lig.prm"  (after forcefield include, before [ moleculetype ])
#   - #include "lig.itp"  (after protein position restraints)
#   - #ifdef POSRES_LIG block for ligand restraints
#   - Add "lig  1" to the [ molecules ] section

# =============================================================================
# STEP 7: Solvation
# =============================================================================

gmx editconf \
  -f complex.gro \
  -o newbox.gro \
  -bt dodecahedron \
  -d 1.0

gmx solvate \
  -cp newbox.gro \
  -cs spc216.gro \
  -p topol.top \
  -o solv.gro

# Verify: Open topol.top — SOL should appear in the [ molecules ] section

# =============================================================================
# STEP 8: Add Ions
# =============================================================================

# Download ions.mdp from http://www.mdtutorials.com/gmx/complex/Files/ions.mdp
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

gmx genion \
  -s ions.tpr \
  -o solv_ions.gro \
  -p topol.top \
  -pname NA \
  -nname CL \
  -neutral \
  -conc 0.15
# Interactive prompt: Select SOL (group 15)

# Verify: Check topol.top — NA and CL should appear in [ molecules ]

# =============================================================================
# STEP 9: Energy Minimization
# =============================================================================

# Download em.mdp from http://www.mdtutorials.com/gmx/complex/Files/em.mdp
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# =============================================================================
# STEP 10: Equilibration — NVT (Constant Volume & Temperature)
# =============================================================================

# Create ligand index group (heavy atoms only)
gmx make_ndx -f lig.gro -o index_lig.ndx
# At prompt, type:
#   0 & ! a H*
#   q

# Generate ligand position restraints
gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
# Interactive prompt: Select System_&_!H* (group 3)

# Create combined Protein+Ligand index
gmx make_ndx -f em.gro -o index.ndx
# At prompt, type:
#   1 | 13      (1 = Protein, 13 = LIG)
#   q

# Download nvt.mdp from http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp
# Edit nvt.mdp: change tc-grps
#   FROM: tc-grps = Protein_JZ4 Water_and_ions
#   TO:   tc-grps = Protein LIG Water_and_ions

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -v -deffnm nvt

# =============================================================================
# STEP 11: Equilibration — NPT (Constant Pressure & Temperature)
# =============================================================================

# Download npt.mdp from http://www.mdtutorials.com/gmx/complex/Files/npt.mdp
# Edit npt.mdp: update tc-grps (same change as nvt.mdp above)

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -v -deffnm npt

# =============================================================================
# STEP 12: Production MD Run
# =============================================================================

# Download md.mdp from mdp_files/ directory in this repo
# Edit md.mdp:
#   - Update tc-grps (same change as above)
#   - Set nsteps to desired simulation length (e.g., 5000000 = 10 ns at 2 fs/step)

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -v -deffnm md_0_10

# To resume a stopped/crashed run:
# gmx mdrun -s md_0_10.tpr -cpi md_0_10.cpt -deffnm md_0_10 -append -v

# =============================================================================
# STEP 13: Trajectory Post-Processing
# =============================================================================

# Rebuild index with Protein+Ligand group
gmx make_ndx -f em.gro -o index.ndx
# At prompt:
#   1 | 13		(1 = Protein, 13 = LIG)
#   q

# Center trajectory and fix periodic boundary conditions
gmx trjconv \
  -s md_0_10.tpr \
  -f md_0_10.xtc \
  -o md_0_10_center.xtc \
  -center -pbc mol -ur compact \
  -n index.ndx
# Prompt: Select "Protein_lig" (21) for centering, "System" (0) for output

# Fit trajectory (remove rotation and translation)
gmx trjconv \
  -s md_0_10.tpr \
  -f md_0_10_center.xtc \
  -o md_0_10_fit.xtc \
  -fit rot+trans \
  -n index.ndx
# Prompt: Select "Backbone" (4) for fitting, "System" (0) for output

# Check total number of frames
gmx check -f md_0_10_fit.xtc | grep "Frames"

# =============================================================================
# STEP 14: RMSD Analysis
# =============================================================================

# Add ligand heavy atoms to index
gmx make_ndx -f em.gro -n index.ndx
# At prompt:
#   13 & ! a H*
#   name lig_heavy
#   q

# RMSD — Ligand
gmx rms \
  -s em.tpr \
  -f md_0_10_center.xtc \
  -n index.ndx \
  -tu ns \
  -o rmsd_lig.xvg
# Prompt: Select "lig_&_!H*" (22) for fitting AND RMSD group
xmgrace rmsd_lig.xvg

# RMSD — Protein Backbone
gmx rms \
  -s em.tpr \
  -f md_0_10_center.xtc \
  -n index.ndx \
  -tu ns \
  -o rmsd_protein.xvg
# Prompt: Select "Backbone" (4) for fitting AND RMSD group
xmgrace rmsd_protein.xvg

# RMSD — Protein + Ligand Complex
gmx rms \
  -s em.tpr \
  -f md_0_10_center.xtc \
  -n index.ndx \
  -tu ns \
  -o rmsd_protein_lig.xvg
# Prompt: Select "Protein_lig" (21) for fitting AND RMSD group
xmgrace rmsd_protein_lig.xvg

# =============================================================================
# STEP 15: RMSF Analysis
# =============================================================================

gmx rmsf \
  -s em.tpr \
  -f md_0_10_center.xtc \
  -n index.ndx \
  -o rmsf_backbone.xvg \
  -res
# Prompt: Select "Backbone" (4)
xmgrace rmsf_backbone.xvg

# =============================================================================
# STEP 16: Radius of Gyration (Rg)
# =============================================================================

gmx gyrate \
  -f md_0_10_center.xtc \
  -s md_0_10.tpr \
  -n index.ndx \
  -tu ns \
  -o rg.xvg
# Prompt: Select "Protein_lig" (21)
xmgrace rg.xvg

# =============================================================================
# STEP 17: Hydrogen Bond Analysis
# =============================================================================

# Number of H-bonds over time
gmx hbond \
  -f md_0_10_center.xtc \
  -s md_0_10.tpr \
  -n index.ndx \
  -num hb.xvg \
  -tu ns
# Prompt: Select "lig" (13) as donor/acceptor, "Protein" (1) as reference
xmgrace hb.xvg

# H-bond distances
gmx hbond \
  -f md_0_10_center.xtc \
  -s md_0_10.tpr \
  -n index.ndx \
  -dist hbond_dist.xvg \
  -tu ns
# Prompt: Select "lig" (13) and "Protein" (1)
xmgrace hbond_dist.xvg

# H-bond angles
gmx hbond \
  -f md_0_10_center.xtc \
  -s md_0_10.tpr \
  -n index.ndx \
  -ang hbond_ang.xvg \
  -tu ns
# Prompt: Select "lig" (13) and "Protein" (1)
xmgrace hbond_ang.xvg

# =============================================================================
# STEP 18: Solvent Accessible Surface Area (SASA)
# =============================================================================

gmx sasa \
  -f md_0_10_center.xtc \
  -s md_0_10.tpr \
  -n index.ndx \
  -tu ns \
  -o sasa_total.xvg \
  -surface "Protein_lig"

xmgrace sasa_total.xvg

# =============================================================================
# STEP 19: Energy Analysis
# =============================================================================

gmx energy -f md_0_10.edr -o energy.xvg

# =============================================================================
# END OF PIPELINE
# =============================================================================
