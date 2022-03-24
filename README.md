# TargetID: Drug Target and Chemotype Identification Pipeline

The TargetID pipeline enables rapid identification and characterization of binding sites in SARS-CoV-2 proteins as well as the core chemical components with which these sites interact. 

The Jupyter notebook residue_clustering.ipynb takes protein-ligand binding data from the PDBspheres program and performs a clustering procedure to identify groups of residues that represent relevant “consensus binding pockets.” 

The Jupyter notebook cluster-composition-ncov-residues.ipynb is used to characterize the consensus pockets in terms of the source organisms of the underlying PDB templates that contribute to each pocket. 

The Python script ligand-clustering.py performs another clustering procedure to identify various groups of ligands that bind to each consensus pocket. 

The remaining files are input data files required to run the three notebooks. 

Note that Python packages NetworkX (2.2), scikit-learn (0.20.2), and RDKit (2020.09.1.0) are required. 

LLNL-CODE-831985