from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolToMolBlock
import subprocess
import glob
import os
import tempfile
from pathlib import Path
import pandas as pd


#quickly put everything into a class...

class InitializeDirectories():
    def __init__(self):
        Path('libgen').mkdir()
        Path('fragment_mol2').mkdir()
        


class SubstituentCreate():
    #Path(f"{os.getcwd()}") / "multixyz.xyz"
    
    def __init__(self, substituentfile = "substituents.csv",
                 base_structs = ["ClC#Cc1c(X)cccc1", "ClC#Cc1cc(X)ccc1", "ClC#Cc1ccc(X)cc1"],
                 sub_list = ["OY", "NY", "FY", "ClY", "YCCCC","YC(F)(F)F", "Y[N+](=O)[O-]", "c2ccccc2", "C2CCCC2"]
                 ):
        self.base_struct = ["ClC#Cc1c(X)cccc1", "ClC#Cc1cc(X)ccc1", "ClC#Cc1ccc(X)cc1"]
        self.substituentfile = Path(f"{os.getcwd()}")/substituentfile
        try:
            self.substituentfile.open()
            self.sub_list = pd.read_csv(SUBLIST_FILE).to_numpy().squeeze()
        except:
            print("no substituent list file found, please place in working directory, reverting to default list of substituents")
            self.sub_list = sub_list
        
        self.library = []
        self.mol2lib = []
    

    def make_substituents(self):
        for base in self.base_struct:
            # Iterate through each substituent and generate the combined structure
            for substituent in self.sub_list:
                # Replace "X" in the base structure with the substituent
                combined_structure = base.replace("X", substituent).replace("Y", "")
                
                # Use RDKit to convert the SMILES string to a molecule object
                mol = Chem.MolFromSmiles(combined_structure)
            
                if mol is not None:
                    # Generate 3D coordinates
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500, nonBondedThresh=200.0)
                    
                    # Extract XYZ coordinates
                    conf = mol.GetConformer(0)
                    #print(atom.GetSymbol())
                    xyz_coordinates = [(atom.GetSymbol(), *tuple(conf.GetAtomPosition(i))) for i, atom in enumerate(mol.GetAtoms())]
                    mol.SetProp("_Name", f"{base}_{substituent}")
                    mol2_block = MolToMolBlock(mol, includeStereo=True, kekulize=True)
                    self.library.append((mol.GetNumAtoms(), combined_structure, xyz_coordinates))\      
                    self.mol2lib.append(mol2_block)
                    
                else:
                    print(f"Invalid SMILES: {self.combined_structure}")        
    
    def write_xyzs(self):
        Path('libgen').mkdir(exist_ok=True)
        #Path('fragment_mol2').mkdir(exist_ok=True)
        Path('fragment_mol2/xyzs').mkdir(parents=True, exist_ok=True)
        
        multixyz = Path(f"{os.getcwd()}") / "multixyz.xyz"
        with multixyz.open(mode="w+") as f:
            for atoms, smiles, xyz in self.library:
                f.write(f"{atoms}\n{smiles}\n{pretty_print_xyz(xyz)}\n")
                print("writing multixyz")
        
        # Lets print all XYZ coordinates to many files and then transform them with obabel into mol2 files
        
        
        for atoms, smiles, xyz in self.library:
            _xyztmp = Path(f"{os.getcwd()}") / f"{smiles.replace('(', '_').replace(')', '_').replace('[', '_').replace(']', '_')}_tmp.xyz"
            with _xyztmp.open(mode="w+") as f:
                print(atoms)
                f.write(f"{atoms}\n{smiles}\n{pretty_print_xyz(xyz)}\n")
                print("writing")
        
        #convert to sybyl mol2 files and erase previous .xyz files
        for i, f in enumerate(glob.glob("*_tmp.xyz")):
            outputfile = Path(f"{os.getcwd()}/fragment_mol2") / f"obabel_\"{f}\"_to.mol2"
            obabel_command = f"obabel \"{f}\" -O {outputfile}"
            subprocess.run(obabel_command, shell=True, check=True)
            #subprocess.run(f"mv \"{f}\"", shell=True, check=True)


    def pretty_print_xyz(cls, xyz_coordinates):
        lines = [f'{symbol} {x:.4f} {y:.4f} {z:.4f}' for symbol, x, y, z in xyz_coordinates]
        return '\n'.join(lines)

    

