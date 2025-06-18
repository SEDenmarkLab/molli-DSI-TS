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

default_substituents = ["YOC",
                        "YOCC",
                        "YOC(C)C",
                        "YOC(C)(C)C",
                        "YOC(c2ccccc2)",
                        "YF",
                        "YCl",
                        "YBr",
                        "YC",
                        "Yc2ccccc2",
                        "YCC",
                        "YC(C)(C)",
                        "YC(C)(C)C",
                        "Y[Si](C)(C)C",
                        "Y[Si](C)(C(C)(C)C)C",
                        "YC(F)(F)F",
                        "Y[N+](=O)[O-]",
                        "YC(c2ccccc2)",
                        "YCC(c2ccccc2)",
                        "Y(c2ccccc2)"]

base_structures = ["ClC#Cc1c(X)cccc1",
                   "ClC#Cc1cc(X)ccc1",
                   "ClC#Cc1ccc(X)cc1"]

class SubstituentCreate():
    #Path(f"{os.getcwd()}") / "multixyz.xyz"
    
    def __init__(self, substituentfile = "substituents.csv",
                 sub_list = default_substituents,
                 base_struct = base_structures
                 ):
        
        
        try:
            self.base_struct = base_struct
        except:
            self.base_struct = base_structures
        self.substituentfile = Path(f"{os.getcwd()}")/substituentfile
        try:
            self.substituentfile.open()
            self.sub_list = pd.read_csv(self.substituentfile).to_numpy().squeeze()
        except:
            print(self.substituentfile)
            print("no substituent list file found, please place in working directory... \n reverting to default list of substituents")
            self.sub_list = sub_list
        
        self.library = []
        self.mol2lib = []
        print(base_struct)
        print(sub_list)
        x=input("sdagasdg")
    

    def make_substituents(self):
        for base in self.base_struct:
            # Iterate through each substituent and generate the combined structure
            for substituent in self.sub_list:
                # Replace "X" in the base structure with the substituent
                combined_structure = base.replace("X", substituent).replace("Y", "")
                print(combined_structure)
                
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
                    self.library.append((mol.GetNumAtoms(), combined_structure, xyz_coordinates))      
                    self.mol2lib.append(mol2_block)
                    
                else:
                    print(f"Invalid SMILES: {combined_structure}")        
    
    def write_xyzs(self):
        Path('libgen').mkdir(exist_ok=True)
        #Path('fragment_mol2').mkdir(exist_ok=True)
        Path('fragment_mol2/xyzs').mkdir(parents=True, exist_ok=True)
        
        multixyz = Path(f"{os.getcwd()}") / "multixyz.xyz"
        with multixyz.open(mode="w+") as f:
            for atoms, smiles, xyz in self.library:
                f.write(f"{atoms}\n{smiles}\n{self.pretty_print_xyz(xyz)}\n")
                print("writing multixyz")
        
        # Lets print all XYZ coordinates to many files and then transform them with obabel into mol2 files
        
        
        for atoms, smiles, xyz in self.library:
            _xyztmp = Path(f"{os.getcwd()}/fragment_mol2/xyzs") / f"{smiles.replace('(', '_').replace(')', '_').replace('[', '_').replace(']', '_')}_tmp.xyz"
            with _xyztmp.open(mode="w+") as f:
                print(atoms)
                f.write(f"{atoms}\n{smiles}\n{self.pretty_print_xyz(xyz)}\n")
                print("writing")
        
        #convert to sybyl mol2 files and erase previous .xyz files
        #for i, f in enumerate(glob.glob("/fragment_mol2/xyzs/*_tmp.xyz")):
        
        
        for i, f in enumerate(Path(f"{os.getcwd()}/fragment_mol2/xyzs").glob('*_tmp.xyz')):
            print("found", i, f)
            #outputfile = Path(f"{os.getcwd()}/fragment_mol2") / f"obabel_\"{f.name}\"_to.mol2"
            outputfile = Path(f"{os.getcwd()}/fragment_mol2") / f"obabel_{f.name}_to.mol2"
            obabel_command = f"obabel \"{f}\" -O {outputfile}"
            print(obabel_command)
            subprocess.run(obabel_command, shell=True, check=True)
            #subprocess.run(f"mv \"{f}\"", shell=True, check=True)


    def pretty_print_xyz(self, xyz_coordinates):
        lines = [f'{symbol} {x:.4f} {y:.4f} {z:.4f}' for symbol, x, y, z in xyz_coordinates]
        return '\n'.join(lines)

    


