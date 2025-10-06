import rdkit
from rdkit import Chem
import networkx as nx

typical_valences = {
    'C': 4,
    'N': 3,
    'O': 2,
    'S': 6,
    'P': 5,
    'H': 1,
    'F': 1,
    'Cl': 1,
    'Br': 1,
    'I': 1,
    'B': 3,
    'Pd': 4
}

#compare the molecule against valences:
def swap_high_valent(mol, swapped_atom = 16, atoms = None):
    assert type(mol)==rdkit.Chem.rdchem.Mol, print("Need to pass an RDKit mol object!")
    tsmol = Chem.Mol(mol)
    TS_atoms = []
    swapped = {} # atom number: what it used to be

    if(atoms is not None):
        assert type(atoms[0])==int, print("Need to pass a list of atom numbers, as integers")
        for atom in atoms:
            tsmol.GetAtoms()[atom].SetAtomicNum(swapped_atom)
            swapped[atom] = tsmol.GetAtomWithIdx(atom).GetSymbol() 
    for atom in tsmol.GetAtoms():
        id = atom.GetIdx()
        symbol = atom.GetSymbol()
        total_valence = atom.GetTotalValence()
        try:
            if(total_valence>typical_valences[symbol]):
                TS_atoms.append(id)
                swapped[id] = symbol #keep track of the swaps
                tsmol.GetAtoms()[id].SetAtomicNum(swapped_atom) #replace with sulfur
        except:
            print(f" {symbol} : atom not in list, update it please")

    return swapped, tsmol


def same_hybrid_check(node1, node2):
    #print("here", node["hybridization"])
    return node1["hybridization"] == node2["hybridization"]

def mol_to_nx(mol):

    G = nx.Graph()  
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic()
                   #also collect the number of bonds
                   )

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    
    return G

def subgraph_matches_nx(full_nx, frag_nx):
    gm = nx.algorithms.isomorphism.GraphMatcher(frag_nx, full_nx, node_match=same_hybrid_check)
    #for subgraph in gm.subgraph_isomorphisms_iter():
    #    print("jfc", subgraph)
    #check matches with: print(sorted(list(dict.fromkeys(item.keys()))))
    return  [i for i in gm.subgraph_isomorphisms_iter()]

def rdkitmol_from_smiles(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string, sanitize=True)
    return mol

#for testing:
DSI_backbone_smi = "O=S1(C2=C(C3=C([H])C([H])=C([H])C([H])=C3C([H])=C2)C4=C5C([H])=C([H])C([H])=C([H])C5=C([H])C=C4S(N1[H])(=O)=O)=O" # "O=S1(C2=C(C3=C(C=CC=C4)C4=CC=C3S([N]1)(=O)=O)C5=CC=CC=C5C=C2)=O"  
DSI_mol = Chem.MolFromSmiles(DSI_backbone_smi, sanitize=False)