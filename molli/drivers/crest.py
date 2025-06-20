from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from math import ceil, pi
from itertools import combinations
import asyncio as aio

def _parse_energies(_n: str):
    """
        Parse energies in the form of
        1   0.000
        2   1.131
    """
    energies = []
    for l in _n.splitlines(keepends=False):
        if l:
            e = float(l.split()[1])
            energies.append(e)
    
    return energies



class AsyncCRESTDriver(AsyncExternalDriver):
    # def __init__(self, name="", scratch_dir="", nprocs=1, encoding="utf8"):
    #     super().__init__(
    #         name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
    #     )
############################################################################################################################################################################    
###################################################### conformer_search_TS #################################################################################################   
 
    async def conformer_search_TS(
        self,
        mol: Molecule,
        method: str = "gfnff",
        ewin: float = 6.0,
        mdlen: float = 15.0,
        mddump: float = 250.0,
        vbdump: float = 1.0,
        chk_topo: bool = False,
        global_lock: list = [], ## atoms to lock (in string)\
        metadyn: list = [], ## atoms to allow metadyn on
        reference_coord: list = [], ## reference XYZ coordinates
        constr_val_angles: list = ['P', 'Cu'],
        force_const: float = 0.05
    ):
        """
        `ewin`: energy window in kcal/mol
        `mdlen`: iMTD-GC molecular dynamics length
        """
        print("checkpoint 1")

        
        g0_xyz = mol.to_xyz()

        nn = mol.name
        print("checkpoint 1.5")
        ### GENERATE VALENCE ANGLE CONSTRAINTS
        atoms_constr = []
        
        cinp = "$constrain\n"
        cinp += f"{global_lock[0][:-6]}\n" #
        cinp += "force constant=6.0\n"
        cinp += f"reference=coord.ref\n"
        cinp += "$metadyn\n"  
        cinp += metadyn
        print("checkpoint 2")
        rinp = "$coord\n"
        rinp += reference_coord
        rinp += "$end\n"
        #print("CHECKING CREST")
        #print(cinp)
        #print(rinp)

        print("CHECKING CINP\n", cinp, "\n CHECKING RINP\n", rinp)
        # command that will be used to execute xtb package
        #_cmd = f"""crest coord {nn}_g0.xyz -{method} -quick -ewin {ewin:0.4f} -mdlen {mdlen:0.4f} -mddump {mddump:0.4f} -vbdump {vbdump:0.4f} -T {self.nprocs}"""
        #_cmd = f"""crest {nn}_g0.xyz --constrain {global_lock[0][:-1]}"""
        _cmd = f"""crest coord --cinp constraints.inp """
        print(_cmd)
        x=input("check _cmd")
        # append


        
        
        #_cmd += f" -cinp {nn}_constraint.inp -fc {force_const:0.4f}"
        
        ### EXECUTE CREST CODE
        
        print("CHECKPOINT 1")
        

        print("CHECKING INP\n", rinp, "\n CHECKING _CONSTRAINTS:", cinp)
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"coord": rinp, f"coord.ref": rinp, f"constraints.inp": cinp},#, f"coord": rinp},
            out_files=["crest_conformers.xyz"],
        )
        print("CHECKPOINT 2", mol.name)
        print("Code", code, "files", files, "stdout", stdout, "stderr", stderr)
        #for item in zip(mol.atoms,mol.geom.coord.tolist()):
        #    print(item)
        try:
            print("LOOKING FOR CONFORMER FILE")
            print("FILES:", files)
            ens1 = files["crest_conformers.xyz"]
        except:
            print("EXCEPTION RAISED")
            raise FileNotFoundError("crest_conformers.xyz")

        geoms = [x for x, _, _ in CartesianGeometry.from_xyz(ens1)]
        mol.embed_conformers(*geoms, mode="w")
        print("CHECKPOINT 3")
        return mol
        
############################################################################################################################################################################    
############################################################################################################################################################################    
    
    async def conformer_search(
        self,
        mol: Molecule,
        method: str = "gfnff",
        ewin: float = 6.0,
        mdlen: float = 15.0,
        mddump: float = 250.0,
        vbdump: float = 1.0,
        chk_topo: bool = False,
        constr_val_angles: list = ['P', 'Cu'],
        force_const: float = 0.05
    ):
        """
        `ewin`: energy window in kcal/mol
        `mdlen`: iMTD-GC molecular dynamics length
        """
        g0_xyz = mol.to_xyz()

        nn = mol.name

        ### GENERATE VALENCE ANGLE CONSTRAINTS
        atoms_constr = []
        for a in constr_val_angles:
            atoms_constr.extend(mol.get_atoms_by_symbol(a))
        cinp = "$constrain\n"
        cinp += self.gen_angle_constraints(mol, atoms_constr)
        cinp += "$end\n"

        # command that will be used to execute xtb package
        _cmd = f"""crest {nn}_g0.xyz -{method} -quick -ewin {ewin:0.4f} -mdlen {mdlen:0.4f} -mddump {mddump:0.4f} -vbdump {vbdump:0.4f} -T {self.nprocs}"""

        # append

        if not chk_topo:
            _cmd += " --noreftopo"
        
        if constr_val_angles:
            _cmd += f" -cinp {nn}_constraint.inp -fc {force_const:0.4f}"
        
        ### EXECUTE CREST CODE
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}_g0.xyz": g0_xyz, f"{nn}_constraint.inp": cinp},
            out_files=["crest_conformers.xyz"],
        )

        try:
            ens1 = files["crest_conformers.xyz"]
        except:
            raise FileNotFoundError("crest_conformers.xyz")

        geoms = [x for x, _, _ in CartesianGeometry.from_xyz(ens1)]
        mol.embed_conformers(*geoms, mode="w")

        return mol

    async def confomer_screen(
        self,
        mol: Molecule,
        method: str = "gfn2",
        ewin: float = 6.0,
    ):
        """
        Any conformer ensemble present in a molecule is reoptimized with the selected method, 
        then pruned and sorted by energies. Useful in eliminating redundancies from deficiencies of GFN-FF, for instance.
        """
        nn = mol.name
        confs = mol.confs_to_multixyz()

        _cmd = f"""crest -screen {nn}_confs.xyz -{method} -ewin {ewin} -T {self.nprocs} """

        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{nn}_confs.xyz": confs},
            out_files=["crest_ensemble.xyz", "crest.energies"]
        )
        
        ens1 = files["crest_ensemble.xyz"]
        nrgs1 = files["crest.energies"]

        geoms = [x for x, _, _ in CartesianGeometry.from_xyz(ens1)]

        _mol = deepcopy(mol)

        _mol.embed_conformers(*geoms, mode="w")
        return _mol

        ### ADD CONFORMER PROPERTIES !!! ###
        # return (code, files, stdout, stderr) 



    @staticmethod
    def gen_angle_constraints(mol: Molecule, atoms: List[Atom]):
        """ Generate constraints for all angles where atom is the middle atom """
        constr = ""
        for a in atoms:
            neigbors = mol.get_connected_atoms(a)
            for a1, a2 in combinations(neigbors, 2):
                i1 = mol.get_atom_idx(a1) + 1
                i2 = mol.get_atom_idx(a) + 1
                i3 = mol.get_atom_idx(a2) + 1
                angle = mol.get_angle(a1, a, a2) * 180 / pi
                constr += f"  angle: {i1}, {i2}, {i3}, {angle:0.4f}\n"

        return constr