from ._core import AsyncExternalDriver
from ..dtypes import Atom, Bond, Molecule, CartesianGeometry
from ..parsing import extract_xtb_atomic_properties
from copy import deepcopy
from datetime import datetime
from glob import glob
from warnings import warn
from typing import List, Callable
from math import ceil, pi
import numpy as np
from itertools import combinations
import asyncio as aio
import re
import asyncio

class AsyncXTBDriver(AsyncExternalDriver):
    def __init__(self, name="", scratch_dir="", nprocs=1, encoding="utf8"):
        super().__init__(
            name=name, scratch_dir=scratch_dir, nprocs=nprocs, encoding=encoding
        )

    async def optimize(
        self,
        mol: Molecule,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        in_place: bool = False,
        get_energies: bool = False,
    ):
        """
        Attempt a geometry optimization with parameters specified
        """

        g0_xyz = mol.to_xyz()
        nn = mol.name
        #print("something's going wrong here")
        #x=input("0")
        
        if(get_energies):
            #print("checkpoint, get energies")
            ##x=input()
            optimized, mol.xtb_energy = await self.xyz_optimize(
                g0_xyz,
                method=method,
                crit=crit,
                xtbinp=xtbinp,
                maxiter=50,
                xyz_name=nn + "_g0.xyz",
                xtb_energies = True,
            )
            #print("OPTIMIZED", optimized)
            ##x=input()
        else:
            optimized = await self.xyz_optimize(
                    g0_xyz,
                    method=method,
                    crit=crit,
                    xtbinp=xtbinp,
                    maxiter=50,
                    xyz_name=nn + "_g0.xyz",
                )

        #print("something's going wrong here")
        ##x=input("1")
        if not in_place:
            #print("something's going wrong here")
            ##x=input("2")
            mol1 = deepcopy(mol)
            mol1.update_geom_from_xyz(optimized, assert_single=True)
            return mol1
        else:
            #print("something's going wrong here")
            ##x=input("3")
            mol.update_geom_from_xyz(optimized, assert_single=True)

    async def xyz_optimize(
        self,
        xyz: str,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        xyz_name: str = "mol",
        xtb_energies: bool = False,
    ):
        # command that will be used to execute xtb package 
        _cmd = f"""xtb {xyz_name}.xyz --{method} --opt {crit} --cycles {maxiter} {"--input param.inp" if xtbinp else ""} -P {self.nprocs}""" 
        print(_cmd)
        # PRINT _cmd
        #print(_cmd)
        ##x=input("something here? 1")
        # pylint: disable=unused-variable
        code, files, stdout, stderr = await self.aexec(
            _cmd,
            inp_files={f"{xyz_name}.xyz": xyz, "param.inp": xtbinp},
            out_files=["xtbopt.xyz"],
        )
        #print(files)
        ##x=input("something here? 2" )
        if xtb_energies:

            for l in stdout.split("\n")[::-1]:
                #print(l)
                if m := re.match(r"\s+\|\s+TOTAL ENERGY\s+(?P<eh>[0-9.-]+)\s+Eh\s+\|.*", l):
                    xtb_energy= float(m["eh"])
                    nxyz = files["xtbopt.xyz"]
                    return nxyz, xtb_energy

                
        if "xtbopt.xyz" in files and xtb_energies == False:
            nxyz = files["xtbopt.xyz"]
            return nxyz
        else:
            raise FileNotFoundError("Could not locate xtb output file.")

    async def optimize_conformers(
        self,
        mol: Molecule,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        in_place: bool = False,
    ):
        """
        Perform the same sort of optimization as used in optimize(),
        but for each conformer instead of for the molecule's main geometry.
        """
        xyzs = mol.confs_to_xyzs()
        nn = mol.name

        optimized_confs = []
        for i, xyz in enumerate(xyzs):
            mol_name = nn + f"_{i}"
            optimized = await self.xyz_optimize(
                xyz,
                method=method,
                crit=crit,
                xtbinp=xtbinp,
                maxiter=maxiter,
                xyz_name=mol_name,
            )
            optimized_confs.append(optimized)

        geoms = [CartesianGeometry.from_xyz(conf)[0][0] for conf in optimized_confs]
        if in_place:
            mol.embed_conformers(*geoms, mode="w")
            # return a value other than None so that it will display to the user as successful
            return True
        else:
            mol1 = deepcopy(mol)
            mol1.embed_conformers(*geoms, mode="w")
            return mol1
        
        """
        self,
        mol: Molecule,
        method: str = "gff",
        crit: str = "normal",
        xtbinp: str = "",
        maxiter: int = 50,
        in_place: bool = False,
        """
    async def optimize_conformers_TS(

        self,
        mol: Molecule,
        method: str="GFN2", 
        force_const: float = 0.5,
        in_place: bool = True, #False, #default is False
        maxcycle: int = 20,
        crit: str = "normal",
        TS_bonds: list = [], # this should contain a list of the bond objects that should be LOCKED, all other long bonds should be FIXED
        TS_fixed: list = [], # this should contain a list of atoms that are to be IMMOBILIZED
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        custom_constraint: list = [],
        fn_suffix: str = 0,
        
    ):
        """
        Perform the same sort of optimization as used in optimize(),
        but for each conformer instead of for the molecule's main geometry.
        """

        xyzs = mol.confs_to_xyzs()
        _xyzs = []
        nn = mol.name

        inp = f"$constrain\n force constant={force_const}\n"      
        inp += TS_bonds[0]
        inp += f"$opt\n  maxcycle={maxcycle}\n"
        inp += "$end\n"       
        print(inp)
        optimized_confs = []
        #print("XYZs", xyzs.strip("\n"))
        #xyzs.strip("\n")
        for i, xyz in enumerate(xyzs):
            mol_name = nn + f"_{i}"
            optimized = await self.xyz_optimize(
                xyz,
                method=method,
                crit=crit,
                xtbinp=inp,
                #maxiter=maxiter,
                xyz_name=mol_name,
            )
            #optimized_confs contains the right number of xyz coordinates after the for loop is complete
            optimized_confs.append(optimized)
        #tgeoms contains the right number of CartesianGeometry objects
        geoms = [CartesianGeometry.from_xyz(conf)[0][0] for conf in optimized_confs]

        
        if in_place:
            mol.embed_conformers(*geoms, mode="w")
            # return a value other than None so that it will display to the user as successful
            return True
        else: #this is default, as inplace=False
            mol1 = deepcopy(mol)
            mol1.embed_conformers(*geoms, mode="w")
            return mol1

    """
    async def constrain_opt_TS(
        self,
        mol: Molecule,
        method: str="GFN2", 
        force_const: float = 0.5,
        in_place: bool = False,
        maxcycle: int = 20,
        TS_bonds: list = [], # this should contain a list of the bond objects that should be LOCKED, all other long bonds should be FIXED
        TS_fixed: list = [], # this should contain a list of atoms that are to be IMMOBILIZED
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        custom_constraint: list = [],
        fn_suffix: str = 0,
    ):

        print("CHECKING CHECKING CHECKING")
        ##x=input(f"1: {mol.name}")
        inp = f"$constrain\n force constant={force_const}\n"      
        print(inp)
        inp += TS_bonds[0]
        print(inp)
        inp += f"$opt\n  maxcycle={maxcycle}\n"
        print(inp)
        inp += "$end\n"       
        print(inp)
        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
            get_energies=True,
        )     

        return m1
    """
    # KSP 2023-10-19
    ## write a function that fixes long bonds for the TS structures. 
    """
    async def fix_long_bonds_TS(
        self,
        mol: Molecule,
        method: str = "gff",
        rss_length_thresh: float = 4.0,
        rss_steps: float = 15,
        rss_maxcycle: int = 20,
        force_const: float = 0.5,
        target_len: float = 1.5,
        in_place: bool = False,
        TS_bonds: list = [], # this should contain a list of the bond objects along the breaking/forming profile (user located)
        TS_bond_lengths: list = [], #this is an ordered list with desired bond lengths (they should be placed in order of bond number!) ## for now...
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        fn_suffix: str = 0,
        ):
        ### test inputs
        TS_bond_lengths = [1.3, 1.4]
        
        ### /test inputs
    """
    #  heres the plan:
    # Need to take in the bond types that are involved in the transition state (user identified)
    # Need to take in all other long bonds 
    
    # Need to generate the TS properly
    """
        
        #print("xtb.py, CHECKING STAGE 1")
        #this finds all long bonds that are related to the TS
        TS_tbf = {}
        for b in mol.bonds:
            if b not in TS_tbf and b in TS_bonds:
                TS_tbf[b] = mol.get_bond_length(b)
        
        
        #print("xtb.py, CHECKING STAGE 2")
        #this collects all the remaining long bonds
        tbf = {}  # to be fixed
        for b in mol.bonds:
            if b not in tbf and mol.get_bond_length(b) >= rss_length_thresh and b not in TS_tbf:
                tbf[b] = mol.get_bond_length(b) # create a dictionary with all bond objects that are "long"
       
        
        #print("xtb.py, CHECKING STAGE 3")
        
        #create input file, constraints need to be set differently for the TS
        #use mode=sequential, optimize the alkyne addition first, then the N--H--C bonds
        
        inp = f"$constrain\n force constant={force_const}\n"
        
        #TS_lb_atoms = set()
        #lb_atoms = set()
      
        for b in TS_tbf:
            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            #print("xtb.py, CHECKING STAGE 4")  
            #print(a1, a2) 
            inp += f"  distance: {a1+1}, {a2+1}, {TS_tbf[b]:0.4f}\n" ## this line is the fucking problem. WTF.

            #TS_lb_atoms.add(b.a1)
            #print("xtb.py, CHECKING STAGE 5")  
        
           
        for c in tbf:
            #print("xtb.py, CHECKING STAGE 6")  
            a1, a2 = mol.get_atom_idx(c.a1), mol.get_atom_idx(c.a2)
            #print("xtb.py, CHECKING STAGE 7") 
            #print(a1, a2)
            inp += f"  distance: {a1+1}, {a2+1}, {tbf[c]:0.4f}\n"
            #print("xtb.py, CHECKING STAGE 8")
            #lb_atoms.add(c.a1)
            #lb_atoms.add(c.a2)
            #print("xtb.py, CHECKING STAGE 9") 
            
        inp += "$scan\n  mode=concerted\n"

        #for i, b in enumerate(TS_tbf):
        #    inp += f"  {i+1}: {TS_tbf[b]:0.4f}, {TS_bond_lengths[i]:0.4f}, {rss_steps}\n"
        for i, b in enumerate(tbf):
            inp += f"  {i+1}: {tbf[b]:0.4f}, {target_len:0.4f}, {rss_steps}\n"        
        inp += "$end\n"

        inp += f"$opt\n  maxcycle={rss_maxcycle}\n"
        inp += "$end\n"       

        
        #print("CHECKING CHECKING INPUT", inp)
        ##x=input()
        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
        )
        #print("xtb.py, CHECKING STAGE 5")
        #print("CHECKING CHECKING CHECKING", m1)
        return m1
    """

    async def constrain_opt_TS(
        self,
        mol: Molecule,
        method: str="GFN2", 
        force_const: float = 0.5,
        in_place: bool = False,
        maxcycle: int = 20,
        TS_bonds: list = [], # this should contain a list of the bond objects that should be LOCKED, all other long bonds should be FIXED
        TS_fixed: list = [], # this should contain a list of atoms that are to be IMMOBILIZED
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        custom_constraint: list = [],
        fn_suffix: str = 0,
    ):

        print("CHECKING CHECKING CHECKING")
        ##x=input(f"1: {mol.name}")
        inp = f"$constrain\n force constant={force_const}\n"      
        print(inp)
        inp += TS_bonds[0]
        print(inp)
        inp += f"$opt\n  maxcycle={maxcycle}\n"
        print(inp)
        inp += "$end\n"       
        print(inp)
        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
            get_energies=True,
        )     

        return m1
        
          
################################################################################################################
    async def fix_long_bonds_TS2(
        self,
        mol: Molecule,
        method: str = "gff",
        rss_length_thresh: float = 4.0,
        rss_steps: float = 20,
        rss_maxcycle: int = 20,
        force_const: float = 0.5,
        target_len: float = 1.5,
        in_place: bool = False,
        TS_bonds: list = [], # this should contain a list of the bond objects that should be LOCKED, all other long bonds should be FIXED
        TS_fixed: list = [], # this should contain a list of atoms that are to be IMMOBILIZED
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        custom_constraint: list = [],
        fn_suffix: str = 0,
    ):
        ### test inputs

        
        ### /test inputs
        """ heres the plan:
        Need to take in the bond types that are involved in the transition state (user identified)
        Need to take in all other long bonds 
        
        Need to generate the TS properly
        """
        
        ##x=input("xtb.py, CHECKING STAGE 1")
        
        #this finds all long bonds that are related to the TS
        TS_tbf = TS_bonds
        """
        TS_tbf = {}
        for b in mol.bonds:
            if b in TS_bonds:
                TS_tbf[b] = mol.get_bond_length(b)
                #print(TS_tbf)
        """
        
        ##x=input("xtb.py, CHECKING STAGE 2")
        #this collects all the remaining long bonds
        tbf = {}  # to be fixed
        for b in mol.bonds:
            if b not in tbf and mol.get_bond_length(b) >= rss_length_thresh and b not in TS_tbf:
                tbf[b] = mol.get_bond_length(b) # create a dictionary with all bond objects that are "long"
                print("mol.bonds check", tbf[b])
                
        print("CHECKING \t\t\t\t\t tbf", tbf)
       
        
        ##x=input("xtb.py, CHECKING STAGE 3")
        
        #create input file, constraints need to be set differently for the TS
        #use mode=sequential, optimize the alkyne addition first, then the N--H--C bonds
        
        inp = f"$constrain\n force constant={force_const}\n"
        
        #TS_lb_atoms = set()
        #lb_atoms = set()
        #print("CHECKING tbf", TS_tbf)

        for c in tbf:
            ##x=input("xtb.py, CHECKING STAGE 4")  
            a1, a2 = mol.get_atom_idx(c.a1), mol.get_atom_idx(c.a2)
            ##x=input("xtb.py, CHECKING STAGE 5") 
            #print(a1, a2)
            inp += f"  distance: {a1+1}, {a2+1}, {tbf[c]:0.4f}\n"
            ##x=input("xtb.py, CHECKING STAGE 6")
            #lb_atoms.add(c.a1)
            #lb_atoms.add(c.a2)
        print("about to check TS_tbf", TS_tbf)
        ##x=input("xtb.py, CHECKING STAGE 7") 
        
        
        
        inp += TS_bonds[0]


        ##x=input("xtb.py, CHECKING STAGE 9")  
        
           

            
        inp += "$scan\n  mode=concerted\n"

        #for i, b in enumerate(TS_tbf):
        #    inp += f"  {i+1}: {TS_tbf[b]:0.4f}, {TS_bond_lengths[i]:0.4f}, {rss_steps}\n"
        ##x=input("xtb.py, CHECKING STAGE 10")  
        for i, b in enumerate(tbf):
            inp += f"  {i+1}: {tbf[b]:0.4f}, {target_len:0.4f}, {rss_steps}\n"        
        ##x=input("xtb.py, CHECKING STAGE 11")  
        inp += "$end\n"

        inp += f"$opt\n  maxcycle={rss_maxcycle}\n"
        inp += "$end\n"       

        
        print(inp)
        ##x=input(f"SENT {molecule.name} to self.optimize")        
        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
        )

        ##x=input("xtb.py, CHECKING STAGE 10")

        #print("CHECKING CHECKING CHECKING", m1)
        ##x=input()
        return m1


    async def fix_long_bonds(
        self,
        mol: Molecule,
        method: str = "gff",
        rss_length_thresh: float = 4.0,
        rss_steps: float = 15,
        rss_maxcycle: int = 20,
        force_const: float = 0.5,
        target_len: float = 1.5,
        in_place: bool = False,
        constrain_bonds: list = ["C-C", "C-H", "C-F", "C-O", "O-H", "N-Cl"],
        fn_suffix: str = 0,
    ):
        """
        Fix all long bonds in the molecule by doing a relaxed surface scan with coordinates constrained
        """
        tbf = {}  # to be fixed
        for b in mol.bonds:
            #print("bond length for", mol.get_bond_length(b))
            if b not in tbf and mol.get_bond_length(b) >= rss_length_thresh:
                tbf[b] = mol.get_bond_length(b) # create a dictionary with all bond objects that are "long"

        if not tbf:
            return mol
        inp = f"$constrain\n  force constant={force_const}\n" # force constant...

        lb_atoms = set() # creates the unique items

        for b in tbf:
            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            inp += f"  distance: {a1+1}, {a2+1}, {tbf[b]:0.4f}\n"
            lb_atoms.add(b.a1)
            lb_atoms.add(b.a2)

        # why is this part even necessary ?!?!
        """
        # generate constraints for C-H bonds
        core_bonds = tuple(mol.yield_bonds(*constrain_bonds))
        inp += self.gen_bond_constraints(mol, core_bonds)
        inp += self.gen_angle_constraints(mol, lb_atoms)
        """
        
        
        inp += "$scan\n  mode=concerted\n"
        for i, b in enumerate(tbf):
            inp += f"  {i+1}: {tbf[b]:0.4f}, {target_len:0.4f}, {rss_steps}\n"

        inp += f"$opt\n  maxcycle={rss_maxcycle}\n"
        inp += "$end\n"

        print("THIS IS THE INPUT FILE:", inp)
        ###x=input()
        m1 = await self.optimize(
            mol,
            method=method,
            crit="crude",
            xtbinp=inp,
            in_place=False,
        )

        return m1

    async def xyz_energy(self, xyz: str, method: str = "gfn2", accuracy: float = 1.0):
        _cmd = f"""xtb struct.xyz --{method} --acc {accuracy:0.2f}"""

        code, files, stdout, stderr = await self.aexec(
            _cmd, inp_files={f"struct.xyz": xyz}
        )

        # This is what we are trying to find in the output file
        # | TOTAL ENERGY             -172.541095318001 Eh   |

        for l in stdout.split("\n")[::-1]:
            if m := re.match(r"\s+\|\s+TOTAL ENERGY\s+(?P<eh>[0-9.-]+)\s+Eh\s+\|.*", l):
                return float(m["eh"])
                
    """ #this is the original function
    async def conformer_energies(
        self, mol: Molecule, method: str = "gfn2", accuracy: float = 1.0
    ):

        xyzs = mol.confs_to_xyzs() # get all conformers associated with the molecule
        nn = mol.name #name of molecule...
        energies = [] #empty array for energies

        for i, xyz in enumerate(xyzs): #iterate through each conformer xyz file
            #print("conformer", i)
            #x=input("check conf #")
            #print(len(xyzs))
            conf_energy = await self.xyz_energy(xyz, method=method, accuracy=accuracy)
            energies.append(conf_energy)

        ref_energy = energies[0]
        #print("ref_energy", ref_energy)
        #print("energies", energies)

        #return (np.array(energies) - ref_energy) * 2625.5  # conversion to kJ/mol
        return (np.array(energies) - ref_energy) * 627.51  # conversion to kcal/mol

    """
    
    #this is the new version
    async def conformer_energies_independent(
        self, mol: Molecule, method: str = "gfn2", accuracy: float = 1.0
    ):
        xyzs = mol.confs_to_xyzs() # get all conformers associated with the molecule
        nn = mol.name #name of molecule...
        energies = [] #empty array for energies
        #tasks = [( _full_mol_copy, _full_mol_frag, random.randrange(0, 10000)) for _ in range(8)]
        tasks = [self.xyz_energy(xyz, method=method, accuracy=accuracy) for xyz in xyzs]
        
        energies = await asyncio.gather(*tasks)
        
        print(energies)
        print(type(energies))
        
        if(None in energies):
            print("WARNING: One of the geometries failed to minimize!")
        energies = [x for x in energies if x != None]
        if(len(energies)<1):
            print("ERROR: NO USEFUL GEOMETRIES. Currently there is no mechanism in place to handle all xtb minimizations failing")
            print("THIS SHOULD BE IMPLEMENTED LATER...")
            return False # KSP addition 6/26/2024
        ref_energy = energies[0]
        
        
        """
        for i, xyz in enumerate(xyzs): #iterate through each conformer xyz file
            tasks = [
            #print("conformer", i)
            #x=input("check conf #")
            #print(len(xyzs))
            conf_energy = await self.xyz_energy(xyz, method=method, accuracy=accuracy)
            energies.append(conf_energy)

        ref_energy = energies[0]
        #print("ref_energy", ref_energy)
        #print("energies", energies)
        """
        #return (np.array(energies) - ref_energy) * 2625.5  # conversion to kJ/mol
        return (np.array(energies) - ref_energy) * 627.51  # conversion to kcal/mol



    async def charges(
        self,
        mol: Molecule,
        method: str = "gfn2",
        accuracy: float = 0.5,
        net_charge: int = 0,
    ):
        """
        Compute atomic charges using XTB methodology

        FIXED 2022-02-23: added at least a partial support for total charge of the molecule.
        DO NOT USE unless you know what you are doing. -SAS
        """
        xyz = mol.to_xyz(n=0)
        _cmd = f"""xtb struct.xyz --sp --{method} --acc {accuracy:0.2f} --chrg {net_charge}"""

        code, files, stdout, stderr = await self.aexec(
            _cmd, inp_files={f"struct.xyz": xyz}, out_files=["charges"]
        )

        charges = np.array(list(map(float, files["charges"].split())), dtype=np.float32)

        return charges

    #### Ian dev on roche atomic properties

    async def atom_prop(self, xyz: str, method: str = "gfn2", accuracy: float = 0.5):
        """
        Compute fukui indices, polarizability, charge, dispersion coeff, and max wiberg bond order for atoms

        """
        _cmd = (
            f"""xtb struct.xyz --{method} --acc {accuracy:0.2f} --vfukui > xtbout.log"""
        )
        code, files, stdout, stderr = await self.aexec(
            _cmd, inp_files={f"struct.xyz": xyz}, out_files=["xtbout.log"]
        )
        # output = stdout
        # print(files)
        # print(stdout)
        outdf = extract_xtb_atomic_properties(files["xtbout.log"])
        return outdf

    async def conformer_atom_props(
        self, mol: Molecule, method: str = "gfn2", accuracy: float = 1.0
    ):
        """
        Compute fukui indices, polarizability, charge, dispersion coeff, and max wiberg bond order for atoms of
        all conformers.
        """
        xyzs = mol.confs_to_xyzs()
        conf_props = []
        for i, xyz in enumerate(xyzs):
            # print(xyz)
            prop_df = await self.atom_prop(xyz, method=method, accuracy=accuracy)
            conf_props.append(prop_df)
        return conf_props

    #### Ian dev end


    # KSP 2023-10-19
    """
    def gen_ts_constraints(mol: Molecule, bonds: List[Bond]):
        constr = """
    
    # /KSP
    
    @staticmethod
    def gen_bond_constraints(mol: Molecule, bonds: List[Bond]):
        """Generate bond distance constraint list"""
        constr = ""
        for b in bonds:
            a1, a2 = mol.get_atom_idx(b.a1), mol.get_atom_idx(b.a2)
            constr += f"  distance: {a1+1}, {a2+1}, {mol.get_bond_length(b):0.4f}\n"
        
        print("CHECKING CONSTRAINTS", constr)
        return constr

    @staticmethod
    def gen_angle_constraints(mol: Molecule, atoms: List[Atom]):
        """Generate constraints for all angles where atom is the middle atom"""
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
