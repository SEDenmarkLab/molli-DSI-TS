#!/usr/bin/env python

import numpy as np
import glob
from copy import deepcopy
import pprint as pp
import pandas as pd

#$% extract important information from NBO files (from orca output)


class NBOExtract():
    def __init__(self, software: str="orca", endswith: str="out"):
        self.software = software
        self.endswith = endswith
        self.flags()
        self.files = self.get_files()
        self.features = ["npa", "nbbp", "e2pert"] ## these are the functions to be called...
        self.feature_methods = []
        self.feature_parsers = []
        for feature in self.features:
            self.feature_methods.append(getattr(self, feature))
        for feature in self.features:
            self.feature_parsers.append(getattr(self,f"{feature}_parse"))
        #the last entry are the elements
        self.all_data_container ={key:[x for x in range(len(self.features)+1)] for key in self.files}# length depends on number of files + keep one for atoms
        self.parsed_all_data_container = deepcopy(self.all_data_container)
        """
        self.parsed_all_data_container organization:
        [[npa data],[nbbp,nbbp_bonds],[e2pert],[atomlist]]
        
        """
        
        
    def flags(self):

        #flag patterns: ["start", "end", "datastart", "dataend"]
        if(self.software == "orca"):
            self.npa_flag = ["  Atom No    Charge        Core      Valence    Rydberg      Total",
                             "Natural Rydberg Basis",
                             "--------------------------------------------------------------------",
                             "===================================================================="
                             ]
            self.nbbp_flag = ["NBBP: NBO bond-bond polarizability matrix",
                              "\n",
                              "------------ ------- ------- ------- ------- ------- ------- ------- -------", # accomodate 8 minimum atoms...
                              "\n", #no clear terminator, only blank lines
                              ]
            self.e2pert_flag = ["SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS",
                                "\n",
                                "===============================================================================",
                                "\n"                              
                                ]
        self.all_flags = [self.npa_flag, self.nbbp_flag, self.e2pert_flag]

    def get_files(self):
        return glob.glob(f"*{self.endswith}")

    def read_all(self):
        for i,file in enumerate(self.files):
            print("... reading file", file)
            self.all_data_container[file] = self.readtext(file)

    def parse_all(self):
        for i, item in enumerate(self.all_data_container): #will be the associated filename
            for j, func in enumerate(self.feature_parsers): #the associated method..
                self.parsed_all_data_container[item][j] = func(self.all_data_container[item][j]) #send the datacontainer for the file
                #print(self.parsed_all_data_container[item][j])

           

            atoms = np.array(self.parsed_all_data_container[item][self.features.index("npa")][0:])

            _atom_collec = []
            for atom in atoms:
                _atom_collec.append(atom[0][0])

            self.parsed_all_data_container[item][len(self.feature_parsers)] = _atom_collec
            #self.parsed_all_data_container[item][len(self.feature_parsers)] = np.array(atoms[:,:0].squeeze()).tolist()
            

    def readtext(self, file):
        flag = -1 # nothing is hit
        lastflag = -1
        collecting = False
        datacontainer = [[x] for x in range(len(self.feature_methods))] #for each file
        with open(file) as f:
            filetext = f.readlines()
        for i,line in enumerate(filetext):
            if(i<100): #skip the first 100 lines, 100 is arbitrary, but it will never contain useful information for the NBO output
                continue
            flag = self.check([line, filetext[i-1]]) #check the current line, and the line before it
            #four states:
            # 1. hit a specific flag (with a unique string)
            # 2. inbetween flags (after hitting positive flag, an initiator)
            # 3. inbetween flags (after hitting a negative flag, a terminator)
            # 4. hit a flag with a generic string (e.g., '\n'), this requires special handling
            if(collecting==True and flag==-1): #self.check(line)==-1): #inbetween flags (after hitting a positive flag)
                datacontainer[lastflag].append(line) 
            elif(collecting==True and flag>=0 and flag<100): #self.check(line)>=0): #hit the terminator flag
                datacontainer[flag].append(line)
                collecting=False
                lastflag = -1
            elif(collecting==False and flag>=0 and flag<100): #self.check(line)>00): #hit the initiator flag
                datacontainer[flag].append(line)
                lastflag = deepcopy(flag)
                collecting = True
            elif(collecting==True and flag>=100): #hit a terminator that is simply a "\n"
                datacontainer[int(flag/100)].append(line)#append the final line and stop collecting
                collecting=False
        return datacontainer
    
    def check(self, line = None):
        assert line[0] != None
        to_return = None
        blank_return = [False, -1] #-1 is arbitrary...
        for f, func in enumerate(self.feature_methods):
            if(func(line) == True): # call each of the functions listed above...
                return f #return the feature_method-associated number .. do not need to scan other functions here...
            if(func(line) == False): # we need to allow it to search through the rest of the functions..
               to_return = -1 #return a "no"
            if(func(line) == -1):
                blank_return = [True, f*100]
        if(blank_return[0] == True):
            return f*100
        return to_return #everything returned -1


    def npa(self,line):
        if(self.npa_flag[0] in line[0] or self.npa_flag[1] in line[0]):
            return True #self.features.index("npa")
        else:
            return False
    def npa_parse(self, datacontainer):
        double_letters = ["Cl", "Br", "Si"]
        npa_raw = (datacontainer[3:])[:-10]
        npa_data = []
        for i, line in enumerate(npa_raw):
            _line = []
            _tmp = line.strip('\n').split()[:3]
            if(len(_tmp[0])>2):
                _line = [_tmp[0], _tmp[1]]
            if(len(_tmp[0])==1):
                _line = [''.join(_tmp[0:2]), _tmp[2]]
            for dl in (double_letters):
                if(dl in _tmp[0]):

                    _line = [_tmp[0], _tmp[1]]

            npa_data.append(_line)
            _line = None

        return npa_data


    #contains hacky legacy code, but works (mostly)
    def e2pert_parse(self, datacontainer):
        

        prnt =1
        energycap = 10.0
        atom_numbers = [1,2,3,4,5,6]
        NBO6=1

        perturbnrg_list = datacontainer


        atoms = ['C', 'H', 'N', 'O', 'F', 'S', 'I']
        alkynewarning = False
        display_by_atoms = []
        print_list = ""
        atom_selected_print_list = ""
        pi_NBO = []
        pruned = []

        e2pert_donor_types = []
        e2pert_donor_atoms = []
        e2pert_acceptor_types = []
        e2pert_acceptor_atoms = []
        e2pert_energies = []
        e2pert_brief = []
        if(int(prnt) == 2 or int(prnt)==1):
            if(NBO6 == 1):

                data_list = []
                for i in range(0, len(perturbnrg_list) - 1):
                    curr_orbs = [-1,-1]
                    curr_nums = []
                    try:
                        pruned = perturbnrg_list[i].split()[:-3]
                    except:
                        pass
                    if(pruned):
                        pass
                    else:
                        continue
                    #first orbital is the easiest
                    curr_orbs[1]=pruned[0]
                    pruned_no_brackets = []

                    for item in pruned:

                        if("BD*" in item):
                            pruned_no_brackets.append(item[:-1])
                        if "(" in item or ")" in item: # or "BD" in item or "NL" in item:
                            pass
                        else:
                            pruned_no_brackets.append(item)
                    
                    _tmpatom = []
                    _types = []
                    for j, item in enumerate(pruned_no_brackets):
                        
                        if ("BD" in item or "BD*" in item) :

                            _types.append(2)
                        if ("LP" in item or "RY" in item):
                            _types.append(1)
                        for atom in atoms:
                            if atom in item:
                                if(len(item)==1 and not pruned_no_brackets[j+1][0].isalpha()):
                                    _tmpatom.append(f"{item}{pruned_no_brackets[j+1]}")
                                elif(len(item)>1):

                                    _tmpatom.append(f"{item}")

                    groups = [[-1,-1],[-1,-1]]

                    try:
                        if(_types[0] == 1):
                            _tmpatom.insert(1," ")
                    except:
                        pass
                    try:
                        if(_types[1] == 1):
                            
                            _tmpatom.append(" ")
                    except:
                        pass
                    for j,count in enumerate(_types):
                        
                        for k in range(count):
                            
                            try:
                                groups[j][k] = _tmpatom[int(f"{j}{k}", 2)] #binary magic!
                            except:
                                pass

                    _grp_nums = []
                    atom_groups = []
                    for group in groups:
                        for item in group:
                            if(item != -1):
                                item.strip("-")
                                _grp = ''.join([char for char in item if char.isdigit()])
                            _grp_nums.append(_grp)
                        atom_groups.append(_grp_nums)
                        _grp_nums =[]
                        _grp = []


                    printvalue = False
                    pi_NBO_temp = []
                    perturb_temp = perturbnrg_list[i].split() 
                    
                    if("within" in perturb_temp[0]):
                        continue
                    if(perturb_temp[0][0].isalpha()):
                        continue

                    eCap = True
                    if(perturb_temp[1] == 'BD'):
                        if(perturb_temp[3] == '2)' or perturb_temp[3] ==  '3)'):
                            if(float(perturb_temp[-3]) > float(energycap)):
                                            

                                            eCap = False
                                            print_list += perturbnrg_list[i]
                                            data_list.append(perturbnrg_list[i])
                                            pi_NBO_temp.append(perturb_temp[0])
                                            pi_NBO_temp.append(perturb_temp[8])
                                            pi_NBO.append(pi_NBO_temp)
                    if(perturb_temp[1] == 'LP'):
                        if(perturb_temp[4] == ' O' and perturb_temp[3] == '2)'):
                            printvalue = True
                        if(perturb_temp[4] == 'N'):
                            printvalue = True
                        if((float(perturb_temp[-3]) > float(energycap)) and (printvalue == True)):

                            eCap = False
                            print_list += perturbnrg_list[i] 
                            data_list.append(perturbnrg_list[i])
                            pi_NBO_temp.append(perturb_temp[0])
                            pi_NBO_temp.append(perturb_temp[6])
                            pi_NBO.append(pi_NBO_temp)
                    if(perturb_temp[3] == '3)'):
                        alkynewarning = True


                    try:
                        _brief = [atom_groups[0][0],atom_groups[0][1],atom_groups[1][0], atom_groups[1][1], float(perturb_temp[-3])]
                    except:
                        continue
                    e2pert_brief.append(_brief)
                    _brief = []
                return e2pert_brief #e2pert_data

    #contains hacky legacy code, but mostly works
    def get_e2pert_data(self, fname, atom_numbers, energycap, header = True, initiators_only = False):
        
        #prnt, energycap,  *atom_numbers
        prnt =1
        energycap = 10.0
        atom_numbers = atom_numbers
        maxval = energycap
        NBO6=1

        atom_selected_print_list = ""
        if(header):
            atom_selected_print_list += "                                                           E(2)  E(j)-E(i) F(i,j)\n"
            atom_selected_print_list += "         Donor NBO (i)           Acceptor NBO (j)        kcal/mol   a.u.    a.u. \n"
            atom_selected_print_list += " ================================================================================\n"

        #cannot used self.parsed_data_container as it contains the abbreviated e2pert data...
        perturbnrg_list = self.all_data_container[fname][2]
        

        atoms = ['C', 'H', 'N', 'O', 'F', 'S', 'I']
        alkynewarning = False
        display_by_atoms = []
        print_list = ""
        
        pi_NBO = []
        pruned = []

        e2pert_donor_types = []
        e2pert_donor_atoms = []
        e2pert_acceptor_types = []
        e2pert_acceptor_atoms = []
        e2pert_energies = []
        e2pert_brief = []

        collected_groups = []
        if(int(prnt) == 2 or int(prnt)==1):
            if(NBO6 == 1):

                data_list = []
                for i in range(0, len(perturbnrg_list) - 1):
                    curr_orbs = [-1,-1]
                    curr_nums = []
                    try:
                        pruned = perturbnrg_list[i].split()[:-3]
                    except:
                        pass
                    if(pruned):
                        pass
                    else:
                        continue
                    #first orbital is the easiest
                    curr_orbs[1]=pruned[0]
                    pruned_no_brackets = []

                    for item in pruned:

                        if("BD*" in item):
                            pruned_no_brackets.append(item[:-1])
                        if "(" in item or ")" in item: # or "BD" in item or "NL" in item:
                            pass
                        else:
                            pruned_no_brackets.append(item)
                    
                    _tmpatom = []
                    _types = []
                    for j, item in enumerate(pruned_no_brackets):
                        
                        if ("BD" in item or "BD*" in item) :

                            _types.append(2)
                        if ("LP" in item or "RY" in item):
                            _types.append(1)
                        for atom in atoms:
                            if atom in item:
                                if(len(item)==1 and not pruned_no_brackets[j+1][0].isalpha()):
                                    _tmpatom.append(f"{item}{pruned_no_brackets[j+1]}")
                                elif(len(item)>1):
                                    _tmpatom.append(f"{item}")
                    groups = [[-1,-1],[-1,-1]]
                    try:
                        if(_types[0] == 1):
                            _tmpatom.insert(1," ")
                    except:
                        pass
                    try:
                        if(_types[1] == 1):
                            
                            _tmpatom.append(" ")
                    except:
                        pass
                    for j,count in enumerate(_types):  
                        for k in range(count):
                            try:

                                groups[j][k] = _tmpatom[int(f"{j}{k}", 2)] #binary magic!
                            except:
                                pass
                    #cleanup
                    _grp_nums = []
                    _grp = []
                    atom_groups = []
                    for group in groups:
                        for item in group:
                            
                            if(item != -1):
                                # print("here")
                                item.strip("-")
                                _grp = ''.join([char for char in item if char.isdigit()])
                            #print(perturbnrg_list)
                            _grp_nums.append(_grp)
                        atom_groups.append(_grp_nums)
                        _grp_nums =[]
                        _grp = []

                    printvalue = False
                    pi_NBO_temp = []
                    perturb_temp = perturbnrg_list[i].split() #slice_string(perturbnrg_list[i])
                    
                    if("within" in perturb_temp[0]):
                        continue
                    if(perturb_temp[0][0].isalpha()):
                        continue

                    eCap = True
                    if(perturb_temp[1] == 'BD'):
                        if(perturb_temp[3] == '2)' or perturb_temp[3] ==  '3)'):
                            if(float(perturb_temp[-3]) > float(energycap)):
                                            eCap = False
                                            print_list += perturbnrg_list[i]
                                            data_list.append(perturbnrg_list[i])
                                            pi_NBO_temp.append(perturb_temp[0])
                                            pi_NBO_temp.append(perturb_temp[8])
                                            pi_NBO.append(pi_NBO_temp)
                    if(perturb_temp[1] == 'LP'):
                        if(perturb_temp[4] == ' O' and perturb_temp[3] == '2)'):
                            printvalue = True
                        if(perturb_temp[4] == 'N'):
                            printvalue = True
                        if((float(perturb_temp[-3]) > float(energycap)) and (printvalue == True)):
                            eCap = False
                            print_list += perturbnrg_list[i] 
                            data_list.append(perturbnrg_list[i])
                            pi_NBO_temp.append(perturb_temp[0])
                            pi_NBO_temp.append(perturb_temp[6])
                            pi_NBO.append(pi_NBO_temp)
                    if(perturb_temp[3] == '3)'):
                        alkynewarning = True
                    allin = True
                    
                    if(eCap == False):
                        if(atom_numbers):
                            for batch in atom_groups:
                                for atom in batch:
                                    try:
                                        if int(atom) not in atom_numbers:
                                            allin = False
                                    except:
                                        pass
                        else:
                            pass  
                        if(allin):
                            print("printing...")
                            collected_groups.append(atom_groups)
                            atom_selected_print_list += perturbnrg_list[i]
                        else:
                            pass
                # if(atom_numbers):
                #     print(atom_selected_print_list)

                else:
                    pass
                atom_selected_print_list += f"Possible Pi NBO's for NBODel method: \n {pi_NBO}"                  
                if(alkynewarning):
                    atom_selected_print_list += "\nWarning: possible alkyne in molecule, potentially need to use BD ( 3) for conjugated alkyne!\n"
        return collected_groups, atom_selected_print_list

    def nbbp_parse(self, datacontainer):
        #nbbp_raw = (deepcopy(datacontainer[self.features.index("nbbp")])[2:])[:-1]

        nbbp_raw = (datacontainer[2:])[:-1]
 
        with open("output_dump", "w+") as f:
            for item in nbbp_raw:
                f.write(item)
        nbbp_processed = []
        maxval = -1
        _final_line_len = -1
        _diff = 0
        nbbp_squeezed = []
        bond_list = []
        for i,line in enumerate(nbbp_raw):
            line = line.strip('\n')
            if " NBO " in line:
                continue
            elif "------------ -------" in line: #------- ------- ------- ------- ------- ------- -------" in line:
                continue
            elif line == '\n' or line =='':
                continue
            else:
                nbbp_squeezed.append(line)
            
            if int(line.split()[0][:-1])>maxval:
                maxval = int(line.split()[0][:-1]) # FUCKING FUCK ARE YOU FUCKING KIDDING ME FUCK
                _tmp = [line.split()[0], ' '.join(line.split()[1:-8])]  #need to collect bond information
                                                                        #piggyback off the "maxval" computation, where we only iterate through the first nbbp block (need to stop after first block
                                                                        # or atleast before the last block, especially if it has <8 items)
                bond_list.append(_tmp)
        
        #print(nbbp_squeezed)

        #attempt #2 (use better parsing + dataframes)
        #create a blank df with cols = [x for x in range(maxval)], index = [x for x in range(maxval)]
        #iterate through nbbp_squeezed, collect all rows until 'maxval' times (one block) are reached, cut that into a tmp_arr. 
        # each iteration through the nbbp_squeezed 'maxval' times increases counter+1 (keeps track of the block number)
        #iterate through tmp_arr
        #split each row first by "     " (.split('     '))
        #split the arr by default (.split()), iterate through each item, make it a float
        #place 
        
        nbbp_df = pd.DataFrame(index = [int(x) for x in range(maxval)])#, columns = [int(x) for x in range(-1, maxval)])
        
        #print(nbbp_df)
        _tmparr = []
        _block = []
        _nbbp_names = []
        for i, line in enumerate(nbbp_squeezed):
            counter = 0
            if(len(_nbbp_names))<maxval:
                _nbbp_names.append(nbbp_squeezed[i].split('    ')[0]) #risky splitting here....
            _tmparr.extend(nbbp_squeezed[i].split('    ')[1].split())
            
            #print(_tmparr)

            _block.append(_tmparr)
            #print(_nbbp_names, _tmparr)
            #x=input()
            _tmparr = []
            #print(i%maxval)
            if((i+1)%maxval == 0 and i>0):

                _range = [int(((i+1)/maxval)*8)-8,int(((i+1)/maxval)*8)-1-(8-len(_block[0]))]# int(((i+1/282)*8)+(len(_block[0])-1))]

                _idx = [x for x in range(int(_range[0]), int(_range[1])+1)]
                #x=input(_idx)
                _block_df = pd.DataFrame(_block, columns = _idx)# [x for x in range(int(_range[0]), int(_range[1])+1)])
                #print("here", _block_df)
                #print("len")
                #x=input()
                nbbp_df = pd.concat([nbbp_df, _block_df], axis=1)
                _block = []
            #print(_nbbp_names[0])
            #exit()
        #print(_nbbp_names)
        nbbp_df.index = _nbbp_names
        return nbbp_df, _nbbp_names
        
    #print output methods for ease of handling

    def get_atom_label(self, atom, fname):
         return self.parsed_all_data_container[fname][-1][atom]
    
    def get_all_atoms(self, fname):
        nums = [i for i,x in enumerate(self.parsed_all_data_container[fname][-1])]
        labels = self.parsed_all_data_container[fname][-1]
        return nums, labels
    
    #return a list of atom indices based on atom label
    def get_atom_nums(self, labels, fname):
        return [i for i,x in enumerate(self.parsed_all_data_container[fname][-1]) if x in labels]
    
    def get_npa(self, fname):
        _batch_atom_npa = []
        for atoms in self.parsed_all_data_container[fname][0]:
            _batch_atom_npa.append(atoms[1])
        return _batch_atom_npa


    def get_bbmatrix(self, fname):
        #return bbmatrix, bbmatrix_indices
        return self.parsed_all_data_container[fname][1][0] # self.parsed_all_data_container[fname][1][1]

    #return the bbmatrix indexes of a particular file
    def get_bbmatrix_idx(self,fname):
        return self.parsed_all_data_container[fname][1][1]

    #find the atoms, return the labels and their respective locations.
    def get_atoms_in_nbbp(self, atoms, fname): #get all atoms, based on symbol (e.g., "S")
        atoms = [x for i,x in enumerate(self.parsed_all_data_container[fname][1][1]) for atom in atoms if f"{atom}" in x]
        atom_locs = [i for i,x in enumerate(self.parsed_all_data_container[fname][1][1]) for atom in atoms if f"{atom}" in x]
        return atoms, atom_locs

    def fname_process(self, fname):
        return ''.join(fname.split('_')[1:3]) #SaxR_2a_3k_HF_def2-svp_out becomes 2a3k, custom for this file batch

    #recall, alist IS A LIST OF ATOMS!
    #convert list of atom#'s to bbmatrix list by "get_Atoms_in_nbbp()"
    def cross_nbbp(self, alist1, alist2, fname, style): #get the intersection of alist1 and alist2, based on atom#
        #first get the 
        alist1_ids = []
        print("uhm", alist1)
        print(alist2)
        alist1_ids = [f"{self.get_atom_label(alist1[i],fname)} {x+1}" for i,x in enumerate(alist1) if int(alist1[i])<10]
        alist1_ids += [f"{self.get_atom_label(alist1[i], fname)}{x+1}" for i,x in enumerate(alist1) if int(alist1[i])>=10]

        alist2_ids = []
        alist2_ids = [f"{self.get_atom_label(alist2[i],fname)} {x+1}" for i,x in enumerate(alist2) if int(alist2[i])<10]
        alist2_ids += [f"{self.get_atom_label(alist2[i],fname)}{x+1}" for i,x in enumerate(alist2) if int(alist2[i])>=10]

        print(alist1_ids)
        print(alist2_ids)

        cols = alist2_ids
        bb_matrix = self.get_bbmatrix(fname)

        print("1",alist2_ids)
        print("2", alist1_ids)
        
        #for item in alist1_ids: 
        #for each item in alist1/2_ids, get the locations in bb_matrix
        #_alist1 is shortlist
        _alist1_locs = [i for i,x in enumerate(bb_matrix.index) for y in alist1_ids if y in x]
        _alist1_items = [x for i,x in enumerate(bb_matrix.index) for y in alist1_ids if y in x]
        #_alist2 is the longlist
        _alist2_locs = [i for i,x in enumerate(bb_matrix.index) for y in alist2_ids if y in x]
        _alist2_items = [x for i,x in enumerate(bb_matrix.index) for y in alist2_ids if y in x]

        final_arr = []
        _tmp = []
        if(style == "mono"):

            if(len(alist1)==1):
                print("inside")
                _tmp.append(self.fname_process(fname))
                for j, _alist2 in enumerate(_alist2_locs):
                    print(alist1,alist2)
                    print(bb_matrix[_alist1_locs[0]][_alist2])
                    print("here')")
                    _tmp.append(bb_matrix[_alist1_locs[0]][_alist2])
                final_arr.append(_tmp)
                _tmp = []
                return final_arr, _alist2_items, _alist1_items
            else:
                style == "square"
                pass
        if(style == "square"):
            for i,alist1 in enumerate(_alist1_locs):
                _tmp.append(_alist1_items[i])
                for j,alist2 in enumerate(_alist2_locs):
                    #print(alist1,alist2)
                    #print(bb_matrix[alist1][alist2])
                    _tmp.append(bb_matrix[alist1][alist2])
                final_arr.append(_tmp)
                _tmp = []
            return final_arr, _alist2_items, _alist1_items



    #attempt #1 (commented out):
    """
        for i, line in enumerate(nbbp_squeezed):
            _collect = []

            if(int((i/maxval))==int(maxval/8) ): #maxval/8 = 20 blocks, then a remainder
                #if(int(line.split()[0][:-1]) == maxval and len(line.split()) != _final_line_len and _final_line_len>0):
                _tmpfix = []# [float(x) for x in line.split() if float(x) < 1]
                for j,item in enumerate(line.split()):
                    #_tmpfix = []
                    try:
                        _x = float(item)
                        if(_x<2.0):
                            #print(line.split())
                            if(j>1):
                                _tmpfix.append(_x)
                                #print("appended", _tmpfix)
                        #x=input("split line")

                    except:
                        pass
                _collect.extend(_tmpfix)

                #_collect.extend(line.split()[(-8+_diff):])
            elif(float((i/maxval))<int(maxval/8)):
                _collect.extend(line.split()[-8:])

            else:
                _collect.extend(line.split()[-8:])

            nbbp_processed.append(_collect)

        nbbp_final = [[x] for x in range(0, maxval)]
        nbbp_noshape = []

        for i in range(0, maxval):
            for j in range(0, int((maxval/8))):
                nbbp_noshape.extend(nbbp_processed[i + maxval*j])

        if(float(maxval/8))>(int(maxval/8)):
            extr_count = int((float(maxval/8)-int(maxval/8))*8)
            extra = nbbp_processed[ -maxval:]

            for item in extra:
                ##nbbp_noshape.extend(item)

                for it in item:
                    #print("here", it, len(extra))
                    nbbp_noshape.append(it)
            
        nbbp_np = np.array(nbbp_noshape)
        nbbp_data = np.reshape(nbbp_np, [maxval,maxval]).squeeze().tolist()

        return nbbp_data, bond_list
    """
    #def e2pert_parse(self, datacontainer):
    #    pass

    def nbbp(self, line):
        if(self.nbbp_flag[0] in line[0]):
            return True
        if(line[0] == self.nbbp_flag[1] and line[1] == self.nbbp_flag[1]):
            return -1
        else:
            return False
    def e2pert(self, line):
        if(self.e2pert_flag[0] in line[0]):
            return True
        if(line[0] == self.e2pert_flag[1] and line[1] == self.e2pert_flag[1]):
            return -1
        else:
            return False


        
"""    
npa_extr = NBOExtract(software="orca")
npa_extr.read_all()
npa_extr.parse_all()
#npa_extr.npa()
# npa_extr.readtext(npa_extr.files[0])
# npa_extr.read_all()
np1 = np.array(npa_extr.parsed_all_data_container[npa_extr.files[0]][0][0:])
print(np1[:,:1].squeeze().tolist())
#print(np1[:,:])
#print(npa_extr.parsed_all_data_container["SaxR_2a_4w_HF_def2-svp_out"][1][44][44]) # checks the polarization of atom 44
"""
