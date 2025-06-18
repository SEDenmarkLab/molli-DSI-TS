from nbotools import nboparse as ne

"""
Structure of nboparse.py

NBOExtract:
    handles bulk file processing of output files 
    
    only orca _specifically_ handled right now,
    though it would not take much effort to incorporate other output files such as gaussian.
    Only the "flags" need to be modified to what is observed in Gaussian output files. 
    The structure of the NBO output should mostly remain the same, so the rest of the parser should
    work without too much effort.

    to enabled rapid reading/writing of large numbers of output files, files are read 

"""



#initialize the data and parse
nbo_extr = ne.NBOExtract(software="orca", endswith="out")
nbo_extr.read_all()
nbo_extr.parse_all()

#print list of files available
print(nbo_extr.files)

##### bulk NPA processing
#get a list of the atom numbers and labels (extracted in order from the npa output)
file0_atom_nums, file0_labels = nbo_extr.get_all_atoms(nbo_extr.files[0])
print(file0_atom_nums)
print(nbo_extr.get_npa(nbo_extr.files[0]))

##### bulk nbbp processing 
#print out the list of nbbp matrix bonds for reference:
print(nbo_extr.get_bbmatrix_idx(nbo_extr.files[0]))
#get the entire nbbp matrix data:
print(nbo_extr.get_bbmatrix(nbo_extr.files[0]))
#get all atoms in nbbp containing sulfur:
print(nbo_extr.get_atoms_in_nbbp(atoms="S", fname=nbo_extr.files[0]))
#get specific cross terms in nbbp based on atom numbers:
print(nbo_extr.cross_nbbp([1,2,3], [4,5,6], nbo_extr.files[0], style = "square"))

#print abbreviated E2PERT data:
print(nbo_extr.parsed_all_data_container)
atoms, outputs = nbo_extr.get_e2pert_data(nbo_extr.files[0], [x for x in range(50)], 10.0, header=True)

print(atoms)
print(outputs)
