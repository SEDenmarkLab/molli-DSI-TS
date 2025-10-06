import molli as ml

#generate the library and initialize the directory structure
libini  = ml.librarygenerate.SubstituentCreate(substituentfile = "substituents.csv", base_struct = ["Clc1ccc(X)cc1", "Clc1cc(X)ccc1"])
libini.make_substituents()
libini.write_xyzs()
