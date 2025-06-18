import molli as ml

#generate the library and initialize the directory structure
libini  = ml.librarygenerate.SubstituentCreate()
libini.make_substituents()
libini.write_xyzs()
