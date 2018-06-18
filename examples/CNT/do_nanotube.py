from ase.build import nanotube

# nanotube(n, m, length=1, bond=1.42, symbol='C', verbose=False, vacuum=None)
atoms = nanotube(2,2, length=5, bond=1.42, symbol='C', verbose=True, vacuum=5.0)
atoms.set_pbc([True,True,True])
atoms.write('CNT.xsf')

