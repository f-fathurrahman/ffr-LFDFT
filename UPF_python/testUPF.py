import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import sys

#tree = ET.parse('C.pw-mt_fhi.UPF')
tree = ET.parse(sys.argv[1])
root = tree.getroot()

IDX_MESH = -1
IDX_V_LOC = -1
IDX_NL = -1
IDX_PSWFC = -1

Nchild = len(root)
print('List of tags:')
for ic in range(Nchild):
    tag = root[ic].tag
    if( tag == 'PP_MESH' ):
        IDX_MESH = ic
    elif( tag == 'PP_LOCAL' ):
        IDX_V_LOC = ic
    elif( tag == 'PP_NONLOCAL' ):
        IDX_NL = ic
    elif( tag == 'PP_PSWFC' ):
        IDX_PSWFC = ic
    print tag

print('')
print('IDX_MESH  = ', IDX_MESH)
print('IDX_V_LOC = ', IDX_V_LOC)
print('IDX_NL    = ', IDX_NL)
print('IDX_PSWFC = ', IDX_PSWFC)
print('')

# root[0] is PP_INFO
# Display the content of PP_INFO
#print(root[0].text)

# Setup radial grid
print(root[IDX_MESH][0].items())
Nradial = int( root[IDX_MESH][0].items()[IDX_MESH][1] )
r = np.zeros(Nradial)
rab = np.zeros(Nradial)

#IDX_MESH = 2
for ir in range(Nradial):
    r[ir] = float( root[IDX_MESH][0].text.split()[ir] )
    rab[ir] = float( root[IDX_MESH][1].text.split()[ir] )

#IDX_V_LOC = 3
V_loc = np.zeros(Nradial)
for ir in range(Nradial):
    V_loc[ir] = float( root[IDX_V_LOC].text.split()[ir] )

#IDX_NL = 5
Nbeta = len(root[IDX_NL]) - 1
beta = np.zeros([Nradial,Nbeta])
beta_r_cut = np.zeros(Nbeta)

for ibeta in range(Nbeta):
    beta_r_cut[ibeta] = float( root[IDX_NL][ibeta].items()[4][1] )
    print("ibeta = %d, r_cut = %f" % (ibeta, beta_r_cut[ibeta]) )
    for ir in range(Nradial):
        beta[ir,ibeta] = float( root[IDX_NL][ibeta].text.split()[ir] )

Npswfc = len(root[IDX_PSWFC])
pswfc = np.zeros( [Nradial,Npswfc] )
pswfc_label = []
for iwfc in range(Npswfc):
    pswfc_label.append( root[IDX_PSWFC][iwfc].items()[3][1] )
    print(pswfc_label[iwfc])
    for ir in range(Nradial):
        pswfc[ir,iwfc] = float( root[IDX_PSWFC][iwfc].text.split()[ir] )


plt.clf()
plt.plot( r, V_loc, marker='o' )
plt.grid()
plt.xlim(0,10.0)
plt.savefig('V_loc.png', dpi=300)

plt.clf()
for ibeta in range(Nbeta):
    plt.plot( r, beta[:,ibeta], marker='o', label='beta-'+str(ibeta+1))
plt.xlim(0, np.max(beta_r_cut) )
plt.grid()
plt.legend()
plt.savefig('beta_NL.png')

plt.clf()
for iwfc in range(Npswfc):
    plt.plot( r, pswfc[:,iwfc], marker='o', label=pswfc_label[iwfc] )
plt.xlim(0, 5.0 )
plt.grid()
plt.legend()
plt.savefig('pswfc.png')
