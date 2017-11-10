from __future__ import print_function
import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import sys

# Read from argument
filename = sys.argv[1]

tree = ET.parse(filename)
root = tree.getroot()

IDX_MESH  = -1
IDX_V_LOC = -1
IDX_NL    = -1
IDX_PSWFC = -1

#
# Loop over children
# search for appropriate index for speficic data
#
Nchild = len(root)
print('List of tags:')
for ic in range(Nchild):
    tag = root[ic].tag
    print(tag)
    #
    if( tag == 'PP_MESH' ):
        IDX_MESH = ic
    elif( tag == 'PP_LOCAL' ):
        IDX_V_LOC = ic
    elif( tag == 'PP_NONLOCAL' ):
        IDX_NL = ic
    elif( tag == 'PP_PSWFC' ):
        IDX_PSWFC = ic


print('')
print('IDX_MESH  = ', IDX_MESH)
print('IDX_V_LOC = ', IDX_V_LOC)
print('IDX_NL    = ', IDX_NL)
print('IDX_PSWFC = ', IDX_PSWFC)
print('')

# root[0] is (and should be) PP_INFO
pp_info = root[0]

# pp_info only contains text
DISPLAY_INFO = False
# Display the content of PP_INFO
if DISPLAY_INFO:
    print(root[0].text)

#
pp_mesh = root[IDX_MESH]
#
# pp_mesh might have attribute
#

# Number of mesh
Nmesh = len(pp_mesh)
IDX_R = 0
IDX_RAB = 0
# Search for pp_r
# It is usually the first mesh, I include it here just in case
# this is not so for the pseudopotential
for i in range(Nmesh):
    if pp_mesh[i].tag == 'PP_R':
        IDX_R = i
    if pp_mesh[i].tag == 'PP_RAB':
        IDX_RAB = i


# Radial grid
pp_r = pp_mesh[IDX_R]
pp_rab = pp_mesh[IDX_RAB]
#print(root[IDX_MESH][0].items())
#Nradial = int( root[IDX_MESH][0].items()[IDX_MESH][1] )
Nradial = int( pp_r.attrib['size'] )
Nradial_ab = int( pp_r.attrib['size'] )
if Nradial != Nradial_ab:
    RuntimeWarning('Nradial and Nradial_ab is different')
r = np.zeros(Nradial)
rab = np.zeros(Nradial)

#IDX_MESH = 2
for ir in range(Nradial):
    r[ir] = float( pp_r.text.split()[ir] )
    rab[ir] = float( pp_rab.text.split()[ir] )


# Local pseudopotential
pp_loc = root[IDX_V_LOC]
# FIXME: check ???
if Nradial != int( pp_loc.attrib['size'] ):
    RuntimeWarning('Nradial and Nradial_ab is different')
V_loc = np.zeros(Nradial)
for ir in range(Nradial):
    V_loc[ir] = float( pp_loc.text.split()[ir] )

"""
#IDX_NL = 5
Nbeta = len(root[IDX_NL]) - 1
beta = np.zeros([Nradial,Nbeta])
beta_r_cut = np.zeros(Nbeta)

for ibeta in range(Nbeta):
    # search for attribute cutoff_radius
    idx = 0
    Lattr = len(root[IDX_NL][ibeta].items())
    for i in range(Lattr):
        if root[IDX_NL][ibeta].items()[i][0] == 'cutoff_radius':
            idx = i
            break
    beta_r_cut[ibeta] = float( root[IDX_NL][ibeta].items()[idx][1] )
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
"""


plt.clf()
plt.plot( r, V_loc, marker='o' )
plt.grid()
plt.xlim(0,8.0)
plt.title(filename)
plt.savefig('V_loc.png', dpi=300)

"""
plt.clf()
for ibeta in range(Nbeta):
    plt.plot( r, beta[:,ibeta], marker='o', label='beta-'+str(ibeta+1))
plt.xlim(0, np.max(beta_r_cut) )
plt.grid()
plt.legend()
plt.title(filename)
plt.savefig('beta_NL.png')

plt.clf()
for iwfc in range(Npswfc):
    plt.plot( r, pswfc[:,iwfc], marker='o', label='pswfc'+pswfc_label[iwfc] )
plt.xlim(0, 5.0 )
plt.grid()
plt.legend()
plt.title(filename)
plt.savefig('pswfc.png')
"""
