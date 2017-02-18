import matplotlib.pyplot as plt
import numpy as np
import math

Nbasis = np.array( [4913, 6859, 9261, 12167, 15625, 19683, 24389, 29791,
	35937, 42875, 50653, 59319, 68921] )

Error = np.array( [0.3858029756E+09, 0.9376990666E+07, 0.6350023804E+06, 0.9420802397E+04,
                   0.2416418199E+03, 0.2922268599E+01, 0.6443719832E-01, 0.2393797987E-02,
									 0.1627382979E-04, 0.4256783261E-04, 0.3225516720E-06, 0.6630791000E-07,
                   0.1268224425E-07] )

plt.clf()
plt.plot( Nbasis, np.log10(Error), marker='o' )
plt.savefig( 'conv_cg.png', dpi=300 )

      
