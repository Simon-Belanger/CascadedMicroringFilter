import numpy as np
from components.MRF import MRF
import matplotlib.pyplot as plt

# Wavelength span for the sweep
lambda_0 = np.linspace(1535, 1550, 1000)*1e-9

# Build the 2 directionnal couplers
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flattop2
k_vec = flattop2(0.3)
couplers = [virtual_DC(k_vec[0],1), virtual_DC(k_vec[1],1), virtual_DC(k_vec[2],1)]


# Create a MRF filter with no thermal crosstalk
mrf = MRF( name='', num_rings=2, R=2.5e-6, couplers=couplers, alpha_wg=3., crosstalk_coeff=[1, 0.5, 0])

# Sweep
#mr.sweep(lambda_0, 1, 0, True)

# Scramble filter
#mrf.manufacturing(np.pi)

# Tuning Map
matsize = 200
V1 = np.linspace(0, np.pi, matsize)
V2 = np.linspace(0, np.pi, matsize)
DDM = np.zeros((matsize,matsize))
for i in range(1,matsize):
    for j in range(1,matsize):
        DDM[i, j] = mrf.test_MRF(1560e-9, [V1[i], V2[j]])
# Plot the tuning map
f = plt.figure()
plt.imshow(DDM, extent=[min(V1),max(V1),min(V2),max(V2)])
plt.xlabel("$\phi_1$ [rad]")
plt.ylabel("$\phi_2$ [rad]")
plt.colorbar()
plt.show()
f.savefig("0_5.pdf", bbox_inches='tight')