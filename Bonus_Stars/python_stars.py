import numpy as np

Array = np.loadtxt(
    "GalaxyFromIllustrisTNG50_Stars_Subhalo521803.txt"
)
Pos = Array[:,0:3]
Mass = Array[:,3]
InitialMass = Array[:,4]
RedshiftFormed = Array[:,5]

x = Array[:, 0]
y = Array[:, 1]
z = Array[:, 2]


tensor = np.zeros( (3,3) )
tensor[0,0] = np.sum(Mass * ( y * y + z * z ))
tensor[1,1] = np.sum(Mass * ( x * x + z * z ))
tensor[2,2] = np.sum(Mass * ( x * x + y * y ))
tensor[0,1] = -np.sum(Mass * x * y )
tensor[1,0] = tensor[0,1]
tensor[0,2] = -np.sum(Mass * x * z )
tensor[2,0] = tensor[0,2]
tensor[1,2] = -np.sum(Mass * y * z )
tensor[2,1] = tensor[1,2]

eigval , eigvec = np.linalg.eig( tensor )

xdir = eigvec[:, 0]
ydir = eigvec[:, 1]
zdir = np.cross(xdir, ydir)

print(xdir)
print(ydir)
print(zdir)