#----------------Wannier2ViDES Updated------------------#
# Version 2. Added correction to atom pairs. Now code considers wrapping around unit cell limit nx, ny, nz.
# Version 3. Added Cutoff parameter
# Version 4. Added sorting feature to list site energy first and then coupling coefficients
# Version 5. Tweaked cutoff distance criteria to see acutal WF centers.
#
# This code produces Hamiiltonian array for NanoTCAD ViDES code.

# Input
############
# Seedname: The seedname of your Wannier calculation.
# Wannier calculation requires write_hr = .true. and  write_xyz = .true.

# expand_grid : How many unit cell you wish to expand the system.
# nx, ny, nz each specify how many times unit cell is going to be repeated.
# The total number of atoms is nx * ny * nz * (unit cell atoms).
# n, Nc parameter used in Hamiltonian command for ViDES is n = nx * ny * (unit cell atoms), Nc = nz

# cutoff : The cutoff distance between unit cells to write Hamiltonian element in.
# If the cutoff distance < unit cell distance, Hamiltonian matrix element is not written.
############

# Output
############
# Wannier2ViDES.H : The Hamiltonian matrix used to create Hamiltonian in ViDES code.
# Wannier2ViDES.xyz : The coordinate of each atom inside the system.
# n, Nc : n, Nc you should write to create the Hamiltonian in ViDES code.
############
#------------------------------------------------------------#

from math import sqrt

seedname = input('Enter a seedname for win file : ')
expand_grid = input("How much do you wish to expand the system? (n_x n_y n_z): ")
cutoff = input("Enter a cutoff distance : ")

expand_coeff = expand_grid.split()
nx = int(expand_coeff[0]); ny = int(expand_coeff[1]); nz = int(expand_coeff[2]);
cutoff = float(cutoff);

#-------Relevant Information Extraction------#
h = open('%s_hr.dat' % seedname, 'r');
H = open('Wannier2ViDES.H', 'w')

h.readline()
data = h.readline()
m = int(data.strip()) # Number of orbitals (or MLWF) inside a unit cell
tot_atom = nx * ny * nz * m # Total number of atoms

win = open('%s.win' % seedname, 'r');
wxyz = open("%s_centres.xyz" % seedname, 'r');

win_data = win.readlines();
wxyz_data = wxyz.readlines();
win.close(); wxyz.close();
atoms_cart = [0] * m # Initialize atoms coordinate
alat = [0] * 3 # Initialize the unit cell lattice vector information (Only orthogonal P allowed)

# Lattice vector information
for i in range(len(win_data)):
    if 'Begin Unit_Cell_Cart' in win_data[i]: # Get the information on unit cell lattice vectors
        for j in range(3):
            data_lst = win_data[i+j+1].split()
            # alat is a length 3 list with length of lattice vectors
            alat[j] = float(data_lst[j])

# Wannier function center information
for i in range(m):
    data_lst = wxyz_data[2+i].split();
    # atoms_cart is a list of list containing coordinates of each WF inside the unit cell
    atoms_cart[i] = data_lst[1:]
    atoms_cart[i] = [float(x) for x in atoms_cart[i]]
#------Information Extraction Done--------#

#------Print out relevant information to stdout------#
print("\n Wannier2ViDES code, Ver. 5");
print("\n Input:");
print("seedname = %s" % seedname);
print("expand_grid = %d %d %d" % (nx, ny, nz));
print("cutoff = %f" % cutoff);
print("\n Information:")
print("Total number of WFs inside Unit cell = %d" % m)
print("Total number of WFs = %d" % tot_atom)
print("n = %d, Nc = %d" % (nx * ny * m, nz)); # Print the relevant information for ViDES code input
print("System size = %f %f %f" %(nx * alat[0], ny * alat[1], nz * alat[2]));
print("\n Output:");
print("Hamiltonian for ViDES code written to : Wannier2ViDES.H");
print("xyz coordinate for atoms in the system written to : Wannier2ViDES.xyz");
#------Print out done-------#

for i in range(2) : data = h.readline() # Read lines until we reach the fourth line, k-degenerate points
while len(data.split()) == 15: # Ingore lines with degenerate k-points
    data = h.readline();
data = h.readline(); # Now we are at the information of the Hamiltonian

H.write("{:>7}{:>7}{:>11}{:>11}\n".format(1,0,0.000000,0.000000)); # Start of Hamiltonian input file
site_E = []; coupling = [];
while len(data) != 0:
    data_lst = data.split();
    for i in range(len(data_lst)): # Change data type to int or float
        if i < 5:
            data_lst[i] = int(data_lst[i]);
        #else:
        #    data_lst[i] = float(data_lst[i])
    # Ignore rows with negative translation
    if (data_lst[0] < 0) or (data_lst[1] < 0) or (data_lst[2] < 0):
        data = h.readline();
        continue
    # To search for the correct atom, compute how many theh number should be translated
    # We list the atoms inside a unit cell, than unit cells translated by x, then x, y, and then x,y,z
    # x tick up : +m, y tick up : +nx*m, z tick up : +ny*nx*m
    tr = data_lst[0] * m + data_lst[1] * nx * m + data_lst[2] * nx * ny * m
    if data_lst[3] > data_lst[4] + tr:
        data = h.readline();
        continue # Frist atom index has to be smaller than the second one

    # Check cutoff distance
    coord1 = atoms_cart[data_lst[3]-1]; # Coordinate of atom 1
    coord2 = atoms_cart[data_lst[4]-1]; # Initial coordinate of atom 2
    coord2 = [coord2[0] + data_lst[0] * alat[0], coord2[1] + data_lst[1] * alat[1],\
            coord2[2] + data_lst[2] * alat[2]]; # Atom translated by unit cell
    dis = sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2);
    if dis > cutoff:
        data = h.readline();
        continue 

    two_units = [[0,0,0], [data_lst[0], data_lst[1], data_lst[2]]] # Two unit cells where we look at atoms
    # Abuse the translational symmetry of the system to obtain H for the expanded system
    while two_units[1][2] < nz:
        # Reset y direction of unit cell
        two_units[0][1] = 0; two_units[1][1] = data_lst[1];
        while two_units[1][1] < ny:
            # Reset x direction of unit cell
            two_units[0][0] = 0; two_units[1][0] = data_lst[0];
            while two_units[1][0] < nx:
                a1 = data_lst[3]; a2 = data_lst[4]; # Initial atom number
                # Translate the atom by the real unit cell
                a1 = a1 + m*two_units[0][0] + nx*m*two_units[0][1] + nx*ny*m*two_units[0][2];
                a2 = a2 + m*two_units[1][0] + nx*m*two_units[1][1] + nx*ny*m*two_units[1][2];
                H_entry = [a1, a2, data_lst[5], data_lst[6]];
                if a1 == a2: # Site energy Hamiltonian entry
                    site_E.append(H_entry);
                else: # Coupling coefficient Hamiltonian entry
                    coupling.append(H_entry);
                # Now Translate the unit cell in x direction
                two_units[0][0] = two_units[0][0] + 1;
                two_units[1][0] = two_units[1][0] + 1;
            # Translate the unit cell in y direction
            two_units[0][1] = two_units[0][1] + 1;
            two_units[1][1] = two_units[1][1] + 1;
        # Translate the unit cell in z direction
        two_units[0][2] = two_units[0][2] + 1;
        two_units[1][2] = two_units[1][2] + 1;
    data = h.readline();

# Sort site energy and coupling terms and concatenate them
site_E.sort(); coupling.sort();
hamiltonian_out = site_E + coupling;
# Now write H file
for c in hamiltonian_out:
    ham_out = '{:>7}{:>7}{:>11}{:>11}\n'.format(c[0], c[1], c[2], c[3]);
    H.write(ham_out);

h.close(); H.close();

xyz = open('Wannier2ViDES.xyz', 'w');

# Now write the atoms coordinate
for j1 in range(nz):
    for j2 in range(ny):
        for j3 in range(nx):
            for i in range(m):
                atoms_coord = [0] * 3
                # Now calculate the atom coordinate of the translated cells
                atoms_coord[0] = atoms_cart[i][0] + j3 * alat[0]
                atoms_coord[1] = atoms_cart[i][1] + j2 * alat[1]
                atoms_coord[2] = atoms_cart[i][2] + j1 * alat[2]
                coord_out = "%f \t %f \t %f \n" % (atoms_coord[0], atoms_coord[1], atoms_coord[2])
                xyz.write(coord_out);
xyz.close();
