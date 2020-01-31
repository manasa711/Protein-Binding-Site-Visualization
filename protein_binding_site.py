#!/usr/bin/env python3

import argparse

def Interface_atoms(pdf_file,thres,chain1,chain2):

    fhand = open(pdf_file)

#creating lists with the coordinates of CA atoms from both the chains

    cds1 = []
    cds2 = []

    for line in fhand:
        line= line.rstrip()
        atom = []
        ch1 = []
        ch2 = []
        if line.startswith("ATOM"):
            atom = line.split()
            if atom[4] == chain1:
                if atom[2] == 'CA': #only considering the C-Alpha atoms for the computation
                    x = atom[6]
                    y = atom[7]
                    z = atom[8]
                    aa_no = atom[5]
                    aa_name = atom[3]
                    ch1.append(x)
                    ch1.append(y)
                    ch1.append(z)
                    ch1.append(aa_no)
                    ch1.append(aa_name)
                    cds1.append(ch1)
            elif atom[4] == chain2:
                if atom[2] == 'CA':
                    x = atom[6]
                    y = atom[7]
                    z = atom[8]
                    aa_no = atom[5]
                    aa_name = atom[3]
                    ch2.append(x)
                    ch2.append(y)
                    ch2.append(z)
                    ch2.append(aa_no)
                    ch2.append(aa_name)
                    cds2.append(ch2)

    #calculating Euclidean Distance between CA atoms of chain 1 and CA atoms of chain2

    euc_vals_1 = [] #list with interface atoms from chain 1
    euc_vals_2 = [] #list with interface atoms from chain 2

    for i in cds1:
        for j in cds2:
            x1 = float(i[0])
            y1 = float(i[1])
            z1 = float(i[2])
            x2 = float(j[0])
            y2 = float(j[1])
            z2 = float(j[2])
            e = ((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2)
            euc = e**0.5
            if euc <= thres:
                op = chain1+":"+str(i[4])+"("+str(i[3])+") interacts with "+chain2+":"+str(j[4])+"("+str(j[3])+")"
                print(op) #prints out the atoms from chain1 which interact with the atoms from chain2
                euc_vals_1.append(int(i[3]))
                euc_vals_2.append(int(j[3]))

    #determining the number of interface CA atoms lying on secondary (heilces or beta-sheets) over all interfaced CA atoms

    helix_no_1 = []
    helix_no_2 = []
    sheet_no_1 = []
    sheet_no_2 = []

    for line in fhand:
        line= line.rstrip()
        helix = []
        sheet = []

        if line.startswith("HELIX"):
            helix = line.split()
            st_h = int(helix[5])
            te_h = int(helix[8])
            if helix[4] == chain1:
                for i in euc_vals_1:
                    if int(i) in range(st_h,te_h+1):
                        if int(i) not in helix_no_1:
                            helix_no_1.append(int(i))
            elif helix[4] == chain2:
                for i in euc_vals_2:
                    if int(i) in range(st_h,te_h+1):
                        if int(i) not in helix_no_2:
                            helix_no_2.append(int(i))

        if line.startswith("SHEET"):
            sheet = line.split()
            st_b = int(sheet[6])
            te_b = int(sheet[9])
            #print(st_b,te_b)
            if sheet[5] == chain1:
                for i in euc_vals_1:
                    if int(i) in range(st_b,te_b+1):
                        if int(i) not in sheet_no_1:
                            sheet_no_1.append(int(i))
            elif sheet[5] == chain2:
                for i in euc_vals_2:
                    if int(i) in range(st_b,te_b+1):
                        if int(i) not in sheet_no_2:
                            sheet_no_2.append(int(i))



    fhand.close()

    interface_atoms_1 = [] #interface atoms from chain1
    interface_atoms_2 = [] #interface atoms from chain2

    #removing repeated CA atoms
    for i in euc_vals_1:
        if int(i) not in interface_atoms_1:
            interface_atoms_1.append(int(i))

    for i in euc_vals_2:
        if int(i) not in interface_atoms_2:
            interface_atoms_2.append(int(i))

    #printing
    print("\n")
    print("Chain Name:", chain1)
    print(str(len(helix_no_1))+"/"+str(len(interface_atoms_1))+" of the interface atoms lying on alpha helices.")
    print(str(len(sheet_no_1))+"/"+str(len(interface_atoms_1))+" of the interface atoms lying on beta sheets.")
    print("\n")
    print("Chain Name:", chain2)
    print(str(len(helix_no_2))+"/"+str(len(interface_atoms_2))+" of the interface atoms lying on alpha helices.")
    print(str(len(sheet_no_2))+"/"+str(len(interface_atoms_2))+" of the interface atoms lying on beta sheets.")
    print("\n")

    interface_atoms_1.sort()
    interface_atoms_2.sort()

    print(chain1)
    for i in range(1,len(interface_atoms_1)):
        r = interface_atoms_1[i]
        g = interface_atoms_1[i-1]
        dis = int(r) - int(g)
        for j in cds1:
            if int(r) == int(j[3]):
                aa2 = j[4]
        for j in cds1:
            if int(g) == int(j[3]):
                aa1 = j[4]
        print(str(aa1)+":"+" closest "+str(aa2)+" at distance "+str(dis))

    print("\n")

    print(chain2)
    for i in range(1,len(interface_atoms_2)):
        r = interface_atoms_2[i]
        g = interface_atoms_2[i-1]
        dis = int(r) - int(g)
        for j in cds2:
            if int(r) == int(j[3]):
                aa2 = j[4]
        for j in cds2:
            if int(g) == int(j[3]):
                aa1 = j[4]
        print(str(aa1)+":"+" closest "+str(aa2)+" at distance "+str(dis))

def main():

    #Argparse code
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help = "This argument is for the protein PDB file")
    parser.add_argument("-t", help = "This argument is for specifying the threshold value")
    parser.add_argument("-c1", help = "The argument is for specifying the first chain name. Usually is given as a single capital letter")
    parser.add_argument("-c2", help = "This argument is for specifying the second chain name. (Single Capital Letter)")
    args = parser.parse_args()

    if args.f:
        print("The input file given is: "+str(args.f))
    if args.t:
        print("The threshold value selected is: "+str(args.t))
    if args.c1:
        print("Specified First Chain: "+str(args.c1))
    if args.c2:
        print("Specified Second Chain: "+str(args.c2))

    print("\n")

    pdb_file = args.f #Variable for storing the pdb file name or path
    thres = int(args.t) #Variable for storing the specified threshold value
    chain1 = args.c1 #Variable for storing the first chain name
    chain2 = args.c2 #Variable for storing the second chain name

    Interface_atoms(pdb_file,thres,chain1,chain2)

if __name__ == "__main__":
    main()
