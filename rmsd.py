from pyrosetta import *
import numpy as np
init('-mute all')
import glob
import os
import argparse

try:
    os.mkdir('results')
except FileExistsError:
    pass

parser = argparse.ArgumentParser()
parser.add_argument( "-inputs", type=str, help='directory with input hotloops' )
args = parser.parse_args()

path = args.inputs
hotloop_pdbs = glob.glob(path + '/*.pdb')

#makes a pairwise distance matrix for 4 consecutive CA
def CA_matrix(pose, fragment_array):
    #initialize CA distance matrix 
    CA_matrix = np.zeros((4,4))
    
    #compute distance between each CA, ex btwn CA1 and CA2 or between CA1 and CA4
    for i in range(0,4):
        CA_atom1 = pose.residue(fragment_array[i]).xyz("CA")
        for j in range(0,4):
            CA_atom2 = pose.residue(fragment_array[j]).xyz("CA")
            CA_matrix[i,j] = CA_atom1.distance(CA_atom2)
    return CA_matrix

#makes a dihedral matrix with phi (col1) and psi(col2) for each residue
def dihedral_vec(pose, fragment_array):
    #initialize dihedral matrix 
    DH_vec = np.zeros((4,2))
    #calculate phi and psi for each residue
    for i in range(0,4):
        DH_vec[i,0] = pose.phi(fragment_array[i])
        DH_vec[i,1] = pose.psi(fragment_array[i])
        
    return DH_vec

#calculates rmsd 
def calculate_rmsd(matrix1, matrix2):
    #squared diff
    squared_diff = (matrix1 - matrix2) ** 2
    #mean
    mean_squared_diff = np.mean(squared_diff)
    #root!
    rmsd = np.sqrt(mean_squared_diff)
    return rmsd

#permute all scaffolds to give all possible 4 AA residue index fragments 
def scaffold_fragment(pose):
    length = pose.total_residue()
    num_list = np.arange(1, length + 1)
    fragments = []
    for i in range(length):
        #this term is for over the "terminal" fragments ex return for length 10 cycpep [10, 1, 2, 3]
        if i > length - 4: 
            fragment = np.concatenate((num_list[i:], num_list[:(i - (length - 4))]))
        else:
            fragment = num_list[i:i+4]
        fragments.append(fragment)
    return fragments

#makes all possible 4 AA residue index fragments from a hotloop 
def hotloop_fragment(pose):
    hl_fragments = []
    length = pose.total_residue()
    num_list = np.arange(1, length + 1)
    for i in range(length - 3):
        fragment = num_list[i:i+4]
        hl_fragments.append(fragment)
    return hl_fragments

def rmsd_compare(hot_loop_pose, hotloop_fragments, scaffold_pose, scaffold_fragments, match):
    for hl_fragment in hotloop_fragments:
        hot_loop_matrix = CA_matrix(hot_loop_pose, hl_fragment)
        hotloop_DH = dihedral_vec(hot_loop_pose, hl_fragment)
        
        for scaff_fragment in scaffold_fragments:
            scaffold_matrix = CA_matrix(scaffold_pose, scaff_fragment)
            scaffold_DH = dihedral_vec(scaffold_pose, scaff_fragment)
            CA_rmsd = calculate_rmsd(scaffold_matrix, hot_loop_matrix)
            DH_rmsd = calculate_rmsd(scaffold_DH, hotloop_DH)
            if CA_rmsd < 1 and DH_rmsd < 10: 
                match = True 
                break
        if match == True:
            break
    return scaffold_pose, scaff_fragment, hl_fragment, CA_rmsd, DH_rmsd, match 

scaffolds = glob.glob("PATH TO SCAFFOLDS")

for hotloop_pdb in hotloop_pdbs:
#    print(hotloop_pdb)
    hot_loop_pose = pose_from_pdb(hotloop_pdb)
    hotloop_fragments = hotloop_fragment(hot_loop_pose)
    for pdb in scaffolds: 
        match = False
        scaffold_pose = pose_from_pdb(pdb)
        scaffold_fragments = scaffold_fragment(scaffold_pose)

        scaffold_pose, scaff_fragment, hl_fragment, CA_rmsd, DH_rmsd, match = rmsd_compare(hot_loop_pose, hotloop_fragments, scaffold_pose, scaffold_fragments, match)

        if match == True:
            name = hotloop_pdb.split("/")[-1].replace(".pdb", "")
            file = open(f"results/{name}.txt", "w")
            lines = [f"scaffold: {pdb}\n" ,f"hotloop: {hotloop_pdb}\n" ,f"scaffold fragment: {scaff_fragment}\n" , f"hotloop fragment: {hl_fragment}\n", f"{CA_rmsd}, {DH_rmsd}"]
            file.writelines(lines)
            file.close()
            break
    
