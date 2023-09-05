#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from pltinteractivelegend import InteractiveLegend as ileg
from tqdm import tqdm

bin_width = 5 # degrees
min = 40
max = 180
bins = np.arange(min,max+bin_width,bin_width)
def histo_line(data):
    histo_height, bin_edges = np.histogram(data, bins=bins, density=True)
    bin_middle = np.diff(bin_edges)/2
    bin_middle = bin_middle + bin_edges[:len(bin_middle)]
    return bin_middle, histo_height

def read_angles(filepath):
    with open(filepath,'r') as f:
        angles = []
        for i,line in enumerate(f):
            if not line:
                print(f'Line {i} is empty.')
                angles.append([])
                continue
            else:
                angles.append([float(x) for x in line.split()])
                
    all_angles = np.array([item for sublist in angles for item in sublist])
    avg_measures = len(all_angles) / len(angles) # avg number of angle measurements per frame
    return all_angles, avg_measures 

# get hydrophobin sequence
with open('sequence.txt', 'r') as file:
    seq = file.read().replace('\n', '')

script_dir = os.path.dirname(__file__) # absolute dir the script is in

# fig, ax = plt.subplots(figsize=(8,7))
# ### PLOT THE BULK WATER 3-BODY ANGLE DISTRIBUTION
# bulk_file = f'{script_dir}/bulk_angles.txt'
# bulk_angles, _ = read_angles(bulk_file)
# bins, bulk_distro = histo_line(bulk_angles)
# plt.plot(bins,bulk_distro,label='bulk',lw=6,color='black')
# plt.ylim(bottom=0,top=0.018)
# plt.xlim(left=40,right=180)
# plt.xlabel(r'Water Triplet Angle ($\theta$)',fontsize=15)
# plt.ylabel(r'$P(\theta)$',fontsize=15)

#%% LOAD AND SAVE DATA
# # get sidechain data
# sidechain_distros = []
# sidechain_names = []
# avg_sidechain_angles_list = []
# sidechain_counter = 0
# print('Processing sidechain data...')
# for resid in tqdm(np.arange(1,71)):
#     # print(f'Looking at residue {resid}')
#     sidechain_file = f'{script_dir}/sidechains/sc_res{resid}_angles.txt'
#     if not os.path.exists(sidechain_file):
#         continue
    
#     sidechain_names.append(f'{seq[resid-1]}{resid}_sidechain')
#     sidechain_angles, avg_sidechain_angles = read_angles(sidechain_file)
#     _, histo = histo_line(sidechain_angles)
#     sidechain_distros.append(histo)
#     avg_sidechain_angles_list.append(avg_sidechain_angles)
#     sidechain_counter += 1

# # turn into DataFrame
# avg_sidechain_angles_list = np.array(avg_sidechain_angles_list).reshape(-1, 1)
# sidechain_data = np.hstack((np.array(sidechain_distros), avg_sidechain_angles_list))
# sidechain_df = pd.DataFrame(data=sidechain_data, index=sidechain_names, columns=np.append(bins, 'avg_sidechain_angles'))
# sidechain_df.to_csv('sidechain_data.csv')

# #%% get backbone data
# backbone_distros = []
# backbone_names = []
# avg_backbone_angles_list = []
# backbone_counter = 0
# print('Processing backbone data...')
# for resid in tqdm(np.arange(1,71)):
#     backbone_file = f'{script_dir}/backbone/bb_res{resid}_angles.txt'
#     if not os.path.exists(backbone_file):
#         continue
    
#     backbone_names.append(f'{seq[resid-1]}{resid}_backbone')
#     backbone_angles,avg_backbone_angles = read_angles(backbone_file)
#     _, histo = histo_line(backbone_angles)
#     backbone_distros.append(histo)
#     avg_backbone_angles_list.append(avg_backbone_angles)
#     backbone_counter += 1

# # turn into DataFrame
# avg_backbone_angles_list = np.array(avg_backbone_angles_list).reshape(-1, 1)
# backbone_data = np.hstack((np.array(backbone_distros), avg_backbone_angles_list))
# backbone_df = pd.DataFrame(data=backbone_data, index=backbone_names, columns=np.append(bins, 'avg_backbone_angles'))
# backbone_df.to_csv('backbone_data.csv')

#%% get residue data
residue_distros = []
residue_names = []
avg_residue_angles_list = []
residue_counter = 0
print('Processing residue data...')
for resid in tqdm(np.arange(333,526+1)):
    residue_file = f'{script_dir}/angles/cv_res{resid}_angles.txt'
    if not os.path.exists(residue_file):
        continue
    
    residue_names.append(f'{seq[resid-333]}{resid}')
    residue_angles,avg_residue_angles = read_angles(residue_file)
    bin_mids, histo = histo_line(residue_angles)
    residue_distros.append(histo)
    avg_residue_angles_list.append(avg_residue_angles)
    residue_counter += 1

# turn into DataFrame
avg_residue_angles_list = np.array(avg_residue_angles_list).reshape(-1, 1)
residue_data = np.hstack((np.array(residue_distros), avg_residue_angles_list))
residue_df = pd.DataFrame(data=residue_data, index=residue_names, columns=np.append(bin_mids, 'avg_residue_angles'))
residue_df.to_csv('residue_data.csv')


#%% PLOT 3-BODY DISTRIBUTION
# # for sidechains
# for i,distro in enumerate(sidechain_distros):
#     if avg_sidechain_angles_list[i] > 1:
#         plt.plot(bins,distro,label=sidechain_names[i])
# plt.plot(bins,bulk_distro,label='bulk',lw=6,color='black')
# plt.ylim(bottom=0,top=0.02)
# plt.xlim(left=40,right=180)
# plt.title('Hydrophobin: Side Chains',fontsize=20)
# plt.xlabel(r'Water Triplet Angle ($\theta$)',fontsize=15)
# plt.ylabel(r'$P(\theta)$',fontsize=15)
# plt.ylim(bottom=0)
# plt.legend(fontsize=4,ncol=4)
# leg = ileg()
# plt.show()

#%%
# # for backbone
# for i,distro in enumerate(backbone_distros):
#     if avg_backbone_angles_list[i] > 1:
#         plt.plot(bins,distro,label=backbone_names[i])
# plt.plot(bins,bulk_distro,label='bulk',lw=6,color='black')
# plt.ylim(bottom=0,top=0.02)
# plt.xlim(left=40,right=180)
# plt.title('Hydrophobin: Backbone',fontsize=20)
# plt.xlabel(r'Water Triplet Angle ($\theta$)',fontsize=15)
# plt.ylabel(r'$P(\theta)$',fontsize=15)
# plt.ylim(bottom=0)
# plt.legend(fontsize=4,ncol=4)
# leg = ileg()
# plt.show()

#%% for full residue
for i,distro in enumerate(residue_distros):
    if avg_residue_angles_list[i] > 1:
        plt.plot(bin_mids,distro,label=residue_names[i])
# plt.plot(bins,bulk_distro,label='bulk',lw=6,color='black')
plt.ylim(bottom=0,top=0.02)
plt.xlim(left=40,right=180)
plt.title('SARS2: Full Residues',fontsize=20)
plt.xlabel(r'Water Triplet Angle ($\theta$)',fontsize=15)
plt.ylabel(r'$P(\theta)$',fontsize=15)
plt.ylim(bottom=0)
plt.legend(fontsize=7,ncol=4)
leg = ileg()
plt.show()

# #%% PLOT DIFFERENCE BETWEEN EACH DISTRIBUTION AND THE BULK DISTRIBUTION
# distros = np.array(distros)
# d_distros = distros - np.array(bulk_distro)

# fig2, ax2 = plt.subplots(figsize=(8,7))
# for i,distro in enumerate(d_distros):
#     plt.plot(bins,distro,label=distro_names[i])


# # for i,distro in enumerate(d_distros):
# #     xnew = np.linspace(40, 180, 80) 
# #     spl = make_interp_spline(bins, distro, k=3)  # type: BSpline
# #     distro_smooth = spl(xnew)
# #     plt.plot(xnew,distro_smooth,label=distro_names[i])
    
# plt.hlines(y = 0, xmin = 40, xmax = 180,color='black')
# plt.xlim(left=40,right=180)
# plt.xlabel(r'3-Body Angle ($\theta$)',fontsize=15)
# plt.ylabel(r'$P(\theta) - P_{bulk}(\theta)$',fontsize=15)
# plt.legend(fontsize=8,ncol=2)
# leg = ileg()

# #%% INTEGRATE FROM 90-120 DEGREES TO FIND % OF TETRAHEDRAL WATER
# tetrah = []
# starti = np.argmin(np.abs(bins-100))
# endi   = np.argmin(np.abs(bins-120))

# # waters around each residue
# tetrahedral_ys = distros[:,starti:endi]
# frac_tetrahedral = np.sum(tetrahedral_ys,axis=1) * (bins[1]-bins[0]) # rectangular integration
# frac_tetrahedral = np.around(frac_tetrahedral,3)

# # bulk waters
# b_tetrahedral_ys = bulk_distro[starti:endi]
# b_frac_tetrahedral = np.sum(b_tetrahedral_ys) * (bins[1]-bins[0]) # rectangular integration
# b_frac_tetrahedral = np.around(b_frac_tetrahedral,3)

# rel_tetrahedrality = frac_tetrahedral / b_frac_tetrahedral


# #%% Plot relative tetrahedrality of sidechains compared to bulk water by residue
# plt.figure(figsize=(14,6))
# residues = [f'{seq[i]}{resid}' for i,resid in enumerate(np.arange(295,313+1))]
# # residues = np.delete(residues,[7,8,9])
# plt.bar(np.delete(residues,[7,8,9]),rel_tetrahedrality[:sidechain_counter])
# plt.ylim(bottom=0.8,top=1.2)
# plt.ylabel('Relative tetrahedrality',fontsize=15)
# plt.xticks(fontsize=12)
# plt.title('Relative tetrahedrality of water around SIDECHAIN heavy atoms',fontsize=20)

# #%% Plot relative tetrahedrality of backbone compared to bulk water by residue
# plt.figure(figsize=(14,6))
# plt.bar(residues,rel_tetrahedrality[sidechain_counter:])
# plt.ylim(bottom=0.8,top=1.2)
# plt.ylabel('Relative tetrahedrality',fontsize=15)
# plt.title('Relative tetrahedrality of water around BACKBONE heavy atoms',fontsize=20)
# plt.xticks(fontsize=12)
# plt.show()

# #%% Print Pymol commands to color heavy atoms by relative tetrahedrality

# # Scale rel_tetrahedrality values (ranging from 0.93 to 1.14) to a color value between 0.10 to 0.99
# # I'm setting the bluest blue to 1.10 rel_tetrahedrality and the reddest red to 0.90 rel_tetrahedrality

# low_val  = 0.90
# high_val = 1.10
# color_val = ( rel_tetrahedrality - low_val ) / (high_val - low_val)

# pymol_object = 'clust5_cutoff30' #name of object in Pymol
# sc_resids = np.concatenate((np.arange(1,8),np.arange(11,20))) # sidechains; ignoring glycines (resid 8,9,10)
# bb_resids = np.arange(1,19+1) #backbone

# # alter the β property to be the relative tetrahedrality around each residue's sidechain
# for i,resid in enumerate(sc_resids):
#     print(f'alter resid {resid} and sidechain and not hydrogen, b={rel_tetrahedrality[i]:.3f}')

# # alter the β property to be the relative tetrahedrality around each residue's backbone 
# for i, resid in enumerate(bb_resids):
#     print(f'alter resid {resid} and not sidechain and not hydrogen, b={rel_tetrahedrality[i+sidechain_counter]:.3f}')
    
#     # adding some if states to fix a Pymol bug (it counts the termini to be sidechains for some reason)
#     if resid==1:
#         print(f'alter resid {resid} and name N, b={rel_tetrahedrality[i+sidechain_counter]:.3f}')
#     if resid==19:
#         print(f'alter resid {resid} and (name OC1 or name OC2), b={rel_tetrahedrality[i+sidechain_counter]:.3f}')

# # color the residues on a red_white_blue spectrum
# print(f'spectrum b, red_white_blue, {pymol_object}, minimum={low_val}, maximum={high_val}')


# #%% Output csv of all the relative tetrahedrality bars
# # columns: hp1/hp2, noSalt/NaCl, cluster#, amino acid name, backbone/sidechain, rel tetrahedrality
# molecule = ['hp1'] * len(rel_tetrahedrality)
# salt_condition = ['noSalt'] * len(rel_tetrahedrality)
# cluster = [1] * len(rel_tetrahedrality)
# residues[6] = 'aa301'
# residues = list(np.delete(residues,[7,8,9])) + residues
# atomsType = ['sc']*sidechain_counter + ['bb']*(len(rel_tetrahedrality)-sidechain_counter)
# rel_tetrahedrality = np.around(rel_tetrahedrality,4)

# df = pd.DataFrame(list(zip(molecule,salt_condition,cluster,residues,atomsType,rel_tetrahedrality)),
#                columns =['molecule','salt condition','cluster#','amino acid',
#                          'backbone/sidechain','relative tetrahedrality'])

# df.to_csv('hp1_noSalt_cluster1.csv',index=False)
# %%
