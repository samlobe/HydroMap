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

