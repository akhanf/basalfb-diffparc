import sklearn
import numpy as np
import nibabel as nib

# Define a function for saving niftis 
def save_label_nii (labels,mask,affine,out_nifti):
    labels_vol = np.zeros(mask.shape)
    labels_vol[mask > 0] = labels+1 #add a 1 so label 0 is diff from bgnd
    labels_nib = nib.Nifti1Image(labels_vol,affine)
    nib.save(labels_nib,out_nifti)

group_data = np.load(snakemake.input.connmap_group_npz)
k = int(snakemake.params.k)
out_nii = snakemake.output.cluster_k_indiv

conn_group = group_data['conn_group']
mask = group_data['mask']
affine = group_data['affine']

# Concat subjects
conn_group_m = np.moveaxis(conn_group,0,2)
conn_concat = conn_group_m.reshape([conn_group_m.shape[0],conn_group_m.shape[1]*conn_group_m.shape[2]])

# get indiv data connectivity
indiv_data = np.load(snakemake.input.connmap_indiv_npz)
conn_indiv = indiv_data['conn']

# Run spectral clustering and save output nifti
from sklearn.cluster import SpectralClustering

#create clustering object with parameters
clustering = SpectralClustering(n_clusters=k, assign_labels="discretize",random_state=0,affinity='cosine').fit(conn_concat)

#get centroids of each cluster 
centroids = np.zeros([k,conn_indiv.shape[1]])
print('shape of conn_concat: {}'.format(conn_concat.shape))
print('shape of conn_indiv: {}'.format(conn_indiv.shape))

print('shape of centroids: {}'.format(centroids.shape))
print('shape of clustering.labels_: {}'.format(clustering.labels_.shape))
for label in range(k):    
    #mean over all voxels with the label to get the centroid
    concat_centroid = np.mean(conn_concat[clustering.labels_==label,:],0)
    print('shape of concat_centroid: {}'.format(concat_centroid.shape))
    #since clustering was done on concatenated features, we want to average those to match indiv number of features
    # to do this, first reshape (unconcatenate):
    unconcat_centroid = concat_centroid.reshape([conn_group_m.shape[1],conn_group_m.shape[2]])
    print('shape of unconcat_centroid: {}'.format(unconcat_centroid.shape))
    # then take the mean:
    centroids[label,:] = np.mean(unconcat_centroid,1).T
    
import matplotlib.pyplot as plt

plt.figure()
plt.plot(centroids.T)
plt.savefig(snakemake.output.centroid_plot)

#apply it to individual
#compute cosine dis between conn_indiv and centroids
from sklearn.metrics.pairwise import cosine_similarity

indiv_sim = cosine_similarity(conn_indiv,centroids)
print('shape of indiv_sim: {}'.format(indiv_sim.shape))

#get index of maximal cosine similarity
indiv_labels = np.argmax(indiv_sim,1)

print('shape of indiv_labels: {}'.format(indiv_labels.shape))

save_label_nii(indiv_labels,mask,affine,out_nii)



