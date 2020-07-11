import sklearn
import numpy as np
import nibabel as nib
import errno
import os
from os.path import join
from sklearn.cluster import SpectralClustering
from sklearn.metrics.pairwise import cosine_similarity
from scipy.io import savemat

# Define a function for saving niftis 
def save_label_nii (labels,mask,affine,out_nifti):
    labels_vol = np.zeros(mask.shape)
    labels_vol[mask > 0] = labels+1 #add a 1 so label 0 is diff from bgnd
    labels_nib = nib.Nifti1Image(labels_vol,affine)
    nib.save(labels_nib,out_nifti)

group_data = np.load(snakemake.input.connmap_group_npz)
out_file_prefix = snakemake.params.out_file_prefix
cluster_indiv_dir = snakemake.output.cluster_indiv_dir
cort_profiles_dir = snakemake.output.cort_profiles_dir
max_k = int(snakemake.params.max_k)


try:
    os.mkdir(cluster_indiv_dir)
    os.mkdir(cort_profiles_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

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


for k in range(2,max_k):

    #create clustering object with parameters
    clustering = SpectralClustering(n_clusters=k, assign_labels="discretize",random_state=0,affinity='cosine').fit(conn_concat)

    #get centroids of each cluster 
    centroids = np.zeros([k,conn_indiv.shape[1]])
    #print('shape of conn_concat: {}'.format(conn_concat.shape))
    #print('shape of conn_indiv: {}'.format(conn_indiv.shape))

    #print('shape of centroids: {}'.format(centroids.shape))
    #print('shape of clustering.labels_: {}'.format(clustering.labels_.shape))
    for label in range(k):    
        #mean over all voxels with the label to get the centroid
        concat_centroid = np.mean(conn_concat[clustering.labels_==label,:],0)
        #print('shape of concat_centroid: {}'.format(concat_centroid.shape))
        #since clustering was done on concatenated features, we want to average those to match indiv number of features
        # to do this, first reshape (unconcatenate):
        unconcat_centroid = concat_centroid.reshape([conn_group_m.shape[1],conn_group_m.shape[2]])
        #print('shape of unconcat_centroid: {}'.format(unconcat_centroid.shape))
        # then take the mean:
        centroids[label,:] = np.mean(unconcat_centroid,1).T
    

    #apply it to individual
    #compute cosine dis between conn_indiv and centroids

    indiv_sim = cosine_similarity(conn_indiv,centroids)
    #print('shape of indiv_sim: {}'.format(indiv_sim.shape))

    #get index of maximal cosine similarity
    indiv_labels = np.argmax(indiv_sim,1)

    #print('shape of indiv_labels: {}'.format(indiv_labels.shape))

    out_nii = join(cluster_indiv_dir,out_file_prefix+'_k-{k}_cluslabels.nii.gz'.format(k=k))

    save_label_nii(indiv_labels,mask,affine,out_nii)

    #for cortical profiles, generate the average profile for each label:
    cort_profiles = np.zeros([k,conn_indiv.shape[1]])

    for label in range(k):
        cort_profiles[label,:] = np.mean( conn_indiv[indiv_labels ==  label,:],0)

    out_cort_profiles_npz = join(cort_profiles_dir,out_file_prefix+'_k-{k}_cortprofiles.npz'.format(k=k))
    out_cort_profiles_mat = join(cort_profiles_dir,out_file_prefix+'_k-{k}_cortprofiles.mat'.format(k=k))
    np.savez(out_cort_profiles_npz,cort_profiles=cort_profiles)

    savemat(out_cort_profiles_mat,{'cort_profiles': cort_profiles})

