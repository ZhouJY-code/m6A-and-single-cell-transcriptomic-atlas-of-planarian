import os
import sys
import numpy as np
from scipy.stats import norm
from scipy.spatial.distance import pdist, squareform
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import loompy
import velocyto as vcy
import pandas as pd


#files = ["cut-ctrl-0d.loom","cut-ctrl-3d.loom"]
#loompy.combine(files, "merge.loom", key="Accession")
vlm = vcy.VelocytoLoom("merge.loom")

velo_list = [x[0:28] for x in list(vlm.ca["CellID"])]
#keep_cell = [cell in annotation.index.tolist() for cell in vlm.ca['CellID']]
annotation = pd.read_csv("cluster_iden_umap.csv", delimiter=',',index_col=0)

keep_cell = [cell in annotation.index.tolist() for cell in velo_list]
vlm.filter_cells(keep_cell)

### Get seurat cells and annotation
#tsneCoord = pd.read_csv(base_dir+"tsneCoordinates.csv", delimiter=',', index_col=0)
#annotation = pd.merge(annotation, tsneCoord, left_index=True, right_index=True)
#annotation = annotation.set_index('type')

annotation = annotation.loc[velo_list]

#color_dict = dict(zip(list(annotation["type"].value_counts().index), list(annotation["type"].value_counts().index)))

# add cluster, color and time as annotation from Seurat object to velocyto object
vlm.set_clusters(cluster_labels=list(np.array(annotation["type"])))
vlm.ca["Clusters"] = vlm.cluster_ix
vlm.ca["time"] = np.array(annotation["orig.ident"])
vlm.ca["type"] = np.array(annotation["type"])
vlm.ts=annotation.loc[:,["UMAP_1","UMAP_2"]].values

plt.figure(figsize=(15,15))
vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:,1], c=vlm.colorandum, s=5)
for i in range(max(vlm.ca["Clusters"])+1):
    ts_m = np.median(vlm.ts[vlm.ca["Clusters"] == i, :], 0)
    plt.text(ts_m[0], ts_m[1], str(vlm.cluster_labels[vlm.ca["Clusters"] == i][0]),
             fontsize=13, bbox={"facecolor":"w", "alpha":0.6})
plt.axis("off")
plt.savefig("umap.png")

vlm.plot_fractions("fractions.png")
plt.savefig("fractions.png",dpi=300)

vlm.normalize("S", size=True,  log=False)
vlm.normalize("U", size=True,  log=False)

plt.figure(figsize=(15,15))
vlm.score_cv_vs_mean(4000, plot=True, max_expr_avg=35)
print(sum(vlm.cv_mean_selected))
vlm.filter_genes(by_cv_vs_mean=True)
plt.savefig("filter_genes.png")

vlm.score_detection_levels(min_expr_counts=3, min_cells_express=3, min_expr_counts_U=3, min_cells_express_U=3)
vlm.filter_genes(by_detection_levels=True)
print("Number of genes to be used:",vlm.S.shape[0])

vlm.normalize_by_total()

plt.figure(figsize=(15,15))
vlm.perform_PCA()
plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
print("number of PCs to be used:",n_comps)
plt.savefig("pca.png")

k = 200
vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=15)

vlm.normalize_median()
vlm.fit_gammas(limit_gamma=True)

vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1)

vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt",n_neighbors=1000, knn_random=True, sampled_fraction=1,threads=15)
#vlm.to_hdf5("estimate_transition_prob.hdf5",pickle_protocol=4)
#vcy.to_hdf5(vlm,"estimate_transition_prob.hdf5")
#vlm=[vlm.encode('utf8') for header in vlm]
vlm.calculate_embedding_shift()
#vlm.to_hdf5("calculate_embedding_shift")


#edgecolorNone,2.reshape(1, -1)

quiver_scale = 110
ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
quiver_kwargs=dict(headaxislength=6, headlength=20, headwidth=7,linewidths=0.2, width=0.0045,edgecolors="k", color=vlm.colorandum[ix_choice], alpha=1)

plt.figure(None,(20,20))
plt.scatter(vlm.embedding[:, 0].reshape(1, -1), vlm.embedding[:, 1].reshape(1, -1),c="0.8", alpha=0.2, s=10, edgecolor=None)
plt.scatter(vlm.embedding[ix_choice, 0].reshape(1, -1), vlm.embedding[ix_choice, 1].reshape(1, -1),
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.1)
plt.quiver(vlm.embedding[ix_choice, 0].reshape(1, -1), vlm.embedding[ix_choice, 1].reshape(1, -1),
           vlm.delta_embedding[ix_choice, 0].reshape(1, -1), vlm.delta_embedding[ix_choice, 1].reshape(1, -1),
           **quiver_kwargs)
plt.axis("off")
plt.savefig("velocity_embedding.png")


plt.figure(figsize=(15,15))
vlm.calculate_grid_arrows(smooth=0.5, steps=(50, 50), n_neighbors=300)
vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":50, "rasterized":True}, min_mass=15, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=1, headwidth=4.8, quiver_scale=0.5)
plt.savefig("velocity_field.png")



#Markov

steps = 200, 200
grs = []
for dim_i in range(vlm.embedding.shape[1]):
    m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
    m = m - 0.025 * np.abs(M - m)
    M = M + 0.025 * np.abs(M - m)
    gr = np.linspace(m, M, steps[dim_i])
    grs.append(gr)

meshes_tuple = np.meshgrid(*grs)
gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

diag_step_dist = np.sqrt((meshes_tuple[0][0,0] - meshes_tuple[0][0,1])**2 + (meshes_tuple[1][0,0] - meshes_tuple[1][1,0])**2)
min_dist = diag_step_dist / 2
ixs = ixs[dist < min_dist]
gridpoints_coordinates = gridpoints_coordinates[dist.flat[:]<min_dist,:]
dist = dist[dist < min_dist]
ixs = np.unique(ixs)

plt.figure(figsize=(15,15))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=vlm.colorandum[ixs], alpha=1, s=30, lw=0.4,
                edgecolor="0.4")

plt.savefig("velocity_Markov_dot.png")


vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='forward', cells_ixs=ixs)
vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=2500)

diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 97)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(figsize=(15,15))
vcy.analysis.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.8, s=10, lw=0.,
                edgecolor=None, cmap="viridis_r")

cm = plt.cm.get_cmap('viridis_r')
plt.colorbar()
plt.axis("off")
plt.savefig("velocity_Markov_end.png")



vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='backwards', cells_ixs=ixs)
vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=3000)

diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 98)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(figsize=(15,15))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.8, s=10, lw=0.,
                edgecolor=None, cmap="viridis_r")
cm = plt.cm.get_cmap('viridis_r')
plt.colorbar()
plt.axis("off")
plt.savefig("velocity_Markov_new.png")

