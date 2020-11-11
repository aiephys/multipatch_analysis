import os, pickle, gc
import numpy as np
import umap
import sklearn.preprocessing, sklearn.decomposition
from aisynphys.stochastic_release_model import StochasticModelRunner, load_cached_model_results
from aisynphys.database import default_db as db
from aisynphys import config
from aisynphys.ui.progressbar import ProgressBar

base_path = os.getcwd()
cache_path = os.path.join(base_path, 'cache')
model_result_cache_path = os.path.join(cache_path, 'stochastic_model_results')

cache_files = [os.path.join(model_result_cache_path, f) for f in os.listdir(model_result_cache_path)]

## Load all model outputs into a single array
agg_result, cache_files, param_space = load_cached_model_results(cache_files[:30], db=db)
agg_shape = agg_result.shape
flat_result = agg_result.reshape(agg_shape[0], np.product(agg_shape[1:]))[:, :100]


# Prescale model data
scaler = sklearn.preprocessing.StandardScaler()
n_obs = flat_result.shape[0]
chunk_size = 10
n_chunks = n_obs // chunk_size

scaler.fit(flat_result)
scaled = scaler.transform(flat_result)

# free up some memory
del agg_result
del flat_result
gc.collect()

# fit sparse PCA
print("Fitting sparse PCA...")
n_pca_components = 50
pca = sklearn.decomposition.MiniBatchSparsePCA(n_components=n_pca_components)
pca.fit(scaled)
print("  Sparse PCA fit complete.")

# run sparse PCA
print("Sparse PCA transform...")
sparse_pca_result = pca.transform(scaled)
sparse_pca_file = os.path.join(cache_path, 'sparse_pca.pkl')
pickle.dump({'result': sparse_pca_result, 'params': param_space, 'cache_files': cache_files, 'sparse_pca': pca}, open(sparse_pca_file, 'wb'))
print("   Sparse PCA transform complete: %s" % sparse_pca_file)

# fit standard PCA
print("Fitting PCA...")
n_pca_components = 30#500
pca = sklearn.decomposition.PCA(n_components=n_pca_components)
pca.fit(scaled)
print("  PCA fit complete.")

# run PCA
print("PCA transform...")
pca_result = pca.transform(scaled)
pca_file = os.path.join(cache_path, 'pca.pkl')
pickle.dump({'result': pca_result, 'params': param_space, 'cache_files': cache_files, 'pca': pca}, open(pca_file, 'wb'))
print("   PCA transform complete: %s" % pca_file)

# umap
n_umap_components = 15#32
reducer = umap.UMAP(
    n_components=n_umap_components,
#     n_neighbors=5,
)

print("Fit UMAP...")
reducer.fit(scaled)
print("   UMAP fit done.")

print("UMAP transform...")
umap_result = reducer.transform(scaled)
umap_file = os.path.join(cache_path, 'umap.pkl')
pickle.dump({'result': umap_result, 'params': param_space, 'cache_files': cache_files, 'umap': reducer}, open(umap_file, 'wb'))
print("   UMAP transform complete: %s" % umap_file)
