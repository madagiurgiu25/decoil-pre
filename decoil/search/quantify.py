"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 4:04 PM 4/18/22

Quantify the likely circles (ecDNA) in the DFS tree
"""
import logging
import sys
import warnings
from collections import defaultdict

from copy import deepcopy
import numpy as np
import pandas as pd
import scipy
from numpy import dot
from numpy.linalg import norm

from sklearn.linear_model import Ridge, Lasso, SGDRegressor, LinearRegression, ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler

from decoil.utils import CLUSTER as cl
from decoil.utils import PATH as pt
from decoil.validate import compare as compare
from decoil.utils import QUAL

warnings.filterwarnings("ignore")

log = logging.getLogger('reconstruct.clean')
log.propagate = False
log.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)


def cosvectors(a, b):
	return dot(a, b) / (norm(a) * norm(b))


def compute_similarity(a, b):
	prod = a * b
	return cosvectors(prod, a), cosvectors(prod, b)


def print_similarity(df):
	"""
	Return similarity of the vectors in the matrix
	"""
	M = np.array(df)
	
	if M.shape[1] > 2:
		for i in range(0, M.shape[1] - 1):
			for j in range(i + 1, M.shape[1]):
				print(i, j, compute_similarity(M[:, i], M[:, j]))
	else:
		print(0, 1, compute_similarity(M[:, 0], M[:, 1]))


def error(u, v):
	"""
	Compute error between 2 vectors.
	Assume u - predicted, v - true
	"""
	diff = np.subtract(u, v)
	return scipy.linalg.norm(diff) / scipy.linalg.norm(v)


def mse(u, v):
	"""
	Compute mean square error between 2 vectors
	Assume u - predicted, v - true
	"""
	return mean_squared_error(v, u)  # true, predicted


def abs_total_error(u, v):
	"""
	Compute absolute max error
	Assume u - predicted, v - true
	"""
	return np.sum(np.abs(np.subtract(u, v)))


def circles_to_matrix(all_cycles, all_fragments):
	"""
	Convert circle and fragments to matrix.
	Consider fragments fid = -fid same (do not distinguish between forward and inverted fragment
	"""

	dict_circles2fragments = defaultdict(list)
	dict_fragments2circles = defaultdict(list)

	# all canonical cycles
	dict_cycles = {}

	df = pd.DataFrame(0, columns=all_fragments, index=list(range(len(all_cycles))))

	for i ,cycle in enumerate(all_cycles):
		dict_cycles[i] = cycle
		# f is a fragment id (if negative the fragment is reversed)
		for f in cycle:
			# ignore strandness
			afid = abs(f)
			dict_fragments2circles[afid].append(i)
			dict_circles2fragments[i].append(afid)
			df.loc[i, afid] += 1

	# remove repetitive fragments in the structure
	for c in dict_circles2fragments:
		dict_circles2fragments[c] = list(set(dict_circles2fragments[c]))

	for f in dict_fragments2circles:
		dict_fragments2circles[f] = list(set(dict_fragments2circles[f]))

	return df, dict_circles2fragments, dict_fragments2circles, dict_cycles


def find_overlaps(dict_fragments2circles, dict_circles2fragments):
	"""
	Find the overlapping circles using the matrix (circles x fragments) = m x n
	"""
	fragment_keys = list(dict_fragments2circles.keys())

	# initialize with max size cluster (number of fragments)
	clusters = defaultdict(dict)
	for c in range(len(fragment_keys)):
		clusters[c]["circles"] = []
		clusters[c]["fragments"] = []

	count = 0
	for f in dict_fragments2circles:
		hit = 0
		# check if fragment can be already assigned to a cluster
		for i in range(0, count):
			found_overlap = False
			if f in clusters[i]:
				found_overlap = True
				hit += 1
			else:
				for c in dict_fragments2circles[f]:
					if c in clusters[i]["circles"]:
						found_overlap = True
						hit += 1
						break

			if found_overlap:
				clusters[i]["fragments"].append(f)
				for c in dict_fragments2circles[f]:
					clusters[i]["circles"].append(c)

		# no overlap found with clusters
		if hit == 0:
			clusters[count]["fragments"].append(f)
			for c in dict_fragments2circles[f]:
				clusters[count]["circles"].append(c)
			count += 1
	
	for c in range(len(clusters)):
		if len(clusters[c]["circles"]) == 0:
			clusters.pop(c)

	# remove redundancy
	for c in clusters:
		clusters[c]["circles"] = list(set(clusters[c]["circles"]))
		clusters[c]["fragments"] = list(set(clusters[c]["fragments"]))
	
	# merge clusters
	merge = True
	while merge:
		
		found_overlap = False
		ck = list(clusters.keys())
		for i in range(0, len(ck)):
			stop = False
			for j in range(i, len(ck)):
				if i != j and np.intersect1d(clusters[ck[i]]["fragments"], clusters[ck[j]]["fragments"]).size > 0:
					stop = True
					found_overlap = True
					break
			if stop == True:
				break
		
		if found_overlap:
			clusters[ck[i]]["circles"] = list(set(clusters[ck[i]]["circles"] + clusters[ck[j]]["circles"]))
			clusters[ck[i]]["fragments"] = list(set(clusters[ck[i]]["fragments"] + clusters[ck[j]]["fragments"]))
			clusters.pop(ck[j])
			merge = True
		else:
			merge = False
	
	return clusters


def compute_cov_per_fragment(df, coverage):
	p = np.array(df) * coverage
	cov_per_fragment = np.sum(p, axis=1)
	return cov_per_fragment


def compute_cutoff(fragments_cov):
	"""
	Compute the cutoff above which we consider likely circle.

	Arguments:
		fragments_cov (list): Coverage list of the fragments composing the circles
	"""
	# m = np.mean(np.array(fragments_cov))
	# std = np.std(np.array(fragments_cov))
	#
	# return (m - std) / 4
	# print("compute cutoff", fragments_cov, np.min(np.array(fragments_cov)))
	# return max(np.min(np.array(fragments_cov)) / 4, 10)
	# return (m - std) / 4
	return max((np.min(np.array(fragments_cov)) / 5), 10)


def weight_fragments(frag_list, dict_frags):
	"""
	Weight fragments using 1/length.

	Arguments:
		frag_list (list): Fragment IDs
		dict_frags (dict): Properties of fragments

	Returns:
		np.array
	"""
	fragment_size = np.array([dict_frags[abs(id)].len for id in frag_list])
	fragment_coverage = np.array([dict_frags[abs(id)].coverage for id in frag_list])
	total_size = np.sum(fragment_size)
	total_coverage = np.sum(fragment_coverage)
	wfragsize =  (- (fragment_size - total_size)) / total_size
	wcovsize  = fragment_coverage / total_coverage

	print(wfragsize)
	print(wcovsize)
	return wfragsize * wcovsize


def add_coverage_covariate_to_lasso(df, y, wgs=QUAL.MEAN_COVERAGE_WGS):
	"""
	Model coverage as feature in Lasso.
	This will allow to account for WGS coverage in the model

	(fragments x circles), where w - whole genome coverage
	   B   ABC  ABBC | ABC
	A  0    1    1   |  1
	B  1    1    2   |  1
	C  0    1    1   |  1
	---------------------
	w  0    0    0   |  1

	y = [A, B, C, D, w] - coverage for each fragments

	Arguments:
		df (pd.DataFrame) : Matrix (fragments x circles)

	Returns:
		Copy of df and y objects, with the coverage as covariate
	"""
	df_temp = deepcopy(df)
	y_temp = deepcopy(y)


	df_temp.loc["wgs"] = np.zeros(df_temp.shape[1])
	df_temp["wgs_circle"] = np.ones(df_temp.shape[0])
	y_temp.append(wgs)
	return df_temp, y_temp


def regress(df_original, y_original, wf, pseudocount=0.001, cutoff=4):
	"""
	Use Lasso to find likely circles that match the coverage profile.
	Lasso Input (fragments x circles) where circles are features.

	Arguments:
		df (pd.DataFrame): Matrix fragments x circles
		y (list): Total coverage per fragment
		wf (list): Fragments weights
		pseudocount (float): Remove 0 values in the matrix
		cutoff (float): Consider likely circles having a coef_ >= cutoff
	"""
	print("Perform lasso regression")

	if df_original.shape[0] > 1:

		# add coverage as covariate
		df, y = add_coverage_covariate_to_lasso(df_original, y_original, wgs=QUAL.MEAN_COVERAGE_WGS)
		print(df)
		print(y)

		# add pseudocounts
		M = np.array(df) + pseudocount

		# weight by 1/fragment length
		# M = M * wf
		# print(wf)
		# print(M)
		# print(y)
		# print(cutoff)
		# total_cov = y
		# regr = ElasticNet(random_state=0, alpha=1, l1_ratio=0.7).fit(M.T, y)
		# filtered = np.ravel(np.argwhere(regr.coef_ >= 2))
		# print(regr.coef_)
		# return df.iloc[filtered, :].index.tolist(), np.take(regr.coef_, filtered)

		# fit
		lasso = Lasso(alpha=0.1, max_iter=1000, fit_intercept=True).fit(M, y)
		# remove wgs coverage covariate from lasso coefs
		lasso_coefs = lasso.coef_[:-1]
		# print(lasso.coef_)
		print(lasso_coefs)
		filtered = np.ravel(np.argwhere(lasso_coefs >= cutoff))

		# something went wrong / multicollinearity
		if len(filtered) == 0:
			filtered = np.ravel(np.argwhere(lasso_coefs >= 0))

		print("after filtering")
		print(np.take(lasso_coefs, filtered))

		return df_original.iloc[:, filtered].columns.tolist(), np.take(lasso_coefs, filtered)

	else:
		# one circle (no need for fitting)
		print(df_original.index.tolist())
		return df_original.index.tolist(), [100]


# return abs_total_error(lasso.coef_, proportions)


def regress_tests(df, y, proportions, pseudocount=0.1, scaled=False):
	"""
	Testing multiple linear regressions
	df - Matrix fragments x circles
	y - total coverage per fragment
	proportions - true circle composition
	pseudocount - remove 0 values in the matrix
	"""
	
	M = np.array(df) + pseudocount
	if scaled == True:
		scaler = StandardScaler()
		scaler.fit(M)
		Mscaled = scaler.transform(M)
	else:
		Mscaled = M
	
	lasso = Lasso(alpha=0.1).fit(Mscaled, y)
	#     lassocv = LassoCV(cv=3, random_state=0).fit(Mscaled, y)
	ridge = Ridge().fit(Mscaled, y)
	linear = LinearRegression().fit(Mscaled, y)
	sgd = SGDRegressor(max_iter=1000, tol=1e-3, average=True).fit(Mscaled, y, coef_init=[100] * Mscaled.shape[1])
	
	return (abs_total_error(lasso.coef_, proportions),
	        abs_total_error(ridge.coef_, proportions),
	        abs_total_error(linear.coef_, proportions),
	        abs_total_error(sgd.coef_, proportions))


def mean_circle(path, fragments):
	"""
	Compute weighted mean of a circular path
	"""
	a = [fragments[item[pt.FRAG]].coverage for item in path]
	alen = [abs(item[pt.START] - item[pt.END]) for item in path]
	alensum = sum(alen)
	alen_w = np.array(alen) / alensum
	return sum(np.array(a) * alen_w)


def compute_proportions(matrix, clusters, fragments, overlapping_circles):
	"""
	Estimate the circles abundance for overlapping circles.

	Arguments:
		matrix (pd.DataFrame): circlex x fragments matrix
		clusters (dict): Contains the representation of all overlapping circles in one cluster
		fragments (decoil.graph.Fragment): Information about the fragments
		all_cycles (dict): All cycles
	"""
	likely_circles = defaultdict()
	
	# quantify proportions per cluster
	for c in clusters:

		# at least 2 overlapping circles
		if len(clusters[c][cl.CIRCLES]) > 1:
			rows = clusters[c][cl.CIRCLES]
			cols = sorted(clusters[c][cl.FRAGMENTS])
			submatrix = matrix.loc[rows, cols]
			# transpose in fragments x circles format (prepare Lasso input)
			submatrix = submatrix.T

			# collapse the total coverage per fragments used for the regression
			coverage_per_fragments = [fragments[id_frag].coverage for id_frag in cols]

			# cutoff above which we consider likely circles
			cutoff = compute_cutoff(coverage_per_fragments)
			wf = weight_fragments(cols, fragments)
			candidates, coef = regress(submatrix, coverage_per_fragments, wf, cutoff=cutoff)

			# add all likely circles
			for i, cand in enumerate(candidates):
				likely_circles[cand] = {}
				likely_circles[cand]["conf"] = overlapping_circles[cand]
				likely_circles[cand]["score"] = coef[i]
		else:
			# singletons
			cid = clusters[c][cl.CIRCLES][0]
			fids = clusters[c][cl.FRAGMENTS]
			likely_circles[cid] = {}
			likely_circles[cid]["conf"] = overlapping_circles[cid]
			likely_circles[cid]["score"] = np.mean([fragments[id_frag].coverage for id_frag in fids])
	
	return likely_circles


# circular_paths_fragments, circular_paths_coverage
# matrix = circles_to_matrix(circles_with_subcircles, coverage_with_subcircles)
# ytrue = compute_cov_per_fragment(matrix, coverage_with_subcircles)
# (elasso, eridge, elinear, esgd), (coef_lasso, coef_ridge, coef_linear, coef_sgd) = regress(matrix, ytrue, coverage_with_subcircles)
def run_quantification(graph):
	"""
	Compute proportions of these circles to find likely true ecDNA
	"""
	# 1. add all unique circles to likely candidates
	log.info("4. Compute likely circles")

	log.info("4.0 Filter out highly similar circles")
	all_cycles = list(graph.get_hash_canonical().values())
	all_fragments = graph.get_fragments()
	# filtered_cycles = compare.filter_out_similar_cycles(all_cycles,graph)

	# 2. convert circles to matrix
	log.info("4.1. Build identity matrix")

	matrix, dict_circles2fragments, dict_fragments2circles, dict_cycles = circles_to_matrix(all_cycles, all_fragments)

	# 3. find overlapping circles
	log.info("4.2 Find clusters")
	clusters = find_overlaps(dict_fragments2circles, dict_circles2fragments)

	# 4. quantify circles
	log.info("4.3. Quantify circles")
	likely_candidates_overlapping = compute_proportions(matrix, clusters, all_fragments, dict_cycles)
	return likely_candidates_overlapping, clusters
