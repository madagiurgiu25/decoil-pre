"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:30 PM 9/19/22

Utils and auxiliary methods
"""

SEPARATOR = "_"

class PROG:

	# decoil modules - for advance use
	VALIDATE = "validate"
	RECONSTRUCT = "reconstruct"
	FILTER = "filter"
	CHECK = "check"

	# pipeline
	RECONSTRUCT_ONLY = "reconstruct-only"
	RECONSTRUCT_PLOT = "reconstruct-plot"
	# PLOT_ONLY = "plot-only"
	SV_ONLY = "sv-only"
	SV_RECONSTRUCT = "sv-reconstruct"
	# SV_RECONSTRUCT_PLOT = "sv-reconstruct-plot"

	DECOIL = "Decoil"

class META_CONFORMATION:

	N_FRAG = 0
	LEN_STR = 1
	SMALL_DEL = 2
	DUP = 3
	INV = 4
	INTERCHR = 5
	MULTI_REGION = 6
	FOLDBACK = 7
	RETURN = 8

	N_FRAG_STR = "#frag"
	LEN_STR = "len"
	SMALL_DEL_STR = "small_del"
	DUP_STR = "dup"
	INV_STR = "inv"
	INTERCHR_STR = "interch"
	MULTI_REGION_STR = "multi+region"
	FOLDBACK_STR = "foldback"
	RETURN_STR = "return"


class TOPOLOGY:

	SIMPLE_EXCISION = 0
	SIMPLE_EVENTS = 1
	MIXED_SIMPLE_EVENTS = 2
	MULTI_REGION_INTRA_CHR = 3
	MULTI_REGION_INTER_CHR = 4
	SIMPLE_DUPLICATIONS = 5
	FOLDBACKS = 6

	SIMPLE_EXCISION_STR = "simple_circle"
	SIMPLE_EVENTS_STR = "simple_sv"
	MIXED_SIMPLE_EVENTS_STR = "mixed_sv"
	MULTI_REGION_INTRA_CHR_STR = "multi_region_intra_chr"
	MULTI_REGION_INTER_CHR_STR = "multi_region_inter_chr"
	SIMPLE_DUPLICATIONS_STR = "simple_duplications"
	FOLDBACKS_STR = "foldbacks"

	DICT = {SIMPLE_EXCISION:SIMPLE_EXCISION_STR,
	        SIMPLE_EVENTS:SIMPLE_EVENTS_STR,
	        MIXED_SIMPLE_EVENTS:MIXED_SIMPLE_EVENTS_STR,
	        MULTI_REGION_INTRA_CHR:MULTI_REGION_INTRA_CHR_STR,
	        MULTI_REGION_INTER_CHR:MULTI_REGION_INTER_CHR_STR,
	        SIMPLE_DUPLICATIONS:SIMPLE_DUPLICATIONS_STR,
	        FOLDBACKS:FOLDBACKS_STR
	        }



class QUAL:
	# MIN_COV = 10
	MIN_COV = 8
	# MIN_VAF = 0.3
	MIN_VAF = 0.01
	# MIN_COV_ALT = 8
	MIN_COV_ALT = 6
	MINIMAL_SV_LEN = 500
	DISTANCE = 50
	# Set threshold for keeping fragments in the graph
	# MINIMAL_FRAGMENT_COVERAGE = 10
	MINIMAL_FRAGMENT_COVERAGE = 5
	# Set maximal distance for 2 fragments to be considered overlapped
	OVERLAP = 250
	# Set WGS mean cov
	MEAN_COVERAGE_WGS = 4
	# Set max coverage
	MAX_COVERAGE = 500
	# Set max coverage default
	MAX_COVERAGE_DEFAULT = 80000
	# Treshold ecDNA size (MB) Deshpande et al 2019. Nat.Comm.
	ECDNA_MINSIZE = 0.1
	# Minimal fragment size (bp)
	MINIMAL_FRAGMENT_SIZE = 500
	# Minimal score (estimated copy-number)
	FILTER_SCORE = 0
	# Filter variants based on a explog function
	EXPLOG_THRESHOLD = 0.1
	# Far
	FAR_FRAGMENTS = 16000


# MIN_COV = 10
# MIN_VAF = 0.2
# MIN_COV_ALT = 4
# DISTANCE = 50
# COVERAGE = 10


class POS:
	MIN = 0
	MAX = 100000000000


class GRAPH_PROP:
	HEAD = 1
	TAIL = 0
	
	FRAGMENT = 0
	DUP = 1
	DEL = 2
	INV = 3
	BND = 4
	INS = 5
	SPATIAL = 6
	REWIRE = 7
	SVTYPE = [FRAGMENT, DEL, DUP, INV, INS, BND, SPATIAL, REWIRE]


class GRAPH_NODE_PROP:
	# id of the fragment for which you create the tail and head nodes
	FRAGMENT_NAME = "parent_fragment"
	VISITED = "visited"
	FRAGMENT = "fragment"
	HEAD = "head"
	TAIL = "tail"
	COLOR = "color"
	CHR = 'CHR'
	POS = 'POS'
	START = 'START'
	END = 'END'


class GRAPH_EDGE_PROP:
	SVTYPE = "svtype"  # {fragment, sv, spatial}
	EDGEID = "edgeid"
	FRAGMENT = "fragment"
	SV = "sv"
	SPATIAL = "spatial"
	
	REPEAT = "repeat"
	
	# FRAGMENT_NAME = "frag"
	FRAGMENT_NAME = "parent_fragment"
	FRAGMENT_LEN = "fraglen"  # fragment len; applies only for SVTYPE == fragment
	LEN = "len"  # distance between breakpoints; applies only for SVTYPE == sv
	WEIGHT = "weight"
	COLOR = "color"
	VISITED = "visisted"
	
	CHR = 'CHR'
	POS = 'POS'
	START = 'START'
	END = 'END'


class TREE_PROP:
	SV_STATE = "sv"
	FRAGMENT_STATE = "fragment"
	DFS_CLASSIC = "classic"
	DFS_FRAGMENT = "fragment-aware"


class VCF_PROP:
	CHR1 = "CHR1"
	CHR2 = "CHR2"
	END = "END"
	START = "START"
	POS = "POS"
	
	DP = "DP"  # total supporting reads at locus
	DV = "DV"  # supporting reads for SV
	DR = "DR"  # supporting reads for reference
	GT = "GT"  # genotype
	GQ = "GQ"
	AF = "AF"
	
	DEL = "DEL"
	DUP = "DUP"
	INS = "INS"
	BND = "BND"
	TRA = "TRA"
	INV = "INV"
	INVDUP = "INVDUP"
	DUPINV = "DUPINV"
	DUPINS = "DUP/INS"
	
	SEQ = "SEQ"
	STRANDS = "STRANDS"
	STRAND = "STRAND"
	SVTYPE = "SVTYPE"
	SVLEN = "SVLEN"
	COVERAGE = "COVERAGE"
	SUPPORT = "SUPPORT"
	
	PRECISE = "PRECISE"
	IMPRECISE = "IMPRECISE"
	SNIFFLES2 = "sniffles2"
	SNIFFLES1 = "sniffles1"  # SV caller name
	CUTESV = "cutesv"
	NANOMONSV = "nanomonsv"
	
	SVCALLERS = [SNIFFLES1, SNIFFLES2, CUTESV, NANOMONSV]

	PASS = "PASS"
	STRANDBIAS = "STRANDBIAS"
	
	SV_COLORS = {DUP: "green",
	             DEL: "red",
	             INV: "blue",
	             BND: "black",
	             INS: "yellow",
	             INVDUP: "green",
	             DUPINS: "green",
				 TRA: "black"}
	
	# allowed SV';s
	SV_COLLECTION = [DEL, DUP, INV, INVDUP, BND, TRA, INS, DUPINS]
	SV_COLLECTION_STRICT = [DEL, DUP, INV, INVDUP, DUPINV, BND, TRA]
	
	ALLOWED_CHR = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
	               '19', '20', '21', '22', 'X', 'Y', 'M', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
	               'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
	               'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM',
	               'CP068255.2', 'CP068256.2', 'CP068257.2', 'CP068258.2', 'CP068259.2',
	               'CP068260.2', 'CP068261.2', 'CP068262.2', 'CP068263.2', 'CP068264.2', 'CP068265.2',
	               'CP068266.2', 'CP068267.2', 'CP068268.2', 'CP068269.2', 'CP068270.2', 'CP068271.2', 'CP068272.2',
	               'CP068273.2', 'CP068274.2', 'CP068275.2', 'CP068276.2', 'CP068277.2', 'CP086569.2', 'CP068254.1']


class CHROMOANASYNTHESIS:

	R_START = 1
	R_END = 2
	R_MOTIF = 3
	R_MOTIF_BINS = 5


class FRAGMENT:
	ID = 0
	STRAND = 1
	NAME = 2
	COV = 3
	
	FRAGID = "id"
	NODE_HEAD = "head_node_id"
	NODE_TAIL = "tail_node_id"
	EDGEID = "edge_id"
	
	HEAD = 1
	TAIL = 0
	
	MB = 1000000
	COVERAGE = "coverage"
	NORM_COV = "norm_cov"
	ORIENTATION = "strand"
	LENGTH = "len"
	CHR = "chr"
	START = "start"
	END = "end"


class CLUSTER:
	CIRCLES = "circles"
	FRAGMENTS = "fragments"
	EDGES = "edges"


class PATH:
	CURRENT = 0
	FROM = 1
	EDGE_STATE = 2
	EDGE_ID = 3
	
	CHR = 0
	START = 1
	END = 2
	STRAND = 3
	FRAG = 4
	CIRCID = 5
	COV = 6
	SCORE = 7

	PATH_FRAGMENTS = "path_fragments"
	PATH_EDGES = "path_edges"
	CIRCID_NAME = "circ_id"


class BED_PROP:
	CHR = "#chr"
	START = "start"
	END = "end"
	CIRC_ID = "circ_id"
	FRAGMENT_ID = "fragment_id"
	STRAND = "strand"
	COVERAGE = "coverage"
	POSITIVE = "+"
	NEGATIVE = "-"
	LEFT = "0"
	RIGHT = "1"
	STRAND_IDX = 5


class LINKS_PROP:
	N1 = "n1"
	N2 = "n2"
	N1SIDE = "n1.side"
	N2SIDE = "n2.side"
	CIRC_ID = "circ_id"
	FRAGMENT_ID = "fragment_id"


