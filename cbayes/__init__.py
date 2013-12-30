"""cmpy_bayes

A set of scripts to make use of Bayesian methods for inference of unifilar HMMs
and epsilon-machines CMPy easier.

"""
from __future__ import absolute_import

from .util_general import check_dir_exists
from .util_general import check_positive_float
from .util_general import check_positive_int
from .util_general import check_probability
from .util_general import check_sr_tuple
from .util_general import correct_neg_zero
from .util_general import create_dir
from .util_general import read_datafile
from .util_general import read_sample_dir

from .util_infer_db import add_topologies_to_db
from .util_infer_db import calc_probs_beta_db
from .util_infer_db import create_sample_summary_file
from .util_infer_db import prior_add_topologies_to_db
from .util_infer_db import sample_db