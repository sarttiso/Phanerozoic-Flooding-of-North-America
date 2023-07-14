import numpy as np
import pandas as pd
import geopandas as gpd
import scipy.stats as stats
from tqdm.notebook import tqdm

"""
LOAD AND PREPROCESS DATA

The following functions load and preprocess Macrostrat data organized by columns and units.
"""
# choose North American Lambert projection
epsg = 2252
proj4 = 'epsg:%d' % epsg

# IN:
# shp: string for file location (e.g., 'data/units.shp')
# index: string for column to use as index (e.g., 'unit_id')
# project_id (=1): default project id is 1 for North America. Set to none to take all columns.
def load_shapefile(shp, index, project_id=1):
    df = gpd.read_file(shp)
    df.set_index(index, inplace=True, drop=False)
    df = df.to_crs(proj4)
    if not project_id==None:
        df = df[df['project_id']==project_id]
    return df

"""
From a series of environmental strings delimited by | (as returned by macrostrat api), get the corresponding environment ids
- currently assumes that all strings in env_str are valid
IN:
env_str: array of the |-separated environment strings given in the units macrostrat API
env: reference dataframe with a key column matching all of the unique environment strings
"""
def env_id(env_str, env):
    n = len(env_str)
    env_ids = []
    for ii in tqdm(range(n), desc='Generating macrostrat environment IDs'):
        # check that string exists, otherwise put nan
        if (env_str[ii] == None) or not isinstance(env_str[ii], str):
            env_ids.append(np.nan)
            continue
        cur_envs = env_str[ii].split('|')
        env_id = []
        for cur_env in cur_envs:
            cur_env = cur_env.rstrip()
            # find matching string in env key (logical indices)
            env_id_idx = (env['key'].str.find(cur_env) != -1).values
            # if it is not present, put nan
            if np.all(env_id_idx == 0):
                env_id.append(np.nan)
                continue
            # convert to linear index of environment in the env dataframe
            env_id_idx = np.argwhere(env_id_idx)[0][0]
            # get the environment id corresponding to this index
            env_id.append(env.iloc[env_id_idx]['environ_id'])
        env_ids.append(env_id)
    return env_ids


"""
Returns array of logical indices into array of environmental indices (i.e. output from env_id) for units that match the query criteria, which allow the user to use the hierarchy of environmental labeling embedded in the Macrostrat data. The user can either specify class/class-type combinations. Choosing a specific name is unnecessary because each name is already uniquely identified by an environmental id.
IN:
env_ids: list of arrays of environment ids as output by env_id()
class_: (default 'marine') either 'marine' or 'non-marine'
type_: type, need not be specified
each: (default False) whether or not each environment, for units that have multiple environments, must satisfy the query criteria, or if only one is sufficient. If true, then all environments must satisfy the class-type combination, if false then only one must.
OUT:

"""
def env_get(env_ids, env, class_='marine', type_=None, each=False):
    assert class_ in ['marine', 'non-marine'], 'class must be marine or non-marine'
    if not type_==None:
        assert type_ in env['type'].unique(), 'typ is not valid string'
    # number of units
    n = len(env_ids)

    # if type not specified, just make it all ones, otherwise find matches in env table and then matches in env_ids
    if type_==None:
        type_idx = np.ones(n).astype(np.bool)
    else:
        type_ids = env.loc[env['type']==type_]['environ_id'].values
        type_idx = np.zeros(n).astype(np.bool)
        if each:
            for ii in range(n):
                type_idx[ii] = np.all(np.in1d(env_ids[ii], type_ids))
        else:
            for ii in range(n):
                type_idx[ii] = np.any(np.in1d(env_ids[ii], type_ids))

    # now do the same for the class
    class_ids = env.loc[env['class']==class_]['environ_id'].values
    class_idx = np.zeros(n).astype(np.bool)
    if each:
        for ii in range(n):
            class_idx[ii] = np.all(np.in1d(env_ids[ii], class_ids))
    else:
        for ii in range(n):
            class_idx[ii] = np.any(np.in1d(env_ids[ii], class_ids))

    # get overlap of resulting indices, which forms the final output array of logical indices
    idx = class_idx & type_idx
    return idx


"""
MACROSTRAT THICKNESS ESTIMATES PROCESSING

The following functions assist in the computation of minimum and maximum sedimentary thickness estimates within columns and specified time intervals.
"""

"""
compute overlap weights for each unit and the given time bin
"""
def bin_overlap_correction(units, t1, t2):
    d = units['b_age'].values - units['t_age'].values
    a = np.vstack((np.zeros((len(units))),
                   units['b_age'].values - t1)).T
    a = np.max(a, axis=1)
    b = np.vstack((np.zeros((len(units))),
                   t2 - units['t_age'].values)).T
    b = np.max(b, axis=1)
    # "bug" here in that I don't account for possibly negative values that would result from (incorrect) Macrostrat database values in which b_age is younger than t_age. needs to be fixed (or just correct database before using this function)
    w = (d - a - b)/d
    w[w < 0] = 0
    return w

"""
idx_groups: dictionary with col_ids and indices of units corresponding to each column
units: dataframe indexed by idx_groups into columns

generate thickness corrections to subtract from the thicknesses naively computed in seds_by_column_and_age
correction are then positive
"""
def unit_overlap_correction(idx_groups, units_in):
    cols = idx_groups.keys()
    # create dataframe of corrections (default zero correction)
    corr = pd.DataFrame(index=cols, columns=['max_corr', 'min_corr'])
    corr[:] = 0
    # loop over columns (might be faster with smart groupby() call)
    for col in cols:
        # units in column
        col_units = units_in.iloc[idx_groups.get(col)][['t_age', 'b_age', 'max_thick', 'min_thick']]
        # loop over units
        n_units = len(col_units)
        ovlap_coords = []
        ovlap_units = []
        # find overlapping units
        for ii in range(n_units):
            for jj in range(ii+1, n_units):
                ti_t, ti_b = col_units.iloc[ii]['t_age'], col_units.iloc[ii]['b_age']
                tj_t, tj_b = col_units.iloc[jj]['t_age'], col_units.iloc[jj]['b_age']
                if (ti_b >= tj_t) & (tj_b >= ti_t):
                    ovlap_coords.append([np.max([ti_t, tj_t]), np.min([ti_b, tj_b])])
                    ovlap_units.append([ii,jj])
        # check if any units overlap, if not, return zero
        if len(ovlap_coords) == 0:
            continue

        # get average (max) sedimentation rate for each unit
        avg_sed_rates_max = col_units['max_thick']/(col_units['b_age']-col_units['t_age'])
        avg_sed_rates_min = col_units['min_thick']/(col_units['b_age']-col_units['t_age'])

        # having found overlapping units, (at II 1), pg. 4 ), find unique coordinates of overlapping segments
        unique_coords = np.unique(np.asarray(ovlap_coords))
        # create list with overlap vertices as indices
        ovlaps_verts_units = np.empty(len(unique_coords), dtype=object)
        for ii in range(len(unique_coords)):
            ovlaps_verts_units[ii] = [np.nan] * 1
        for coords, units in zip(ovlap_coords, ovlap_units):
            idx_coords = (unique_coords >= coords[0]) & (unique_coords <= coords[1])
            cur_ovlaps = ovlaps_verts_units[idx_coords]
            for ii in range(len(cur_ovlaps)):
                cur_ovlaps[ii] = cur_ovlaps[ii] + units
            ovlaps_verts_units[idx_coords] = cur_ovlaps
        for ii in range(len(ovlaps_verts_units)):
            # get unique unit indices
            ovlaps_verts_units[ii] = np.unique(np.asarray(ovlaps_verts_units[ii]))
            # get rid of nans
            not_na_idx = np.logical_not(np.isnan(ovlaps_verts_units[ii]))
            ovlaps_verts_units[ii] = ovlaps_verts_units[ii][not_na_idx]

        # now with array with overlapping vertices denoted by overlapping unit indices at each vertex, create
        # vector with entries denoting segments, not endpoints (pg.3 figure at bottom)
        n_ovlap_segs = len(unique_coords)-1
        ovlap_segs_units = np.empty(n_ovlap_segs, dtype=object)
        ovlap_segs_wids = np.zeros(n_ovlap_segs)
        for ii in range(n_ovlap_segs):
            ovlap_segs_units[ii] = list(set(ovlaps_verts_units[ii+1]) & set(ovlaps_verts_units[ii]))
            ovlap_segs_wids[ii] = unique_coords[ii+1] - unique_coords[ii]

        # subtract thinnest units of the overlapping ones (i.e. n-1 slowest avg sed rates times overlap width)
        # loop over vertices
        for ii in range(n_ovlap_segs):
            # get units in current overlapping segment
            cur_units = ovlap_segs_units[ii]
            # for each unit, get thicknesses (sed rates x durations)
            ovlap_thicks_max = np.sort(avg_sed_rates_max.iloc[cur_units].values * ovlap_segs_wids[ii])
            ovlap_thicks_min = np.sort(avg_sed_rates_min.iloc[cur_units].values * ovlap_segs_wids[ii])
            # take sum of smallest/largest thicknesses as maximum/minimum corrections
            corr.loc[col, 'max_corr'] = corr.loc[col, 'max_corr'] + np.sum(ovlap_thicks_max[0:-1])
            corr.loc[col, 'min_corr'] = corr.loc[col, 'min_corr'] + np.sum(ovlap_thicks_min[1:])

    return corr


def seds_by_column_and_age_range(units, columns, t1, t2, method='thickness'):
    """
    units: dataframe of units in which to estimate thicknesses
    t[1] > t[0]
    method='thickness' (default), 'presence' : whether to estimate columnwise thicknesses or just presence/absence a la Peters/Husson.
    returns:
    tmp_df: dataframe indexed by col_id with column geometry. if method='presence', then the singular data column is boolean indicating presence of sediments within the requested time interval. if method='thickness', then two data columns are returned containing the maximum and minimum estimated accumulations of sediment within the requested time interval
    """
    assert t1 > t2, 't1 must be larger than t2'
    assert method in ['thickness', 'presence'], 'invalid value for method'
    # copy data so we don't mess with originals
    units = units.copy()
    columns_tmp = columns.copy()

    # correct for units that span beyond t1 or t2
    w = bin_overlap_correction(units, t1, t2)
    # indices of units within time bounds
    idx_units_in_bounds = w != 0
    # if just looking for presence, then this is enough
    if method=='presence':
        tmp_df = columns_tmp[['geometry']].copy()
        units['presence'] = idx_units_in_bounds.astype(int)
        tmp_df['presence'] = units.groupby(['col_id']).sum(numeric_only=True)['presence']
        tmp_df.loc[tmp_df['presence'] > 0, 'presence'] = 1
        return tmp_df
    elif method=='thickness':
        # get just these units
        w = w[idx_units_in_bounds]
        units = units[idx_units_in_bounds]
        # apply weights
        max_thicks = w * units['max_thick']
        min_thicks = w * units['min_thick']
        units['max_thick'] = max_thicks
        units['min_thick'] = min_thicks
        # sum all seds in each columns
        grouped = units.groupby(['col_id'])
        thicknesses_cols = grouped[['min_thick', 'max_thick']].sum()
        # reset time bounds on units that straddle time bin
        units.loc[units['t_age'] < t2, 't_age'] = t2
        units.loc[units['b_age'] > t1, 'b_age'] = t1
        # get and apply thickness correction
        corr = unit_overlap_correction(grouped.indices, units)
        thicknesses_cols['max_thick'] = thicknesses_cols['max_thick'] - corr['max_corr']
        thicknesses_cols['min_thick'] = thicknesses_cols['min_thick'] - corr['min_corr']
        # combine thicknesses with column geometry
        tmp_df = columns_tmp[['geometry']].merge(thicknesses_cols, left_index=True, right_index=True, how='outer')
        # set nans to zeros
        nan_idx = np.isnan(pd.to_numeric(tmp_df['max_thick']))
        tmp_df.loc[nan_idx, 'max_thick'] = 0
        nan_idx = np.isnan(pd.to_numeric(tmp_df['min_thick']))
        tmp_df.loc[nan_idx, 'min_thick'] = 0
        # make numeric
        tmp_df['max_thick'] = pd.to_numeric(tmp_df['max_thick'])
        tmp_df['min_thick'] = pd.to_numeric(tmp_df['min_thick'])
        return tmp_df


"""
Implementation of York 1969 10.1016/S0012-821X(68)80059-7
IN:
x: mean x-values
y: mean y-values
wx: weights for x-values (typically 1/sigma^2)
wy: weights for y-values (typically 1/sigma^2)
r: correlation coefficient between sigma_x and sigma_y
OUT:
b: maximum likelihood estimate for slope of line
a: maximum likelihood estimate for intercept of lin
b_sig: standard deviation of slope for line
a_sig: standard deviation of intercept for line
mswd: reduced chi-squared statistic for residuals with respect to the maximum likelihood linear model
"""
def yorkfit(x, y, wx, wy, r, thres=1e-3):
    n = len(x)
    # get first guess for b
    b = stats.linregress(x, y)[0]

    # initialize various quantities
    alpha = np.sqrt(wx*wy)

    # now iterate as per manuscript to improve b
    delta = thres+1
    while delta > thres:
        # update values from current value of b
        W = (wx*wy)/(wx + b**2*wy - 2*b*r*alpha)
        xbar = np.sum(x*W)/np.sum(W)
        ybar = np.sum(y*W)/np.sum(W)
        U = x - xbar
        V = y - ybar
        beta = W * (U/wy + b*V/wx - (b*U+V)*r/alpha)
        # update b
        b_new = np.sum(W*beta*V)/np.sum(W*beta*U)
        delta = np.abs(b_new - b)
        b = b_new

    # compute a
    a = ybar - b*xbar
    # compute adjusted x, y
    x_adj = xbar + beta
    x_adj_bar = np.sum(W*x_adj)/np.sum(W)
    u = x_adj - x_adj_bar
    # compute parameter uncertainties
    b_sig = 1/np.sum(W*u**2)
    a_sig = 1/np.sum(W) + x_adj_bar**2*b_sig**2
    # compute goodness of fit (reduced chi-squared statistic)
    mswd = np.sum(W*(y-b*x-a)**2)/(n-2)

    return b, a, b_sig, a_sig, mswd


def findseq(x, val, noteq=False):
    """
    Find sequences of a given value within an input vector.
    IN:
     x: vector of values in which to find sequences
     val: scalar value to find sequences of in x
     noteq: (false) whether to find sequences equal or not equal to the supplied
        value
    OUT:
     idx: array that contains in rows the number of total sequences of val, with
       the first column containing the begin indices of each sequence, the second
       column containing the end indices of sequences, and the third column
       contains the length of the sequence.
    """
    x = x.copy().squeeze()
    assert len(x.shape) == 1, "x must be vector"
    # indices of value in x, and
    # compute differences of x, since subsequent occurences of val in x will
    # produce zeros after differencing. append nonzero value at end to make
    # x and difx the same size
    if noteq:
        validx = np.argwhere(x != val).squeeze()
        x[validx] = val+1
        difx = np.append(np.diff(x),1)
    else:
        validx = np.argwhere(x == val).squeeze()
        difx = np.append(np.diff(x),1)
    nval = len(validx)
    # if val not in x, warn user
    if nval == 0:
        warnings.warn("value val not found in x")
        return 0

    # now, where validx is one and difx is zero, we know that we have
    # neighboring values of val in x. Where validx is one and difx is nonzero,
    # we have end of a sequence

    # now loop over all occurrences of val in x and construct idx
    c1 = 0
    idx = np.empty((1,3))
    while c1 < nval:
        curidx = np.array([[validx[c1],validx[c1],1]])
        c2 = 0
        while difx[validx[c1]+c2] == 0:
            curidx[0,1] += 1
            curidx[0,2] += 1
            c2 += 1
        idx = np.append(idx,curidx,axis=0)
        c1 = c1+c2+1
    idx = idx[1:,:].astype(int)
    return idx


from scipy.signal import windows

"""
multi-taper spectral density estimation using DPSS sequences as described by Percival and Walden.

IN:
y: time series of evenly spaced observations in time
dt: sample spacing in sample coordinate
nw: (default 4) time-halfbandwidth product

OUT:
S_est: estimate for the power spectral density estimation
f: frequency axis

TO DO:
- allow user to specify which/how many frequencies to get
"""
def multitaper(y, dt, nw=4):
    N = len(y)
    W = nw/N
    K = int(2*nw-1)
    # frequencies
    f = freq(N, dt)
    Nf = len(f)
    # time indices
    t = np.arange(N)
    # get discrete prolate spheroidal sequence (dpss) windows
    wins = windows.dpss(N, nw, Kmax=K)
    # get tapered spectra
    Sk = np.zeros((Nf, K))
    for ii in range(K):
        # loop over frequencies
        for ff in range(Nf):
            # compute spectral density estimate
            Sk[ff, ii] = dt*np.abs(np.sum(wins[ii,:]*y*np.exp(-1j*2*np.pi*f[ff]*t*dt)))**2

    # get eigenvalues for N, W
    evals = dpss_evals(N, W)
    # implement adaptive multitaper spectral estimator (Percival and Walden, pg. 370)
    # start with eqn 369a (i.e. no estimate for weights bk)
    K_cur = 1
    S_est = np.sum(np.tile(evals[0:K_cur], (Nf,1))*Sk[:, 0:K_cur], axis=1)/np.sum(evals[0:K_cur])
    bk = np.zeros((Nf, K))
    # make convenient tiled version of eigenvalues
    evals_tile = np.tile(evals[0:K], (Nf,1))
    # iterate over equations 368a and 370a
    for ii in range(5):
        # get weights
        bk = np.tile(S_est, (K, 1)).T/(evals_tile*np.tile(S_est, (K, 1)).T + (1-evals_tile)*np.var(y))
        # update spectrum
        S_est = np.sum(bk**2*evals_tile*Sk[:, 0:K], axis=1)/np.sum(bk**2*evals_tile, axis=1)

    return S_est, f

# get dpss eigenvalues for given N, W
def dpss_evals(N, W):
    t = np.arange(N)
    t1, t2 = np.meshgrid(t, t)
    dt = t1-t2
    dt[np.diag_indices(N)] = 1
    # construct matrix
    A = np.sin(2*np.pi*W*dt)/(np.pi*dt)
    # set diagonal manually (l'hopital)
    A[np.diag_indices(N)] = 2*W
    # compute eigenvalues (should all be real)
    evals = np.real(np.linalg.eig(A)[0])
    # sort by magnitude
    evals[::-1].sort()
    return evals

# generate one-sided frequency axis for given N data and dt sample spacing
def freq(N, dt):
    fs = 1/dt
    fi = fs/N
    fx = np.arange(0, fs/2+fi, fi)
    return fx
