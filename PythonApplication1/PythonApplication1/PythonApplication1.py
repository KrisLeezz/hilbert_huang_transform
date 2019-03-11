# -*- coding=utf-8 -*-
import numpy as np

def find_extrema_simple(T, S):
    """
    Performs extrema detection, where extremum is defined as a point,
    that is above/below its neighbours.

    See :meth:`EMD.find_extrema`.
    """

    # Finds indexes of zero-crossings
    S1, S2 = S[:-1], S[1:]
    indzer = np.nonzero(S1*S2<0)[0]
    if np.any(S==0):
        iz = np.nonzero(S==0)[0]
        if np.any(np.diff(iz)==1):
            zer = (S==0)
            dz = np.diff(np.append(np.append(0, zer), 0))
            debz = np.nonzero(dz==1)[0]
            finz = np.nonzero(dz==-1)[0]-1
            indz = np.round((debz+finz)/2.)
        else:
            indz = iz

        indzer = np.sort(np.append(indzer, indz))

    # Finds local extrema
    d = np.diff(S)
    d1, d2 = d[:-1], d[1:]
    indmin = np.nonzero(np.r_[d1*d2<0] & np.r_[d1<0])[0]+1#获得非零索引的位置
    indmax = np.nonzero(np.r_[d1*d2<0] & np.r_[d1>0])[0]+1

    # When two or more points have the same value
    if np.any(d==0):

        imax, imin = [], []

        bad = (d==0)
        dd = np.diff(np.append(np.append(0, bad), 0))
        debs = np.nonzero(dd==1)[0]
        fins = np.nonzero(dd==-1)[0]
        if debs[0] == 1:
            if len(debs)>1:
                debs, fins = debs[1:], fins[1:]
            else:
                debs, fins = [], []

        if len(debs) > 0:
            if fins[-1] == len(S)-1:
                if len(debs) > 1:
                    debs, fins = debs[:-1], fins[:-1]
                else:
                    debs, fins = [], []

        lc = len(debs)
        if lc > 0:
            for k in range(lc):
                if d[debs[k]-1] > 0:
                    if d[fins[k]] < 0:
                        imax.append(np.round((fins[k]+debs[k])/2.))
                else:
                    if d[fins[k]] > 0:
                        imin.append(np.round((fins[k]+debs[k])/2.))

        if len(imax) > 0:
            indmax = indmax.tolist()
            for x in imax: indmax.append(int(x))
            indmax.sort()

        if len(imin) > 0:
            indmin = indmin.tolist()
            for x in imin: indmin.append(int(x))
            indmin.sort()

    local_max_pos = T[indmax]
    local_max_val = S[indmax]
    local_min_pos = T[indmin]
    local_min_val = S[indmin]

    return local_max_pos, local_max_val, local_min_pos, local_min_val, indzer

s=np.array([2,0,5,3,1,8,8,7])
t=np.arange(0,len(s),1)
max_pos,max_val,min_pos,min_val,zero=find_extrema_simple(t,s)
print ('')

