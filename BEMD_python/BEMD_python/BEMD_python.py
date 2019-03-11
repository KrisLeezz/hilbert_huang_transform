# -*- coding= utf-8 -*-
import numpy as np
import math
## This program was designed for BEMD, the reference matlab code is from ‘Rehman and Mandic’
#There are 4 paramters at most:
#	X,your bivariable data;
#	num_directions,the number of directions you set to project, num_directions=64 is default.
#	'stop_criteria',your stop criteria,two were approved:
#		-'stop', which uses the standard stooping criterion specified in reference [2]. 'stop' is default.
#		-'fix_h', which uses the modified version of the stooping criteria specified in reference[3]
#	criteria_value, this is decided by stop criteria
#		-'stop' is correspond to three threshold and tolerance values, see reference[2]
#		the defult three value is 0.075,0.75,0.075
#		-'fix_h' is correspond to n_iter, which is an integer specifying the number of consecutive iterations.
#		when the number of extrema and the number of zero crossings differ at most by one, see reference[3]
#references:
	#[1]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010
	#[2]  G. Rilling, P. Flandrin and P. Goncalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
	#     Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003
	#[3]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
	#     Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003
def memd(x,*args):
	#参数初始化
	seq,ndir,N,N_dim,stop_crit,sd,sd2,tol=set_value(x)
	#创建x的索引
	t=np.arange(0,len(x),1)
	print ('the total num_direction is %d'%ndir)#seq就是生成的方向序列
	n_imf=1
	while stop_emd(x,seq,ndir,N_dim)!=0:
		if stop_crit=='stop':
			stop_sifting(x,t,sd,sd2,tol,seq,ndir,N,N_dim)
def stop_sifting(x,t,sd,sd2,tol,seq,ndir,N,N_dim):#停止提取一个IMF循环的停止条件
	envelope_mean(x,t,seq,ndir,N,N_dim)
def envelope_mean(x,t,seq,ndir,N,N_dim):
	#计算包络线的均值和mode振幅的估计
	NBSYM=2
	count=0
	env_mean=np.zeros((len(t),N_dim))
	amp=np.zeros((len(t),1))
	nem=np.zeros((ndir,1))
	nzm=np.zerps((ndir,1))
	for it in range (0,64):
		print(it)
		#获取hammersely 序列对应的角度
		#投影
		#计算投影后信号的极值
		#边界条件约束
		#极值插值成廓线
	#所有方向的阔线取均值
def stop_emd(x,seq,ndir,N_dim):#判断是否满足stop_emd条件,当提取的极值extream小于三个的时候，停止第一层循环
	#对于每一个方向，变换hammeresely序列为角度获得值
	for it in range(0,N_dim):
		#根据生成的投影方向进行投影
		#t=np.arange(0,len(x),1)#创造x的位置索引
		max_pos,max_val,min_pos,min_val,zero=local_peaks(t,x)#寻找极大值和极小值
		ner[it]=len(max_pos)+len(min_pos)#记录每个方向极大值和极小值的数目
	if ner.all()<3:#判断是否停止循环
		stp=0#循环停止
	else:stp=1#循环不停止
	#stp=1#测试条件
	return stp

def local_peaks(T, X):#
    """
	a sampler method, another is to find by parabol, I will add it in the future.
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
    indmin = np.nonzero(np.r_[d1*d2<0] & np.r_[d1<0])[0]+1
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

def set_value(x,num_direction=64,stop_criteria='stop',criteria_value=[0.075,0.75,0.075]):
	#validate the number of parameter!!!!!!!!!!!!!!!!!!!但是这个好像没有什么作用
	if (x.any()==None or num_direction==None or stop_criteria==None or criteria_value==None ):
		print ('you miss some parameter')
	else: print('The parameter is OK')
	
	ndir=num_direction#number of projecting direction
	stp_crit=stop_criteria#迭代停止条件
	stp_val=criteria_value
	stp_vec=[]# 'stop'条件下的阈值
	stp_cnt=[]# 'fix_h'条件下的阈值
	sd=[]
	sd2=[]
	tol=[]
	if (stp_crit=='stop' and len(stp_val)==3):
		stp_vec=stp_val
		#将三个值赋给sd,sd2,tol
		sd=stp_val[0]
		sd2=stp_val[1]
		tol=stp_val[2]
	elif (stp_crit=='fix_h' and len(stp_val)==1):
		stp_cnt=stp_val#注意这个值是数组，这里暂时没有考虑‘fix_h’条件下的默认值
	else: print ('there is a misktake of stop criteria and stop value')
	maxiterations=[]#最大的迭代条件
	
	prm=[2,3,5,7,11,13,17,19,23,29,31,
	  37,41,43,47,53,59,61,67,71,73,79,
	  83,89,97,101,103,107,109,113,127,131,137,139,149]
	q=x
	N=q.shape[0]#数据的记录数
	N_dim=q.shape[1]#数据的维数
	if q.shape[0]<q.shape[1]:
		print ('The records is less than dimension')
	else: 
		print('The total number of data is ',q.shape[0])
		print('The dimension is ',q.shape[1])	
	#Specifies maximum number of channels that can be processed by the code
	Max_channels=32
	nbit=0
	maxiterations=1000
	#生成均匀的方向，准备后边投影用
	#生成Hammersley Sequence,注意这里需要的是64个方向的位置，所以样本数i=64，
	base=[]#用于变换的基底
	base.append(-ndir)#hammersley sequence的第一个方向是均匀分布i/N
	if (N_dim==2):#2维=圆,3维=2 sphere=球，4维=3 sphere,..., n维=n-1 sphere，
					#但是我不知道hammersley对不同维度的数据处理方式是否一样
		#调用hammersely函数，return seq()
		seq=[]
		interval=2*math.pi/ndir#????用pi表示？应该是默认弧度，这里偷懒了reference-BEMD的方法，按照MEMD应该是hammersely序列
		for i in range (0,ndir):
			seq.append(i*interval)
	return seq,ndir,N,N_dim,stp_crit,sd,sd2,tol

if __name__ == '__main__':
	x=np.array([[1,2],[2,3],[3,4]])
	memd(x)
