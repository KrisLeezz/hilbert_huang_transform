# -*- coding= utf-8 -*-
import numpy as np
import math
import matlab
import matlab.engine
import scipy.io
import csv
import matplotlib.pyplot as plt
#the most several important parameters are n_dir,tol,sd,sd2，其中tol的推荐值是0.05，sd2=10*sd.
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

global eng
eng=matlab.engine.start_matlab()
def memd(x,*args):
	#参数初始化,产生方向sampling,seq即产生的方向
	seq,ndir,N,N_dim,stop_crit,sd,sd2,tol,MAXITERATIONS,nbit=set_value(x)
	print ('-----------------------------------------setvalue success-----------------------------------')
	#创建x的索引
	t=np.arange(0,len(x),1)
	if (len(t)!=len(x)):
		print('the index is not equal to the record')
	#print ('the total num_direction is %d'%ndir)#seq就是生成的方向序列
	n_imf=0
	#判断是否需要进行emd分解
	q = np.zeros((N_dim,20,N))
	r=x
	j=1
	while stop_emd(r,seq,ndir,N_dim)!=0:#外侧的大循环只靠是否有充足的极值点，按理说这里应该是判断r是否足够小，但是这个没有标准，只能用这种方法来代替，实际是判断r是不是单调的
		print('----------------------------------the %d emd should continue-------------------------------------------'%j)
		m=r
		#进行IMF的提取
		#判断是否继续进行IMF的提取，即提取的IMF是否满足条件
		if stop_crit=='stop':#三个标准
			stop_sift,env_mean=stop_sifting(m,t,sd,sd2,tol,seq,ndir,N,N_dim)
		#防止当前imf太小,计算精度达不到,这里有点bug不符合常理
		if (np.max(abs(m))<np.max(abs(x)*0.0000000001)):
			if stop_sift==1:
				print ('warning:forced stop of EMD : too small amplitude,but the shifting is continuing')
			else: 
				print('forced stop of EMD : too small amplitude')
			break
		while stop_sift==1 and nbit<MAXITERATIONS:
			m=m-env_mean#1087,4
			if stop_crit=='stop':
				stop_sift,env_mean=stop_sifting(m,t,sd,sd2,tol,seq,ndir,N,N_dim)
			nbit=nbit+1
		m=m.T
		for v in range(0,N_dim):
			#m_imf=m[v].reshape((1,-1))
			if (n_imf<20):
				q[v][n_imf]=m[v].T.reshape((1,-1))
			else:print('please add the size of q')
		#q[n_imf]=m
		n_imf = n_imf+1
		r=r-m.T#提取了IMF的残差
		nbit=0
		j=j+1
	q[:,n_imf,:]=r.T#最后一个imf记录r
	return q
def stop_sifting(x,t,sd,sd2,tol,seq,ndir,N,N_dim):
	print('-------------------------------sifting:calculate env_mean----------------------------------')
	#计算mean
	env_mean,nem,nzm,amp=envelope_mean(x,t,seq,ndir,N,N_dim)
	sx= (np.sum(env_mean**2,axis=1))**0.5#因为这个是个二维的
	if (amp.any()!=0):
		sx=sx.reshape(-1,)
		sx=sx/amp #Flandrin那篇文章为什么要和振幅比？我猜想是因为避免局部有信息，局部没信息，导致IMF提取不完整
	a=sx>sd#注意这里的sd都是取了开方的
	aa=np.mean(a)#python的运算精度没有matlab高
	if (((aa>tol)or any(sx>sd2))and any(nem>2)):#tol相当于文章的1-阿尔法，越大越严格；sd2也是越大越严格,
		stp=1#不满足停止条件
	else:stp=0
	#D=x-mean
	#判断D是否满足IMF条件，是的话IMF=D,return stp_sift=0;否的话R=X-mean，即D，再进行shifting
	return stp,env_mean
def envelope_mean(x,t,seq,ndir,N,N_dim):#计算每个方向的包络线
	#计算包络线的均值和mode振幅的估计
	NBSYM=2
	count=0
	env_mean=np.zeros((len(t),N_dim))
	amp=np.zeros((len(t),))
	nem=np.zeros((ndir,1))
	nzm=np.zeros((ndir,1))
	for it in range (0,ndir):
		#投影
		print('the %d direction'%(it+1))
		y=math.cos(seq[it])*x[:,0]-math.sin(seq[it])*x[:,1]############################################
		#t=np.arange(0,len(y),1)
		#计算投影后信号的极值
		max_pos,max_val,min_pos,min_val,zero=local_peaks(t,y)#返回的下标是按照matlab的形式
		nem[it]=len(max_pos)+len(min_pos)
		nzm[it]=len(zero)
		#边界条件约束
		tmin,tmax,zmin,zmax,mode=boundary_conditions(min_pos,max_pos,t,y,x,NBSYM)
		#极值插值成廓线,三次样条插值
		if (mode):
			tmin=np.array(tmin)#savemat那个函数貌似只支持array，matlab出来的函数必须这么转换一下子
			tmin=tmin.astype(float)
			tmax=np.array(tmax)
			tmax=tmax.astype(float)
			zmin=np.array(zmin)
			zmin=zmin.astype(float)
			zmax=np.array(zmax)
			zmax=zmax.astype(float)
			tt=np.array(t+1)#因为这个变量要传到matlab里
			tt=tt.astype(float)
			scipy.io.savemat('tmin.mat',mdict={'tmin':tmin})
			tminpath='tmin.mat'
			scipy.io.savemat('tamx.mat',mdict={'tamx':tmax})
			tmaxpath='tamx.mat'
			scipy.io.savemat('zmin.mat',mdict={'zmin':zmin})
			zminpath='zmin.mat'
			scipy.io.savemat('zmax.mat',mdict={'zmax':zmax})
			zmaxpath='zmax.mat'
			scipy.io.savemat('tt.mat',mdict={'t':tt})
			tpath='tt.mat'
			#eng=matlab.engine.start_matlab()
			env_min=eng.myspline(tminpath,zminpath,tpath)#
			env_min=(np.array(env_min)).T
			env_max=eng.myspline(tmaxpath,zmaxpath,tpath)#
			env_max=(np.array(env_max)).T
			#计算env_mean
			mean=env_max-env_min#二维的，每个变量的振幅
			mean2=mean**2
			meansum=np.sum(mean2,axis=1)#分量的振幅平方和=原向量振幅的平方
			meansqrt=meansum**0.5
			amp=amp+meansqrt/2#除以2是因为默认是对称的，只取振幅的一半与mean的形状进行比较
			env_mean = env_mean + (env_max+env_min)/2#两个方向各自的mean部分
			print('this is the direction %d th is finished'%it)
		else: count=count+1#有的方向可能没有信息了
	if (ndir>count):#所有方向的阔线和振幅取均值
		env_mean=env_mean/(ndir-count)
		amp=amp/(ndir-count)
	else:print('wrong in envelope_mean')
	return env_mean,nem,nzm,amp
	


def boundary_conditions(indmin,indmax,t,xx,z,nbsym):#镜像对称处理两端的值
	#这里懒得写调用matlab函数
	#python的变量不能直接传送到matlab中
	scipy.io.savemat('indmin.mat',mdict={'indimn':indmin})
	scipy.io.savemat('indmax.mat',mdict={'indmax':indmax})
	scipy.io.savemat('t.mat',mdict={'t':t+1})
	scipy.io.savemat('x.mat',mdict={'x':xx})
	scipy.io.savemat('z.mat',mdict={'z':z})
	scipy.io.savemat('nbsym.mat',mdict={'nbsym':nbsym})
	#eng=matlab.engine.start_matlab()
	timn,tmax,zmin,zmax,mode=eng.myboundary_conditions(nargout=5)
	#eng.exit()
	return timn,tmax,zmin,zmax,mode


def stop_emd(x,seq,ndir,N_dim):#判断是否满足stop_emd条件,当每个方向提取的极值extream小于三个的时候，说明信息不够了，停止进行emd分解
	#对于每一个方向，根据方向sampling的结果进行投影
	print('whether enough extreme exist?')
	ner=np.zeros((ndir,1))
	for it in range(0,ndir):
		#根据生成的投影方向进行投影
		y=math.cos(seq[it])*x[:,0]-math.sin(seq[it])*x[:,1]#################################################
		t=np.arange(0,len(y)+0,1)#创造y的位置索引,用于寻找极值的位置,第一个元素的位置是0
		max_pos,max_val,min_pos,min_val,zero=local_peaks(t,y)#寻找极大值和极小值，返回索引和值，返回的索引下标是从1开始
		ner[it]=len(max_pos)+len(min_pos)#记录每个方向极大值和极小值的数目
	if (ner<3).all():#判断是否有足够的极值进行提取,all()是否全为非零元素，是的话返回1
		stp=0
		print('there is not enough extreme. The EMD should be ended.')
	else:
		stp=1#循环不停止
		print ('Yes')
	#stp=1#测试条件
	return stp

def local_peaks(T, S):
     # 寻找过0的次数
    S1, S2 = S[:-1], S[1:]#制造计算差分的序列
    indzer = np.nonzero(S1*S2<0)[0]
    if np.any(S==0):
        iz = (np.nonzero(S==0)[0])
        if np.any(np.diff(iz)==1):#连续出现两个0,求穿过x轴的点
            zer = (S==0)
            dz = np.diff(np.append(np.append(0, zer), 0))#前后补0，false是不是0（0），true是1原来的值是0
            debz = (np.nonzero(dz==1)[0])
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

    return local_max_pos+1, local_max_val, local_min_pos+1, local_min_val, indzer+1

def set_value(x,num_direction=4,stop_criteria='stop',criteria_value=[0.05,0.5,0.05]):
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
	
	#prm=[2,3,5,7,11,13,17,19,23,29,31,
	#  37,41,43,47,53,59,61,67,71,73,79,
	#  83,89,97,101,103,107,109,113,127,131,137,139,149]
	q=x
	N=q.shape[0]#数据的记录数
	N_dim=q.shape[1]#数据的维数
	if q.shape[0]<q.shape[1]:
		print ('The records is less than dimension')
	else: 
		print('The total number of data is ',q.shape[0])
		print('The dimension is ',q.shape[1])	
	#Specifies maximum number of channels that can be processed by the code
	Max_channels=2
	nbit=0
	maxiterations=1000
	#生成均匀的方向，准备后边投影用
	##生成Hammersley Sequence,注意这里需要的是64个方向的位置，所以样本数i=64，
	#base=[]#用于变换的基底
	#base.append(-ndir)#hammersley sequence的第一个方向是均匀分布i/N
	if (N_dim==2):#2维=圆,3维=2 sphere=球，4维=3 sphere,..., n维=n-1 sphere，
					#但是我不知道hammersley对不同维度的数据处理方式是否一样,这里仅支持2维
		#调用hammersely函数，return seq()
		seq=[]
		interval=2*math.pi/ndir#????用pi表示？应该是默认弧度，这里偷懒了reference-BEMD的方法，按照MEMD应该是hammersely序列，但是2维的h序列也是均匀分布
		for i in range (0,ndir):
			seq.append(i*interval)
	return seq,ndir,N,N_dim,stp_crit,sd,sd2,tol,maxiterations,nbit
def visula(imf):
	numRows=len(imf[0])
	t=np.arange(0,len(imf[0][0])+0,1)
	numCols=2

	plt.figure(12)
	n=1
	ii=0
	for i in range(0,numRows):
		if (imf[ii][i].all()==0):
			break
		else:
			plt.subplot(numRows,numCols,n)
			plt.plot(t,imf[ii][i])
			n=n+1		
			plt.subplot(numRows,numCols,n)
			plt.plot(t,imf[ii+1][i])
			n=n+1
	plt.show()
if __name__ == '__main__':
	x=[]
	csvfile=open('sample.csv','r')
	csvreader=csv.reader(csvfile)
	for i in csvreader:
		x.append(i)
	x=np.array(x)
	x=x.astype(float)
	print('read data successfully')
	IMF=memd(x)
	visula(IMF)
	eng.exit()
