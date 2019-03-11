# -*- coding= utf-8 -*-
import numpy
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
	seq,ndir=set_value(x)
	#判断是否满足stop_emd条件
	stop_emd(x,seq,ndir)
def stop_emd(x,seq,ndir):
	print (seq,ndir)

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
		print ('The record is less than dimension')
	else: 
		print('The total number of data is ',q.shape[0])
		print('The dimension is ',q.shape[1])	
	#Specifies maximum number of channels that can be processed by the code
	Max_channels=32
	nbit=0
	maxiterations=1000
	#生成Hammersley Sequence,注意这里需要的是64个方向的位置，所以样本数i=64，
	base=[]#用于变换的基底
	base.append(-ndir)#hammersley sequence的第一个方向是均匀分布i/N
	if (N_dim==2):#2维=圆,3维=2 sphere=球，4维=3 sphere,..., n维=n-1 sphere
		#调用hammersely函数，return seq()
		seq=[]
		interval=2*180/ndir#????用pi表示？？？？？？？？reference-BEMD
		for i in range (0,ndir):
			seq.append(i*interval)
	return seq,ndir

if __name__ == '__main__':
	x=numpy.array([[1,2],[2,3],[3,4]])
	memd(x)
