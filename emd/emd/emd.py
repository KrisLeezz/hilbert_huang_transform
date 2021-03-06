# -*- coding=utf-8 -*-
#test
if __name__ == "__main__":
	import numpy as np
	from PyEMD import EMD,EEMD,Visualisation
	import matplotlib.pyplot as plt
	import csv

	def emd_decompose(a,t):
		emd=EMD()
		emd.emd(a)
		imfs,res=emd.get_imfs_and_residue()
		plt.plot(t,a)
		plt.title('origin sequence')
		vis=Visualisation(emd)
		vis.plot_imfs(t=t)
		vis.plot_instant_freq(t)
		vis.show()
		plt.show()
		return imfs, res
	def eemd_decompose(a,t):
		eemd = EEMD(trials=200,noise_width=0.4)(a)
		eimfs,res = eemd[:-1],eemd[-1]#这的参数有问题！！！！！！！！！！！！！！！！！！！！！
		vis = Visualisation()
		vis.plot_imfs(imfs=eimfs, residue=res, t=t, include_residue=True)
		vis.plot_instant_freq(t, imfs=eimfs) # 
		vis.show()
		return eimfs,res

	#--------------------------reading data---------------------------------#
	csvfile1=open('51-126.csv')
	csvread1=csv.reader(csvfile1)
	x=[]
	n=0
	for row in csvread1:
		if row[0]!='TAVG':
			x.append(row[0])
			n=n+1
	x=np.array(x)
	x=x.astype(float)
	csvfile1.close()

	y=x

	t=np.arange(0,len(y),1)

	choice='emd'
	if choice == 'emd':
		IMFS=emd_decompose(y,t)
		csvfile2=open('imfs.csv','wb')
		csv_write2=csv.writer(csvfile2)
		for i in range (0,len(IMFS[0]-1)):
			csv_write2.writerow(IMFS[0][i])
		#csv_write2.writerow(IMFS[-1])

		csvfile2.close
		print '-------------------------emd decompose finish-------------------------------'
	elif choice =='eemd':
		EIMFS,ERES=eemd_decompose(y,t)
		csvfile2=open('imfs.csv','wb')
		csv_write2=csv.writer(csvfile2)
		csv_write2.writerows(EIMFS)
		csvfile2.close
		print EIMFS.shape,ERES.shape
		print '-------------------------eemd decompose finish-------------------------------'