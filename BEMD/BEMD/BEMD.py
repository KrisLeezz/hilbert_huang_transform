# -*- coding= utf-8 -*-
import matlab#进行数据转换，将matlab数组转换为python支持的格式
import matlab.engine#调用函数
import csv
import scipy.io
import numpy

def readdata(path):
	csvfile=open(path+'.csv','r')
	csvread=csv.reader(csvfile)
	data=[]
	for line in csvread:
		data.append(line)
	data=numpy.array(data)
	data=data.astype(float)
	scipy.io.savemat(path+'.mat',mdict={'data':data})
	csvfile.close()
#=========================================================================================
path='test3'#指定文件名，注意不要加上后缀
readdata(path)
eng=matlab.engine.start_matlab()
IMF=eng.myexample(path+'.mat')
#print (type(IMF))#matlab传回来的是一个object,通过numpy转换。注：可以指定返回的参数的个数
IMF=numpy.array(IMF)
#=========================================================================================
for i in range(0,len(IMF)):
	csvfile=open('IMF%i.csv'%i,'w',newline='')
	csvwriter=csv.writer(csvfile)
	csvwriter.writerows(IMF[i])
	csvfile.close()

	

#eng.exit()#加上这句的话调用的matlab的eng会自动关闭
