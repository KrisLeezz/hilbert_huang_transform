# -*- coding= utf-8 -*-
import matlab#进行数据转换，将matlab数组转换为python支持的格式
import matlab.engine#调用函数
import csv
import numpy


def readdata(path):
	path=path
	data=[]
	csvfile=open(path,'r')
	csvread=csv.reader(csvfile)
	for line in csvread:
		data.append(line)
	return data

x=readdata('test3.csv')
x=numpy.array(x)
x=x.astype(float)
eng=matlab.engine.start_matlab()
IMF=eng.memd(x)

print (IMF)