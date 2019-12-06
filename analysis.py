###################################
import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.mlab as mlab
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import matplotlib.cm as cm
import random

import pandas as pd
###############################################
#options:

readCSV = False
correlate = False
tertiles = True
i6_to_jc_cutoffs = False
i6_to_i12_cutoffs = True
i6_to_c24_cutoffs = True


# Functions
###############################################
def stdDev(x):
	'''function to compute standard deviation'''
	xAvg=np.mean(x)
	xOut=0.
	for k in range(len(x)):
		xOut=xOut+(x[k]-xAvg)**2
	xOut=xOut/(k+1)
	xOut=sqrt(xOut)
	return xOut

def Variance(x):
	'''function to compute the variance (std dev squared)'''
	xAvg=np.mean(x)
	xOut=0.
	for k in range(len(x)):
		xOut=xOut+(x[k]-xAvg)**2
	xOut=xOut/(k+1)
	return xOut

def SumOfSquares(x):
	'''function to compute the sum of squares'''
	xOut=0.
	for k in range(len(x)):
		xOut=xOut+x[k]**2
	return xOut

def corr(x,y):
	''' function to find the correlation of two arrays'''
	xAvg=np.mean(x)
	Avgy=np.mean(y)
	rxy=0.
	n=min(len(x),len(y))
	if(n==0):
	    print(x)
	    print(y)
	for k in range(n):
		rxy=rxy+(x[k]-xAvg)*(y[k]-Avgy)
	rxy=rxy/(k+1)
	stdDevx=np.std(x)
	stdDevy=np.std(y)
	rxy=rxy/(stdDevx*stdDevy)
	return rxy

def clean(x):
    for i in range(0, len(x)):
        if x[i]!='$$.$' and x[i]!='NA' and x[i]!='$' and x[i]!='$$' and x[i]!='999':
            x[i] = float(x[i])

def makeMask(x1,x2,mask):
    for i in range(0, len(x1)):
        if(x1[i]!=-999.0 and x2[i]!=-999.0):
            mask[i]=0

def correlateSurveys(KEY1, KEY2):
    for firstList in KEY1:
        for secondList in KEY2:
            List1 = dataDictionary[firstList]
            List2 = dataDictionary[secondList]
            mask=np.ones(shape=[len(dataDictionary[firstList])])

            makeMask(List1,List2,mask)

            tmp = np.asarray(List1)
            tmp2 = np.asarray(List2)

            tmp = np.ma.masked_array(tmp, mask)
            tmp2 = np.ma.masked_array(tmp2, mask)

            if(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2))>0.4):
                print(firstList)
                print(secondList)
                print(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2)))
####################################################

wddata = ('/Users/garyk/Documents/code/Jivita/data/')

if readCSV:
	reader = csv.DictReader(open(wddata+'garyk1.csv','r'))

	result = {}
	for row in reader:
	    for column, value in row.items():
	        result.setdefault(column, []).append(value)

	pickle.dump(result, open(wddata+"save.p", "wb" ))

dataDictionary = pickle.load(open(wddata+"save.p", "rb"))

nSurveys = len(dataDictionary['jcslphngry'])

passed = 0
failed = 0


JCsurvey = []
I6survey = []
I12survey = []
C24survey = []
SSsurvey=[]
M3survey=[]

for key in dataDictionary.keys():
    if ("date" not in key) and ("id" not in key) and ("age" not in key):
        if(key[:2]=="jc" and key!="JCHHMEMB"):
            JCsurvey.append(key)
        elif(key[:2]=="i6"):
            I6survey.append(key)
        elif(key[:3]=="i12"):
            I12survey.append(key)
        elif(key[:3]=="c24"):
            C24survey.append(key)
        elif(key[:2]=="ss" and key!="ssmicros"):
            SSsurvey.append(key)
        elif(key[:2]=="m3"):
            M3survey.append(key)

i6dict = {}
for key in I6survey:
	i6dict[key]=[]
	for i in range(len(dataDictionary[key])):
		if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"):
			dataDictionary[key][i]='-999'
		i6dict[key].append(float(dataDictionary[key][i]))
		dataDictionary[key][i]=float(dataDictionary[key][i])

i12dict = {}
for key in I12survey:
	i12dict[key]=[]
	for i in range(len(dataDictionary[key])):
		if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"):
			dataDictionary[key][i]='-999'
		i12dict[key].append(float(dataDictionary[key][i]))
		dataDictionary[key][i]=float(dataDictionary[key][i])

c24dict = {}
for key in C24survey:
	c24dict[key]=[]
	for i in range(len(dataDictionary[key])):
		if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"):
			dataDictionary[key][i]='-999'
		c24dict[key].append(float(dataDictionary[key][i]))
		dataDictionary[key][i]=float(dataDictionary[key][i])

jcdict = {}
for key in JCsurvey:
	jcdict[key]=[]
	for i in range(len(dataDictionary[key])):
		if(key!="jcmage" and dataDictionary[key][i]=='9'):
			dataDictionary[key][i]='-999'
		if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"):
			dataDictionary[key][i]='-999'
		jcdict[key].append(float(dataDictionary[key][i]))
		dataDictionary[key][i]=float(dataDictionary[key][i])

for key in SSsurvey:
    for i in range(len(dataDictionary[key])):
        if(key in ['ssclass','sshclass','ssempn','sshempn','sshempn2','sscattle','ssgoats','sschickens','ssducks','ssmango','sspapaya','ssbanana','sscoconut','ssjfruit','ssricest']):
            if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"  or dataDictionary[key][i]=="$$" or dataDictionary[key][i]=="99"):
                dataDictionary[key][i]='-999'
        elif(key in ['sscropnu','ssgrovenu','sspondnu','sshomenu','sslandnu','ssgardennu']):
            if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"  or dataDictionary[key][i]=="$$" or dataDictionary[key][i]=="999"):
                dataDictionary[key][i]='-999'
        elif(key == 'SSLOANAMT'):
            if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"  or dataDictionary[key][i]=="$$" or dataDictionary[key][i]=="999999"):
                dataDictionary[key][i]='-999'
        else:
            if(dataDictionary[key][i]=="NA" or dataDictionary[key][i]=="$"  or dataDictionary[key][i]=="$$"  or dataDictionary[key][i]=="9"):
                dataDictionary[key][i]='-999'
        dataDictionary[key][i]=float(dataDictionary[key][i])

if correlate:
	for firstList in JCsurvey:
	    for secondList in SSsurvey:
	        List1 = dataDictionary[firstList]
	        List2 = dataDictionary[secondList]
	        mask=np.ones(shape=[len(List1)])

	        makeMask(List1,List2,mask)

	        tmp = np.asarray(List1)
	        tmp2 = np.asarray(List2)

	        tmp = np.ma.masked_array(tmp, mask)
	        tmp2 = np.ma.masked_array(tmp2, mask)

	        if(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2))>0.4):
	            print(firstList)
	            print(secondList)
	            print(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2)))

if tertiles:
	jcPD = pd.DataFrame(jcdict)
	i6PD = pd.DataFrame(i6dict)
	i12PD = pd.DataFrame(i12dict)
	c24PD = pd.DataFrame(c24dict)

	cleani6=np.where((i6PD> -1).all(1))
	cleani12=np.where((i12PD> -1).all(1))
	cleanc24=np.where((c24PD> -1).all(1))
	cleanjc=np.where((jcPD> -1).all(1))

	#compare 6 month and 12 month survey results
	if i6_to_i12_cutoffs:
		intersection = np.intersect1d(cleani6, cleani12)
		i6PD_T=i6PD.T
		i12PD_T=i12PD.T

		i6scores = np.zeros(intersection.shape)
		i12scores = np.zeros(intersection.shape)

		counter = -1
		for i in intersection:
			counter += 1
			for key in I6survey:
				if(key == "i6fssqmls"):
					i6scores[counter] -= i6PD_T[i][key]
				else:
					i6scores[counter] += i6PD_T[i][key]
			for key in I12survey:
				if(key == "i12fssqmls"):
					i12scores[counter] -= i12PD_T[i][key]
				else:
					i12scores[counter] += i12PD_T[i][key]

		scores_i6_pd = pd.DataFrame(data = i6scores[:])
		print(scores_i6_pd.keys())
		# scores_i6_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()

		scores_i12_pd = pd.DataFrame(data = i12scores[:])
		# scores_i12_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()

		secure_i6 = np.where(scores_i6_pd[0] == 4)[0]
		insecure_i6 =  np.where(scores_i6_pd[0] > 4)[0]

		secure_i12 = np.where(scores_i12_pd[0] == 4)[0]
		insecure_i12 =  np.where(scores_i12_pd[0] > 4)[0]

		intersection = np.intersect1d(secure_i6, secure_i12)

		accuracy = intersection.shape[0]/len(secure_i6)
		print("i6 secure that are in i12: "+ str(accuracy))

		insecure_i6_pd = scores_i6_pd[scores_i6_pd[0] > 4]
		i6_cutoffs = insecure_i6_pd.quantile([0.33,0.66])

		insecure_i12_pd = scores_i12_pd[scores_i12_pd[0] > 4]
		i12_cutoffs = insecure_i12_pd.quantile([0.33,0.66])

		i6_insecure_s = np.where(insecure_i6_pd <= i6_cutoffs[0][0.33])[0]
		i12_insecure_s = np.where(insecure_i12_pd <= i12_cutoffs[0][0.33])[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])
		print(i6_cutoffs)
		print(i12_cutoffs)
		insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
		print(insecure_accuracy)

		i6_insecure_s = np.where((insecure_i6_pd > i6_cutoffs[0][0.33]) & (insecure_i6_pd <= i6_cutoffs[0][0.66]))[0]
		i12_insecure_s = np.where((insecure_i12_pd > i12_cutoffs[0][0.33]) & (insecure_i12_pd <= i12_cutoffs[0][0.66]))[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])

		insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
		print(insecure_accuracy)

		i6_insecure_s = np.where(insecure_i6_pd > i6_cutoffs[0][0.66])[0]
		i12_insecure_s = np.where(insecure_i12_pd > i12_cutoffs[0][0.66])[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])

		insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
		print(insecure_accuracy)

		print(len(i6_insecure_s))

	if i6_to_c24_cutoffs:
	  intersection = np.intersect1d(cleani6, cleanc24)
	  i6PD_T=i6PD.T
	  c24PD_T=c24PD.T

	  i6scores = np.zeros(intersection.shape)
	  c24scores = np.zeros(intersection.shape)

	  counter = -1
	  for i in intersection:
	    counter += 1
	    for key in I6survey:
	      if(key == "i6fssqmls"):
	        i6scores[counter] -= i6PD_T[i][key]
	      else:
	        i6scores[counter] += i6PD_T[i][key]
	    for key in C24survey:
	      if(key == "c24fssqmls"):
	        c24scores[counter] -= c24PD_T[i][key]
	      else:
	        c24scores[counter] += c24PD_T[i][key]

	  scores_i6_pd = pd.DataFrame(data = i6scores[:])
	  print(scores_i6_pd.keys())
	  # scores_i6_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()

	  scores_c24_pd = pd.DataFrame(data = c24scores[:])
	  # scores_c24_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()

	  secure_i6 = np.where(scores_i6_pd[0] == 4)[0]
	  insecure_i6 =  np.where(scores_i6_pd[0] > 4)[0]

	  secure_c24 = np.where(scores_c24_pd[0] == 4)[0]
	  insecure_c24 =  np.where(scores_c24_pd[0] > 4)[0]

	  intersection = np.intersect1d(secure_i6, secure_c24)

	  accuracy = intersection.shape[0]/len(secure_i6)
	  print("i6 secure that are in c24: "+ str(accuracy))

	  insecure_i6_pd = scores_i6_pd[scores_i6_pd[0] > 4]
	  i6_cutoffs = insecure_i6_pd.quantile([0.33,0.66])

	  insecure_c24_pd = scores_c24_pd[scores_c24_pd[0] > 4]
	  c24_cutoffs = insecure_c24_pd.quantile([0.33,0.66])

	  i6_insecure_s = np.where(insecure_i6_pd <= i6_cutoffs[0][0.33])[0]
	  c24_insecure_s = np.where(insecure_c24_pd <= c24_cutoffs[0][0.33])[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])
	  print(i6_cutoffs)
	  print(c24_cutoffs)
	  insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
	  print(insecure_accuracy)

	  i6_insecure_s = np.where((insecure_i6_pd > i6_cutoffs[0][0.33]) & (insecure_i6_pd <= i6_cutoffs[0][0.66]))[0]
	  c24_insecure_s = np.where((insecure_c24_pd > c24_cutoffs[0][0.33]) & (insecure_c24_pd <= c24_cutoffs[0][0.66]))[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])

	  insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
	  print(insecure_accuracy)

	  i6_insecure_s = np.where(insecure_i6_pd > i6_cutoffs[0][0.66])[0]
	  c24_insecure_s = np.where(insecure_c24_pd > c24_cutoffs[0][0.66])[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])

	  insecure_accuracy = intersection.shape[0]/len(i6_insecure_s)
	  print(insecure_accuracy)
	  print(len(i6_insecure_s))

	if i6_to_jc_cutoffs:
		intersection = np.intersect1d(cleani6, cleanjc)
		jcPD_T=jcPD.T
		i6PD_T=i6PD.T

		i6scores = np.zeros(intersection.shape)
		jcscores = np.zeros(intersection.shape)

		counter = -1
		for i in intersection:
			counter += 1
			for key in I6survey:
				if(key == "i6fssqmls" or key == "i6fsrice"):
					i6scores[counter] -= i6PD_T[i][key]
				else:
					i6scores[counter] += i6PD_T[i][key]
			for key in JCsurvey:
				if(key == "jcswor" or key == "jcnofood" or key == "jcslphngry" or key == "jcdayhngry"):
					jcscores[counter] -= jcPD_T[i][key]
				else:
					jcscores[counter] += jcPD_T[i][key]

		scores_i6_pd = pd.DataFrame(data = i6scores[:])
		scores_jc_pd = pd.DataFrame(data = jcscores[:])

		secure_i6_pd = scores_i6_pd[scores_i6_pd[0] == -2]
		insecure_i6_pd =  scores_i6_pd[scores_i6_pd[0] > -2]

		i6_cutoffs = scores_i6_pd.quantile([0.33,0.66])
		jc_cutoffs = scores_jc_pd.quantile([0.33,0.66])

		i6_insecure = np.where(scores_i6_pd < i6_cutoffs[0][0.33])
		jc_insecure = np.where(scores_jc_pd < jc_cutoffs[0][0.33])

		scores_i6_pd = insecure_i6_pd

		intersection = np.intersect1d(i6_insecure[0], jc_insecure)

		insecure_accuracy = intersection.shape[0]/scores_i6_pd.shape[0]
		print(insecure_accuracy)
		#
		# insecure_accuracy = intersection.shape[0]/scores_jc_pd.shape[0]
		# print(insecure_accuracy)

		i6_insecure = np.where((scores_i6_pd >= i6_cutoffs[0][0.33]) & (scores_jc_pd <= i6_cutoffs[0][0.66]))
		jc_insecure = np.where((scores_jc_pd >= jc_cutoffs[0][0.33]) & (scores_jc_pd <= jc_cutoffs[0][0.66]))

		intersection = np.intersect1d(i6_insecure, jc_insecure)

		insecure_accuracy = intersection.shape[0]/len(i6_insecure[0])
		print(insecure_accuracy)

		insecure_accuracy = intersection.shape[0]/len(jc_insecure[0])
		print(insecure_accuracy)

		i6_insecure = np.where(scores_i6_pd > i6_cutoffs[0][0.66])
		jc_insecure = np.where(scores_jc_pd > jc_cutoffs[0][0.66])

		intersection = np.intersect1d(i6_insecure, jc_insecure)

		insecure_accuracy = intersection.shape[0]/len(i6_insecure[0])
		print(insecure_accuracy)
		print(len(i6_insecure_s))

		# insecure_accuracy = intersection.shape[0]/len(jc_insecure[0])
		# print(insecure_accuracy)

	# for key in I6survey:
	# 	i6PD[key]
	# print(pd.Series(i6PD(I6survey)==0))
	# for i in range(len(dataDictionary['sectorid_j3'])):
	# 	print(i)
	# 	print(jcPD[i])
	# 	if((jcPD[i] > -1).all(jcPD.dtypes) and (i6PD[i] > -1).all(i6PD.dtypes)):
	# 		user_ids.append(i)
	#
	# ci6PD = jc_i6_pd.where()
	# cjcPD = jcPD.where(i6PD[I6survey]> -1 and i6PD > -1)
