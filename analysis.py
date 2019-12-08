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
#Configuration options:

readCSV = False
correlate = False

tertiles = True #true = calculate the tertiles
Rquestion = True #true = remove the fsmortg indicator

i6_to_jc_cutoffs = True #true = calculate the tertiles for i6 to jc
i12_to_jc_cutoffs = True #true = calculate the tertiles for i12 to jc

i6_to_i12_cutoffs = True #true = calculate the tertiles for i6 to i12
i6_to_c24_cutoffs = True #true = calculate the tertiles for i6 to c24

#data folder source
wddata = ('/Users/garyk/Documents/code/Jivita/data/')

# useful functions:
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

#reads from the original data file, and creates a dictionary
if readCSV:
	reader = csv.DictReader(open(wddata+'garyk1.csv','r'))

	result = {}
	for row in reader:
	    for column, value in row.items():
	        result.setdefault(column, []).append(value)

	pickle.dump(result, open(wddata+"save.p", "wb" ))

dataDictionary = pickle.load(open(wddata+"save.p", "rb"))
#get total length for other purposes
nSurveys = len(dataDictionary['jcslphngry'])

#seperate keys based on their survey
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

#using the corresponding keys in the dictionary, create 'cleaned' dictionaries by survey (only keeping valid survey responses)
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

# if correlate is true, will correlate each of the indicators
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

	        if(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2))>0.5):
	            print(firstList)
	            print(secondList)
	            print(corr(np.ma.compressed(tmp),np.ma.compressed(tmp2)))

# if tertiles is true, will create tertiles of data
if tertiles:
	#pandas dataframes made from dictionaries
	jcPD = pd.DataFrame(jcdict)
	i6PD = pd.DataFrame(i6dict)
	i12PD = pd.DataFrame(i12dict)
	c24PD = pd.DataFrame(c24dict)
	#pandas dataframes cleaned (remove -999's)
	cleani6=np.where((i6PD> -1).all(1))
	cleani12=np.where((i12PD> -1).all(1))
	cleanc24=np.where((c24PD> -1).all(1))
	cleanjc=np.where((jcPD> -1).all(1))

	#compare 6 month and 12 month survey results
	if i6_to_i12_cutoffs:

		#get the intersection: those who responded to both surveys
		intersection = np.intersect1d(cleani6, cleani12)
		i6PD_T=i6PD.T
		i12PD_T=i12PD.T

		i6scores = np.zeros(intersection.shape)
		i12scores = np.zeros(intersection.shape)

		#set the score for what a food secure score is (depends on how many questions are used based on the Rquestion config)
		if(Rquestion):
			secure_score = 3
		else:
			secure_score = 4

		#summing each indicator to get the household food insecurity score for each
		counter = -1
		for i in intersection:
			counter += 1
			if(Rquestion):
				for key in I6survey:
					if(key != "i6fsmortg"):
						if(key == "i6fssqmls"):
							i6scores[counter] -= i6PD_T[i][key]
						else:
							i6scores[counter] += i6PD_T[i][key]
				for key in I12survey:
					if(key != "i12fsmortg"):
						if(key == "i12fssqmls"):
							i12scores[counter] -= i12PD_T[i][key]
						else:
							i12scores[counter] += i12PD_T[i][key]
			else:
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
		# scores_i6_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()
		scores_i12_pd = pd.DataFrame(data = i12scores[:])
		# scores_i12_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()

		#seperate food secure and insecure by HFI score
		secure_i6 = np.where(scores_i6_pd[0] == secure_score)[0]
		insecure_i6 =  np.where(scores_i6_pd[0] > secure_score)[0]

		secure_i12 = np.where(scores_i12_pd[0] == secure_score)[0]
		insecure_i12 =  np.where(scores_i12_pd[0] > secure_score)[0]

		#determine number classified consistently as secure between i6 and i12
		intersection = np.intersect1d(secure_i6, secure_i12)

		secure = intersection.shape[0]/len(secure_i6)

		#determine number consistently as insecure between i6 and i12
		intersection = np.intersect1d(insecure_i6, insecure_i12)
		insecure = intersection.shape[0]/len(insecure_i6)

		#taking the food insecure group, get tertile scores
		insecure_i6_pd = scores_i6_pd[scores_i6_pd[0] > secure_score]
		i6_cutoffs = insecure_i6_pd.quantile([0.33,0.66])

		insecure_i12_pd = scores_i12_pd[scores_i12_pd[0] > secure_score]
		i12_cutoffs = insecure_i12_pd.quantile([0.33,0.66])

		#find mild insecure tertile, calculate consistency between i6 and i12
		i6_insecure_s = np.where(insecure_i6_pd <= i6_cutoffs[0][0.33])[0]
		i12_insecure_s = np.where(insecure_i12_pd <= i12_cutoffs[0][0.33])[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])

		insecure_accuracy1 = intersection.shape[0]/len(i6_insecure_s)

		#find moderate insecure tertile, calculate consistency between i6 and i12
		i6_insecure_s = np.where((insecure_i6_pd > i6_cutoffs[0][0.33]) & (insecure_i6_pd <= i6_cutoffs[0][0.66]))[0]
		i12_insecure_s = np.where((insecure_i12_pd > i12_cutoffs[0][0.33]) & (insecure_i12_pd <= i12_cutoffs[0][0.66]))[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])

		insecure_accuracy2 = intersection.shape[0]/len(i6_insecure_s)

		#find severe insecure tertile, calculate consistency between i6 and i12
		i6_insecure_s = np.where(insecure_i6_pd > i6_cutoffs[0][0.66])[0]
		i12_insecure_s = np.where(insecure_i12_pd > i12_cutoffs[0][0.66])[0]

		intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_i12[i12_insecure_s])

		insecure_accuracy3 = intersection.shape[0]/len(i6_insecure_s)

		#write results to csv
		with open(wddata + 'i6_to_i12_comparison', 'w') as writefile:
			writefile  = csv.writer(writefile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			writefile.writerow(['category', 'percent', 'cutoff_i6', 'cutoff_i12'])
			writefile.writerow(['Food_secure', secure])
			writefile.writerow(['Food_insecure', insecure])
			writefile.writerow(['tertile1', insecure_accuracy1, i6_cutoffs[0][0.33], i12_cutoffs[0][0.33]])
			writefile.writerow(['tertile2', insecure_accuracy2, i6_cutoffs[0][0.66], i12_cutoffs[0][0.66]])
			writefile.writerow(['tertile3', insecure_accuracy3])

	#compare 6 month and 24 month survey results
	if i6_to_c24_cutoffs:
	  intersection = np.intersect1d(cleani6, cleanc24)
	  i6PD_T=i6PD.T
	  c24PD_T=c24PD.T

	  #get the intersection: those who responded to both surveys
	  i6scores = np.zeros(intersection.shape)
	  c24scores = np.zeros(intersection.shape)

	  #set the score for what a food secure score is (depends on how many questions are used based on the Rquestion config)
	  if(Rquestion):
		  secure_score = 3
	  else:
		  secure_score = 4

	  #summing each indicator to get the household food insecurity score for each
	  counter = -1
	  for i in intersection:
		  counter += 1
		  if(Rquestion):
			  for key in I6survey:
				  if(key!='i6fsmortg'):
					  if(key == "i6fssqmls"):
						  i6scores[counter] -= i6PD_T[i][key]
					  else:
						  i6scores[counter] += i6PD_T[i][key]
			  for key in C24survey:
				  if(key!='c24fsmortg'):
					  if(key == "c24fssqmls"):
						  c24scores[counter] -= c24PD_T[i][key]
					  else:
						  c24scores[counter] += c24PD_T[i][key]
		  else:
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
	  # scores_i6_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()

	  scores_c24_pd = pd.DataFrame(data = c24scores[:])
	  # scores_c24_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()

	  #seperate food secure and insecure by HFI score
	  secure_i6 = np.where(scores_i6_pd[0] == secure_score)[0]
	  insecure_i6 =  np.where(scores_i6_pd[0] > secure_score)[0]

	  secure_c24 = np.where(scores_c24_pd[0] == secure_score)[0]
	  insecure_c24 =  np.where(scores_c24_pd[0] > secure_score)[0]

	  #determine number classified consistently as secure between i6 and c24
	  intersection = np.intersect1d(secure_i6, secure_c24)
	  secure = intersection.shape[0]/len(secure_i6)

	  #determine number consistently as insecure between i6 and c24
	  intersection = np.intersect1d(insecure_i6, insecure_c24)
	  insecure = intersection.shape[0]/len(insecure_i6)

	  #taking the food insecure group, get tertile scores
	  insecure_i6_pd = scores_i6_pd[scores_i6_pd[0] > secure_score]
	  i6_cutoffs = insecure_i6_pd.quantile([0.33,0.66])

	  insecure_c24_pd = scores_c24_pd[scores_c24_pd[0] > secure_score]
	  c24_cutoffs = insecure_c24_pd.quantile([0.33,0.66])

	  #find mild insecure tertile, calculate consistency between i6 and c24
	  i6_insecure_s = np.where(insecure_i6_pd <= i6_cutoffs[0][0.33])[0]
	  c24_insecure_s = np.where(insecure_c24_pd <= c24_cutoffs[0][0.33])[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])

	  insecure_accuracy1 = intersection.shape[0]/len(i6_insecure_s)

	  #find moderate insecure tertile, calculate consistency between i6 and i12
	  i6_insecure_s = np.where((insecure_i6_pd > i6_cutoffs[0][0.33]) & (insecure_i6_pd <= i6_cutoffs[0][0.66]))[0]
	  c24_insecure_s = np.where((insecure_c24_pd > c24_cutoffs[0][0.33]) & (insecure_c24_pd <= c24_cutoffs[0][0.66]))[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])

	  insecure_accuracy2 = intersection.shape[0]/len(i6_insecure_s)

	  #find severe insecure tertile, calculate consistency between i6 and i12
	  i6_insecure_s = np.where(insecure_i6_pd > i6_cutoffs[0][0.66])[0]
	  c24_insecure_s = np.where(insecure_c24_pd > c24_cutoffs[0][0.66])[0]

	  intersection = np.intersect1d(insecure_i6[i6_insecure_s], insecure_c24[c24_insecure_s])

	  insecure_accuracy3 = intersection.shape[0]/len(i6_insecure_s)

	#write results to csv
	with open(wddata + 'i6_to_c24_comparison', 'w') as writefile:
		writefile  = csv.writer(writefile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		writefile.writerow(['category', 'percent', 'cutoff_i6' , 'cutoff_c24'])
		writefile.writerow(['Food_secure', secure])
		writefile.writerow(['Food_insecure', insecure])
		writefile.writerow(['tertile1', insecure_accuracy1, i6_cutoffs[0][0.33], c24_cutoffs[0][0.33]])
		writefile.writerow(['tertile2', insecure_accuracy2, i6_cutoffs[0][0.66], c24_cutoffs[0][0.66]])
		writefile.writerow(['tertile3', insecure_accuracy3])

	#compare the full HFI scores to a score from fewer questions on the JC survey
	if i6_to_jc_cutoffs:
		#find intersection to represent those who logged responses for both
		intersection = np.intersect1d(cleani6, cleanjc)
		jcPD_T=jcPD.T
		i6PD_T=i6PD.T

		i6scores = np.zeros(intersection.shape)
		jcscores = np.zeros(intersection.shape)

		counter = -1
		for i in intersection:
			counter += 1
			for key in I6survey:
				if(key!='i6fsmortg'):
					if(key == "i6fssqmls"):
						i6scores[counter] -= i6PD_T[i][key]
					else:
					    i6scores[counter] += i6PD_T[i][key]
			# using the food security questions in JC survey
			#"jcfsrice"
			for key in ["jcfsrice", "jcfswor", "jcnofood", "jcslphngry", "jcdayhngry"]:
					jcscores[counter] += jcPD_T[i][key]

		scores_i6_pd = pd.DataFrame(data = i6scores[:])
		# scores_i6_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()
		scores_jc_pd = pd.DataFrame(data = jcscores[:])
		# scores_jc_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
		# plt.show()

		i6_secure_score = 3
		jc_secure_score = 0

		secure_i6 = np.where(scores_i6_pd[0] == i6_secure_score)[0]
		insecure_i6 =  np.where(scores_i6_pd[0] > i6_secure_score)[0]

		secure_jc = np.where(scores_jc_pd[0] == jc_secure_score)[0]
		insecure_jc =  np.where(scores_jc_pd[0] > jc_secure_score)[0]

		#determine number classified consistently as secure between i6 and jc
		intersection = np.intersect1d(secure_i6, secure_jc)
		secure = intersection.shape[0]/len(secure_i6)

		#determine number consistently as insecure between i6 and jc
		intersection = np.intersect1d(insecure_i6, insecure_jc)
		insecure = intersection.shape[0]/len(insecure_i6)

		insecure_i6_pd = scores_i6_pd[scores_i6_pd[0] > i6_secure_score]
		i6_cutoffs = insecure_i6_pd.quantile([0.33,0.66])

		insecure_jc_pd = scores_jc_pd[scores_jc_pd[0] > jc_secure_score]
		jc_cutoffs = insecure_jc_pd.quantile([0.33,0.66])

		#mild insecurity group, calculate consistency:
		i6_insecure = np.where(scores_i6_pd < i6_cutoffs[0][0.33])
		jc_insecure = np.where(scores_jc_pd <= jc_cutoffs[0][0.33])

		scores_i6_pd = insecure_i6_pd

		intersection = np.intersect1d(i6_insecure, jc_insecure)

		insecure_accuracy1 = intersection.shape[0]/len(i6_insecure[0])

		#moderate insecurity group, calculate consistency:
		i6_insecure = np.where((scores_i6_pd >= i6_cutoffs[0][0.33]) & (scores_i6_pd <= i6_cutoffs[0][0.66]))
		jc_insecure = np.where((scores_jc_pd > jc_cutoffs[0][0.33]) & (scores_jc_pd <= jc_cutoffs[0][0.66]))

		intersection = np.intersect1d(i6_insecure, jc_insecure)

		insecure_accuracy2 = intersection.shape[0]/len(i6_insecure[0])

		#severe insecurity group, calculate consistency:
		i6_insecure = np.where(scores_i6_pd > i6_cutoffs[0][0.66])
		jc_insecure = np.where(scores_jc_pd > jc_cutoffs[0][0.66])

		intersection = np.intersect1d(i6_insecure, jc_insecure)

		insecure_accuracy3 = intersection.shape[0]/len(i6_insecure[0])

		#write results to csv
		with open(wddata + 'i6_to_jc_comparison', 'w') as writefile:
			writefile  = csv.writer(writefile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			writefile.writerow(['category', 'percent', 'cutoff_i6', 'cutoff_jc'])
			writefile.writerow(['Food_secure', secure])
			writefile.writerow(['Food_insecure', insecure])
			writefile.writerow(['tertile1', insecure_accuracy1, i6_cutoffs[0][0.33], jc_cutoffs[0][0.33]])
			writefile.writerow(['tertile2', insecure_accuracy2, i6_cutoffs[0][0.66], jc_cutoffs[0][0.66]])
			writefile.writerow(['tertile3', insecure_accuracy3])
	#compare the full HFI scores to a score from fewer questions on the JC survey
	if i12_to_jc_cutoffs:
	  #find intersection to represent those who logged responses for both
	  intersection = np.intersect1d(cleani12, cleanjc)
	  jcPD_T=jcPD.T
	  i12PD_T=i12PD.T

	  i12scores = np.zeros(intersection.shape)
	  jcscores = np.zeros(intersection.shape)

	  counter = -1
	  for i in intersection:
	    counter += 1
	    for key in I12survey:
	      if(key!='i12fsmortg'):
	        if(key == "i12fssqmls"):
	          i12scores[counter] -= i12PD_T[i][key]
	        else:
	            i12scores[counter] += i12PD_T[i][key]
	    # using the food security questions in JC survey
	    #"jcfsrice"
	    for key in ["jcfsrice", "jcfswor", "jcnofood", "jcslphngry", "jcdayhngry"]:
	        jcscores[counter] += jcPD_T[i][key]

	  scores_i12_pd = pd.DataFrame(data = i12scores[:])
	  # scores_i12_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()
	  scores_jc_pd = pd.DataFrame(data = jcscores[:])
	  # scores_jc_pd.plot(kind='hist', bins=[-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
	  # plt.show()

	  i12_secure_score = 3
	  jc_secure_score = 0

	  secure_i12 = np.where(scores_i12_pd[0] == i12_secure_score)[0]
	  insecure_i12 =  np.where(scores_i12_pd[0] > i12_secure_score)[0]

	  secure_jc = np.where(scores_jc_pd[0] == jc_secure_score)[0]
	  insecure_jc =  np.where(scores_jc_pd[0] > jc_secure_score)[0]

	  #determine number classified consistently as secure between i12 and jc
	  intersection = np.intersect1d(secure_i12, secure_jc)
	  secure = intersection.shape[0]/len(secure_i12)

	  #determine number consistently as insecure between i12 and jc
	  intersection = np.intersect1d(insecure_i12, insecure_jc)
	  insecure = intersection.shape[0]/len(insecure_i12)

	  insecure_i12_pd = scores_i12_pd[scores_i12_pd[0] > i12_secure_score]
	  i12_cutoffs = insecure_i12_pd.quantile([0.33,0.66])

	  insecure_jc_pd = scores_jc_pd[scores_jc_pd[0] > jc_secure_score]
	  jc_cutoffs = insecure_jc_pd.quantile([0.33,0.66])

	  #mild insecurity group, calculate consistency:
	  i12_insecure = np.where(scores_i12_pd < i12_cutoffs[0][0.33])
	  jc_insecure = np.where(scores_jc_pd <= jc_cutoffs[0][0.33])

	  scores_i12_pd = insecure_i12_pd

	  intersection = np.intersect1d(i12_insecure, jc_insecure)

	  insecure_accuracy1 = intersection.shape[0]/len(i12_insecure[0])

	  #moderate insecurity group, calculate consistency:
	  i12_insecure = np.where((scores_i12_pd >= i12_cutoffs[0][0.33]) & (scores_i12_pd <= i12_cutoffs[0][0.66]))
	  jc_insecure = np.where((scores_jc_pd > jc_cutoffs[0][0.33]) & (scores_jc_pd <= jc_cutoffs[0][0.66]))

	  intersection = np.intersect1d(i12_insecure, jc_insecure)

	  insecure_accuracy2 = intersection.shape[0]/len(i12_insecure[0])

	  #severe insecurity group, calculate consistency:
	  i12_insecure = np.where(scores_i12_pd > i12_cutoffs[0][0.66])
	  jc_insecure = np.where(scores_jc_pd > jc_cutoffs[0][0.66])

	  intersection = np.intersect1d(i12_insecure, jc_insecure)

	  insecure_accuracy3 = intersection.shape[0]/len(i12_insecure[0])

	  #write results to csv
	  with open(wddata + 'i12_to_jc_comparison', 'w') as writefile:
	    writefile  = csv.writer(writefile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	    writefile.writerow(['category', 'percent', 'cutoff_i12', 'cutoff_jc'])
	    writefile.writerow(['Food_secure', secure])
	    writefile.writerow(['Food_insecure', insecure])
	    writefile.writerow(['tertile1', insecure_accuracy1, i12_cutoffs[0][0.33], jc_cutoffs[0][0.33]])
	    writefile.writerow(['tertile2', insecure_accuracy2, i12_cutoffs[0][0.66], jc_cutoffs[0][0.66]])
	    writefile.writerow(['tertile3', insecure_accuracy3])
