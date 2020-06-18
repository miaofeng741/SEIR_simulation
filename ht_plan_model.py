
# -*- coding: utf-8 -*-
# Before running this model, you need to store the initial value in the system variable (from ht_plan_model_initial.py)
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy 
import scipy.stats as st
import scipy.special as sps
from scipy.stats import gamma
from scipy.stats import lognorm
from scipy.stats import norm
from scipy.stats import uniform
from collections import Counter
import datetime

def one_back_to_past_forE(distribution, a, b, N, back_len):#Random number generator
###distribution: Distribution type，a\b：Distributed parameter，N:Number of people
###back_len:Length of time

    ### uniform
    N = int(N)
    cdf_0_1 = numpy.random.uniform(0.0001,distribution.cdf(back_len, a, scale = b),N) 
    #print('Random 0-1 CDF generated')
    #print(cdf_0_1)
    #############################Inverse method##################################
    random_N = []
    for i in range(len(cdf_0_1)):
        random = distribution.ppf(cdf_0_1[i], a, scale = b)
        random_N.append(math.ceil(random))
    #print(random_N)
    #########################list#############
    back_all = [0 for x in range(1, (back_len+1))]
    for i in range(len(random_N)):
        back_list = [0 for m in range(1,(back_len+1))]###produce back list
        back_list[int(random_N[i])-1] = 1 
        if i ==0:
            #Initialize
            back_all = back_list
        else:
            for n in range(len(back_list)):
                back_all[n] = back_all[n] + back_list[n]
    # print('Random CDF\'s x')
    # print(random_N)
    return back_all

def one_back_to_past(distribution, mu, sigma, N, begin,back_len):#Random number generator
###distribution: Distribution type，a\b：Distributed parameter，N:Number of people
###back_len:Length of time
    
    ##
    ##
    N = int(N)
    cdf_0_1 = numpy.random.uniform(distribution.cdf(begin-1+0.0001, mu, scale =sigma),distribution.cdf(back_len, mu, scale =sigma),N) 
    #print('Random 0-1 CDF generated')
    #print(cdf_0_1)
    #############################Inverse method##################################
    random_N = []
    for i in range(len(cdf_0_1)):
        random = distribution.ppf(cdf_0_1[i], mu, scale =sigma)
        random_N.append(math.ceil(random))
    #print(random_N)
    #########################list#############
    back_all = [0 for x in range(1, (back_len+1))]
    for i in range(len(random_N)):
        back_list = [0 for m in range(1,(back_len+1))]
        back_list[int(random_N[i])-1] = 1 
        if i ==0:
            #Initialize
            back_all = back_list
        else:
            for n in range(len(back_list)):
                back_all[n] = back_all[n] + back_list[n]
    # print('Random CDF\'s x')
    # print(random_N)
    return back_all


def DQ_to_DQ_wait(DQ_new):#When can a person waiting to be treated heal himself
    back_data = [i for i in range(1,29)]#title
    back_list = [] 
    back_list = one_back_to_past(gamma,22.563, 0.842, DQ_new, 13,28) 
    back_data = numpy.row_stack((back_data,back_list))
    back_data= back_data[1,:] 
    return back_data

def SDQ_New_wait(SQ_New,DQ_wait): # Simulate the process of receiving patients in hospital, those who get sick are admitted first
    N = int(SQ_New) 
    DQ = sum(DQ_wait)
    if N > DQ:
        rest = N - DQ
        N = DQ 
    else:
        rest = 0 
    #N patients were admitted
    t = 0 #count
    back_all = [0 for x in range(len(DQ_wait))]
    while True:
        if t < N: 
            X = math.ceil(numpy.random.uniform(0,28,1)) 
            back_list = [0 for m in range(len(DQ_wait))]
            back_list[X-1] = 1 
            if (back_list[X-1]+back_all[X-1]) <=DQ_wait[X-1]: 
                for n in range(len(back_list)): 
                        back_all[n] = back_all[n] + back_list[n]
                t = t+1  
            if t >= N:           
                break
        else:
            break
    #update DQ_wait
    for t in range(len(DQ_wait)):
        DQ_wait[t] = DQ_wait[t] - back_all[t]
    return DQ_wait,rest


def Rule_new_FF(DQ_new_wait_max,SQ_new): #Update the queue for admission
    SQ_new_remained = SQ_new #rest number 
    for t in range(1,DQ_new_wait_max.shape[0]):
        if SQ_new_remained>0:
            DQ_new_wait_max[t,:],SQ_new_remained = SDQ_New_wait(SQ_new_remained,DQ_new_wait_max[t,:]) 
        else:
            break
    return DQ_new_wait_max

def SDZ_New_wait(SZ_new,DZ_wait):
    SZ_new_wait = [0]*len(DZ_wait)
    t = 0 
    N = 0 
    for t in range(len(DZ_wait)):
        N = N + DZ_wait[t]
        if N <= SZ_new:
            SZ_new_wait[t] = DZ_wait[t]
        else:
            SZ_new_wait[t] = SZ_new-(N-DZ_wait[t])
            break
    return SZ_new_wait

def E_to_E_wait(E_new):  #The state of the patient during the incubation period
    ###read lognorm  distribution###############
    back_data = [i for i in range(1,22)]#title
    back_list = [] 
    back_list = one_back_to_past_forE(lognorm,0.640037, 4.2369, E_new,21)
    back_data = numpy.row_stack((back_data,back_list)) 
    back_data= back_data[1:,:] 
    return back_data
########
##########

def E_to_E_wait_times(E_new,times):
    E_new_wait_max = [i for i in range(21)]
    E_new_wait_mean = [i for i in range(21)]
    for i in range(times):
        E_new_wait= E_to_E_wait(E_new)
        E_new_wait_max = numpy.row_stack((E_new_wait_max,E_new_wait))
    E_new_wait_max = E_new_wait_max[1:,:]
    for i in range(len(E_new_wait_mean)):
        E_new_wait_mean[i] = numpy.mean(E_new_wait_max[:,i])
    return E_new_wait_mean

def Et_out(Et,Eout): #Calculate the incubation period of the city patient status
    N = int(Eout) 
    #start 
    t = 0 
    back_all = [0 for x in range(len(Et))]
    while True:
        if t < N: 
            X = math.ceil(numpy.random.uniform(0,21,1))  
            back_list = [0 for m in range(len(Et))]
            back_list[X-1] = 1
            if (back_list[X-1]+back_all[X-1]) <=Et[X-1]:
                for n in range(len(back_list)):
                        back_all[n] = back_all[n] + back_list[n]
                t = t+1  
            if t >= N:           
                break
        else:
            break
    #up date DQ_wait 
    for t in range(len(Et)):
        Et[t] = Et[t] - back_all[t]
    return Et



def SQ_to_SQ_wait(SQ_new):# Time to cure mild patients
    back_data = [i for i in range(1,27)]#title
    back_list = [] 
    back_list = one_back_to_past(gamma,18.062,0.941, SQ_new,11,26) 
    back_data = numpy.row_stack((back_data,back_list)) 
    back_data= back_data[1,:]  
    return back_data

def SQ_to_SQ_wait_times(SQ_new,times): # Time to cure mild patients_MTKL
    SQ_new_wait_max = [i for i in range(26)]
    SQ_new_wait_mean = [i for i in range(26)]
    for i in range(times):
        SQ_new_wait= SQ_to_SQ_wait(SQ_new)
        SQ_new_wait_max = numpy.row_stack((SQ_new_wait_max,SQ_new_wait))
    SQ_new_wait_max = SQ_new_wait_max[1:,:]
    for i in range(len(SQ_new_wait_mean)):
        SQ_new_wait_mean[i] = numpy.mean(SQ_new_wait_max[:,i])
    SQ_new_wait_mean=numpy.rint(SQ_new_wait_mean)
    return SQ_new_wait_mean

def SZ_to_SZ_wait(SZ_new):#Time to cure critically ill patients
    back_data = [i for i in range(1,42)]
    back_list = [] 
    back_list = one_back_to_past(gamma,14.951, 1.551, SZ_new, 12,41) 
    back_data = numpy.row_stack((back_data,back_list)) 
    back_data= back_data[1,:] 
    return back_data

def SZ_to_SZ_wait_times(SZ_new,times):
    SZ_new_wait_max = [i for i in range(41)]
    SZ_new_wait_mean = [i for i in range(41)]
    for i in range(times):
        SZ_new_wait= SZ_to_SZ_wait(SZ_new)
        SZ_new_wait_max = numpy.row_stack((SZ_new_wait_max,SZ_new_wait))
    SZ_new_wait_max = SZ_new_wait_max[1:,:]
    for i in range(len(SZ_new_wait_mean)):
        SZ_new_wait_mean[i] = numpy.mean(SZ_new_wait_max[:,i])
    SZ_new_wait_mean=numpy.rint(SZ_new_wait_mean)
    return SZ_new_wait_mean

def DQ_to_DQ_wait_times(DQ_new,times):
    DQ_new_wait_max = [i for i in range(28)]
    DQ_new_wait_mean = [i for i in range(28)]
    for i in range(times):
        DQ_new_wait= DQ_to_DQ_wait(DQ_new)
        DQ_new_wait_max = numpy.row_stack((DQ_new_wait_max,DQ_new_wait))
    DQ_new_wait_max = DQ_new_wait_max[1:,:]
    for i in range(len(DQ_new_wait_mean)):
        DQ_new_wait_mean[i] = numpy.mean(DQ_new_wait_max[:,i])
    DQ_new_wait_mean=numpy.rint(DQ_new_wait_mean)
    return DQ_new_wait_mean


if __name__ == '__main__':
    bed_list = numpy.genfromtxt("different_ht_plan.csv", 
                          dtype=str,  delimiter=",", skip_header = 1)
    BQ_back_list_all =  numpy.array(bed_list[:,3], dtype=numpy.float32) 
    BZ_back_list_all =  numpy.array(bed_list[:,4], dtype=numpy.float32) 
    BZ_back_list = BZ_back_list_all[20:]
    BQ_back_list = BQ_back_list_all[20:]
    αQ=0.81
    αZ=0.19

    times = 100
    for t in range(25): #There will be 25 days from December 27 to January 20   It could be different lengths depending on what you want, like 21, 26
       print(t)
   
       E_new = round((S/N)*beta_list[t]*(Et[0]+Et[1]+I_list[20+t]))  
       E_new_t.append(E_new) 
       if E_new < 0:
           E_new = 0
       print(t)
       I_new_temp_list = [] 
       Et_temp_max = [i for i in range(len(Et))] 
       S_temp_list =[]
       N_temp_list =[]
       DQ_temp_list = [] 
       RQ_new_temp_list =[]
       DQ_new_wait_max_temp_all = 0 
       SQ_new_temp_list = []
       DZ_temp_list = []
       RZ_new_temp_list =[]
       DZ_wait_temp_max = [i for i in range(len(DZ_wait))] 
       SZ_new_temp_list =[]
       SQD_new_temp_list =[] 
       ht_temp_list =[] 
       I_temp_list =[] 
       SQ_wait_temp_max = [i for i in range(len(SQ_wait))] 
       SZ_wait_temp_max = [i for i in range(len(SZ_wait))] 
       SQD_wait_temp_max = [i for i in range(len(SQD_wait))] 
       
       LQ_new_temp_list =[]
       LZ_new_temp_list =[]
       LQD_new_temp_list =[]
       BQ_temp_list = []
       BZ_temp_list = []
       for x in range(times):# I'm going to loop around times to get the mean of the states
           ###############
           Et_temp = Et*1
           I_new_list_temp = I_new_list[:]
           S_temp = S
           N_temp = N
           DQ_temp=DQ
           RQ_new_temp = RQ_new 
           BQ_new_list_temp = BQ_new_list[:]
           DQ_new_wait_max_temp = DQ_new_wait_max*1
           RZ_new_temp = RZ_new
           DZ_temp = DZ
           BQ_temp = BQ
           BZ_temp = BZ
           
           
           BZ_new_list_temp =BZ_new_list[:]
           DZ_wait_temp = DZ_wait[:]
           SQ_new_list_temp = SQ_new_list[:]
           SZ_new_list_temp = SZ_new_list[:]
           SQD_new_list_temp = SQD_new_list[:]
           RQ_new_list_temp = RQ_new_list[:]
           RZ_new_list_temp = RZ_new_list[:]
           I_temp = I
           SQ_wait_temp = SQ_wait[:]
           SZ_wait_temp = SZ_wait[:]
           SQD_wait_temp = SQD_wait[:]
           
           LQ_new_list_temp = LQ_new_list[:]
           LQD_new_list_temp = LQD_new_list[:]
           LZ_new_list_temp= LZ_new_list[:]
           
           #Step 3 use E_new to get E_new_wait
           E_new_wait_temp = E_to_E_wait_times(E_new,1) 
           E_new_wait_temp=numpy.rint(E_new_wait_temp)  
           
           #E_new_wait_MAX = numpy.row_stack((E_new_wait_MAX,E_new_wait)) 
       
           
           I_new_temp=Et_temp[0]+E_new_wait_temp[0]
           I_new_list_temp.append(I_new_temp)
           
           
           for i in range(20):
               Et_temp[i]=Et_temp[i+1]+E_new_wait_temp[i+1]######前20个
           Et_temp[20]=0
           
           #Consider population mobility
           #only consider E and S 
           E_temp = sum(Et_temp) 
           S_temp = S_temp- E_new 
           
           Eout_temp = round((E_temp/N_temp)*population_out[t])
           Sout_temp = round((S_temp/N_temp)*population_out[t])
           
           
           Et_temp = Et_out(Et_temp,Eout_temp)
           
           S_temp=S_temp+population_in[t]-E_new-Sout_temp
           N_temp=N_temp+population_in[t]-population_out[t]
           
           
           #######ht###########
           DQ_new_temp = round(αQ * I_new_list_temp[t]) 
           DQ_temp = DQ_temp + DQ_new_temp - RQ_new_temp 
           BQ_temp = BQ_new_list_temp[t]   
           DQ_new_wait_temp= DQ_to_DQ_wait_times(DQ_new_temp,1) 
       
           DQ_new_wait_max_temp = numpy.row_stack((DQ_new_wait_max_temp,DQ_new_wait_temp))
           SQ_new_temp = min(BQ_temp,DQ_temp) 
           DQ_new_wait_max_temp = Rule_new_FF(DQ_new_wait_max_temp,SQ_new_temp) 
           DQ_temp = DQ_temp - SQ_new_temp
           BQ_temp = BQ_temp - SQ_new_temp
           DZ_new_temp = round(αZ * I_new_list_temp[t]) 
           DZ_temp = DZ_temp+  DZ_new_temp - RZ_new_temp 
           BZ_temp = BZ_new_list_temp[t]   
       
           DZ_wait_temp[13] = DZ_wait_temp[13] + DZ_new_temp
           SZ_new_temp = min(DZ_temp,BZ_temp) 
           SDZ_new_wait_temp = SDZ_New_wait(SZ_new_temp,DZ_wait_temp)
           DZ_wait_temp=[DZ_wait_temp[i]-SDZ_new_wait_temp[i] for i in range(len(DZ_wait_temp))]

           DZ_temp = DZ_temp - SZ_new_temp
           BZ_temp = BZ_temp - SZ_new_temp


           if BZ_temp*DQ_temp > 0:
               SQD_new_temp = min (BZ_temp,DQ_temp)
           else:
               SQD_new_temp = 0
           
          
           DQ_new_wait_max_temp = Rule_new_FF(DQ_new_wait_max_temp,SQD_new_temp)

           DQ_temp = DQ_temp - SQD_new_temp
           BZ_temp = BZ_temp - SQD_new_temp

           ht_temp= SQ_new_temp + SZ_new_temp + SQD_new_temp 
           ht_temp_list.append(ht_temp) 

           RQ_new_temp = sum(DQ_new_wait_max_temp[1:,0]) 
           DQ_new_wait_max_temp = DQ_new_wait_max_temp[:,1:]
   
           back_0 = numpy.ones((DQ_new_wait_max_temp.shape[0],1))*0
           DQ_new_wait_max_temp = numpy.column_stack((DQ_new_wait_max_temp,back_0))

           RZ_new_temp = DZ_wait_temp[0]
           del DZ_wait_temp[0]
           DZ_wait_temp = DZ_wait_temp +[0]
       
       
           SQ_new_list_temp.append(SQ_new_temp) 
           SZ_new_list_temp.append(SZ_new_temp) 
           SQD_new_list_temp.append(SQD_new_temp)
           RQ_new_list_temp.append(RQ_new_temp) 
           RZ_new_list_temp.append(RZ_new_temp) 
           
               

           SQ_new_wait_temp= SQ_to_SQ_wait_times(SQ_new_temp,1)
           SQ_wait_temp=[SQ_wait_temp[i]+SQ_new_wait_temp[i] for i in range(len(SQ_wait_temp))]
           #SQ_wait_max_1 = numpy.row_stack((SQ_wait_max_1,SQ_wait))  
           SZ_new_wait_temp= SZ_to_SZ_wait_times(SZ_new_temp,1) 
           SZ_wait_temp=[SZ_wait_temp[i]+SZ_new_wait_temp[i] for i in range(len(SZ_wait_temp))]


           SQD_new_wait_temp = SQ_to_SQ_wait_times(SQD_new_temp,1) 
           SQD_wait_temp=[SQD_wait_temp[i]+SQD_new_wait_temp[i] for i in range(len(SQD_wait_temp))]
 

           LQ_new_temp = SQ_wait_temp[0]
           del SQ_wait_temp[0]
           SQ_wait_temp = SQ_wait_temp +[0]
           LQD_new_temp = SQD_wait_temp[0]
           del SQD_wait_temp[0]
           SQD_wait_temp = SQD_wait_temp +[0]
           LZ_new_temp = SZ_wait_temp[0]
           del SZ_wait_temp[0]
           SZ_wait_temp = SZ_wait_temp +[0]
           
           LQ_new_list_temp.append(LQ_new_temp)
           LQD_new_list_temp.append(LQD_new_temp)
           LZ_new_list_temp.append(LZ_new_temp)
           BQ_temp_list.append(BQ_temp)
           BZ_temp_list.append(BZ_temp)
           #save data
           I_new_temp_list.append(I_new_temp)  
           Et_temp_max = numpy.row_stack((Et_temp_max,Et_temp)) 
           S_temp_list.append(S_temp)
           N_temp_list.append(N_temp)
           DQ_temp_list.append(DQ_temp)
           RQ_new_temp_list.append(RQ_new_temp)
           DQ_new_wait_max_temp_all =  DQ_new_wait_max_temp_all  + DQ_new_wait_max_temp
           DZ_temp_list.append(DZ_temp)
           RZ_new_temp_list.append(RZ_new_temp)
           DZ_wait_temp_max = numpy.row_stack((DZ_wait_temp_max,DZ_wait_temp))
           SQ_new_temp_list.append(SQ_new_temp) 
           SZ_new_temp_list.append(SZ_new_temp)
           SQD_new_temp_list.append(SQD_new_temp)
           I_temp_list.append(I_temp)
           SQ_wait_temp_max = numpy.row_stack((SQ_wait_temp_max,SQ_wait_temp))
           SZ_wait_temp_max = numpy.row_stack((SZ_wait_temp_max,SZ_wait_temp))
           SQD_wait_temp_max = numpy.row_stack((SQD_wait_temp_max,SQD_wait_temp))
           
           LQ_new_temp_list.append(LQ_new_temp)
           LZ_new_temp_list.append(LZ_new_temp)
           LQD_new_temp_list.append(LQD_new_temp)

       I_new = round(numpy.mean(I_new_temp_list))
       I_new_list.append(I_new)
       Et_temp_max = Et_temp_max[1:] 
       Et = [i for i in range(len(Et_temp))]  
       for i in range(len(Et_temp)):
           Et[i] = round(numpy.mean(Et_temp_max[:,i]))
       
       S= numpy.mean(S_temp_list) 
       N= numpy.mean(N_temp_list)  
       RQ_new_list_100 = numpy.row_stack((RQ_new_list_100,RQ_new_temp_list)) 
       RQ_new =  round(numpy.mean(RQ_new_temp_list))
       RQ_new_list.append(RQ_new)
       DQ_new_wait_max = DQ_new_wait_max_temp_all/times
       DQ_new_wait_max = numpy.rint(DQ_new_wait_max)
       DQ_sum  = DQ_new_wait_max[1:]
       DQ = DQ_sum.sum()
       DQ_list.append(DQ)
       I_new_list_100 = numpy.row_stack((I_new_list_100,I_new_temp_list))
       DQ_list_100 = numpy.row_stack((DQ_list_100,DQ_temp_list)) 
       DZ_list_100 = numpy.row_stack((DZ_list_100,DZ_temp_list)) 
       RZ_new_list_100= numpy.row_stack((RZ_new_list_100,RZ_new_temp_list)) 
       RZ_new =  round(numpy.mean(RZ_new_temp_list)) 
       RZ_new_list.append(RZ_new) 
       DZ_wait_temp_max = DZ_wait_temp_max[1:] 
       DZ_wait = [i for i in range(len(DZ_wait_temp))] 
       for i in range(len(DZ_wait)):
           DZ_wait[i] = round(numpy.mean(DZ_wait_temp_max[:,i]))
       
       DZ = 0 
       for ele in range(0, len(DZ_wait)): 
           DZ = DZ + DZ_wait[ele] 
       DZ_list.append(DZ)    
           
       SQ_new = round(numpy.mean(SQ_new_temp_list))
       SQ_new_list.append(SQ_new)
       SQ_new_list_100 = numpy.row_stack((SQ_new_list_100,SQ_new_temp_list)) 
       SZ_new =  round(numpy.mean(SZ_new_temp_list)) 
       SZ_new_list.append(SZ_new)
       SZ_new_list_100= numpy.row_stack((SZ_new_list_100,SZ_new_temp_list))
       
       SQD_new =  round(numpy.mean(SQD_new_temp_list)) 
       SQD_new_list.append(SQD_new)
       SQD_new_list_100= numpy.row_stack((SQD_new_list_100,SQD_new_temp_list))
       
       ht =  round(numpy.mean(ht_temp_list)) 
       ht_list.append(ht)
       ht_list_100= numpy.row_stack((ht_list_100,ht_temp_list))
       
       if t <= 25:
               I = I +  I_new_list[t] - ht_list[t+20] - RQ_new_list[t+20] - RZ_new_list[t+20]
       else:
               I = DQ+DZ     
       I_list.append(I)
       
       I_list_100= numpy.row_stack((I_list_100,I_temp_list))     
           
       SQ_wait_temp_max = SQ_wait_temp_max[1:] 
       SQ_wait = [i for i in range(len(SQ_wait_temp))]  
       for i in range(len(SQ_wait)):
           SQ_wait[i] = round(numpy.mean(SQ_wait_temp_max[:,i]))
           
       SZ_wait_temp_max = SZ_wait_temp_max[1:] 
       SZ_wait = [i for i in range(len(SZ_wait_temp))]  
       for i in range(len(SZ_wait)):
           SZ_wait[i] = round(numpy.mean(SZ_wait_temp_max[:,i]))
       
       SQD_wait_temp_max = SQD_wait_temp_max[1:] 
       SQD_wait = [i for i in range(len(SQD_wait_temp))]  
       for i in range(len(SQD_wait)):
           SQD_wait[i] = round(numpy.mean(SQD_wait_temp_max[:,i])  )    
               
       LQ_new_list_100= numpy.row_stack((LQ_new_list_100,LQ_new_temp_list)) 
       LQ_new =  round(numpy.mean(LQ_new_temp_list)) 
       LQ_new_list.append(LQ_new) 
       
       LZ_new_list_100= numpy.row_stack((LZ_new_list_100,LZ_new_temp_list)) 
       LZ_new =  round(numpy.mean(LZ_new_temp_list)) 
       LZ_new_list.append(LZ_new) 
       
       
       LQD_new_list_100= numpy.row_stack((LQD_new_list_100,LQD_new_temp_list)) 
       LQD_new =  round(numpy.mean(LQD_new_temp_list)) 
       LQD_new_list.append(LQD_new) 
       BQ =  round(numpy.mean(BQ_temp_list)) 
       BZ =  round(numpy.mean(BZ_temp_list)) 