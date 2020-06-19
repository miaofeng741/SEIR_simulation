# -*- coding: utf-8 -*-
#This procedure is used to simulate hospital admission and disease 
#development under the condition of known transmission rate,
# and calculate the value of each state variable at key time nodes 
#such as 1.20 1.23.1.25

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
    back_list = one_back_to_past(gamma,18.063,0.941, SQ_new,11,26) 
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
    back_list = one_back_to_past(gamma,14.951, 1.552, SZ_new, 12,41) 
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
    ### Read known data### 
    bed_number_data = numpy.genfromtxt("bed_number.csv",  dtype = str , 
                                delimiter=",", skip_header = 1) #Read known Inewt data and bed data
    WU_I_new_list_all = bed_number_data[:,1] #The first column is the Inewt data known in the paper
    BZ_new_list_all = bed_number_data[:,4] # The fourth column is the number of newly admitted patients in designated hospitals (admitted for serious diseases) every day
    BQ_new_list_all = bed_number_data[:,5] # The fifth column is the number of new patients admitted to the makeshift hospitals (for mild cases) every day
    BZ_back_list_all = bed_number_data[:,2] # The second column is the number of new beds per day in designated hospitals
    BQ_back_list_all = bed_number_data[:,3] # The third column is the number of new beds per day in fangcang hospitals
    #Data type conversion
    WU_I_new_list_all = numpy.array(WU_I_new_list_all, dtype=numpy.float32)
    BZ_new_list_all = numpy.array(BZ_new_list_all, dtype=numpy.float32)
    BQ_new_list_all = numpy.array(BQ_new_list_all, dtype=numpy.float32)
    BZ_back_list_all = numpy.array(BZ_back_list_all, dtype=numpy.float32)
    BQ_back_list_all = numpy.array(BQ_back_list_all, dtype=numpy.float32)
    ##########Reads the value of a lognorma distribution#############
    lognormal_data = numpy.genfromtxt("Incubation_period_lognormal_list.csv", 
                          dtype=str,  delimiter=",", skip_header = 1)
    lognormal =  numpy.array(lognormal_data[:,0], dtype=numpy.float32)
    ##########Read the floating population data################
    population_data = numpy.genfromtxt("movement of population.csv", 
                          dtype=str,  delimiter=",", skip_header = 1)
    population_in_all = population_data[:,1]  
    population_out_all = population_data[:,2] 
    population_in_all =  numpy.array(population_in_all, dtype=numpy.float32)
    population_out_all =  numpy.array(population_out_all, dtype=numpy.float32)
    ##########Read beta data
    beta_list = numpy.genfromtxt("beta_list.csv", 
                          dtype=str,  delimiter=",", skip_header = 1)
    #plan1
    beta_list =  numpy.array(beta_list[:,2], dtype=numpy.float32)
    
    ###########Read initial Et#####
    Et_random_data = numpy.genfromtxt("Et_random.csv", 
                          dtype=str,  delimiter=",", skip_header = 1)
    Et_random =  numpy.array(Et_random_data[:,1:], dtype=numpy.float32)
    n= int(numpy.random.uniform(0,50-0.00001,1))
    Et = Et_random[0]
    
    
    #Store the initial value
    times = 1000#Number of experimental cycles
    αQ = 0.81
    αZ = 0.19 #Proportion of severe cases in Inewt
    DQ = 0 #The number of patients to be treated at the initial time is 0
    RQ_new = 0 #The number of mild cases self-healing at the initial time was 0
    DZ = 0 
    RZ_new = 0 
    DZ_wait = [0]*14 
    I = 0 
    ht = 0
    N = 11212000 
    S = 11211861  
    SQ_new_list= [] #Save the data
    SZ_new_list= []  #Save the data
    SQD_new_list= [] #Save the data
    RQ_new_list= [0]#Save the data
    RZ_new_list= [0] #Save the data
    LQ_new_list = [0]#Save the data
    LQD_new_list = [0]#Save the data
    LZ_new_list = [0]#Save the data
    ht_list = [] #Save the data
    I_list = [0] 
    E_new_t = [] 
    record_Eout = []
    record_Sout = []
    record_S = [11211861] 
    record_N = [11212000]
    E= 0
    E_new_wait_MAX=[i for i in range(1,22)]
    DQ_new_wait_max = [i for i in range(28) ]
    Et_Max=[i for i in range(1,22)] 
    SQ_wait_max_1 = [[i for i in range(1,27)]] 
    SZ_wait_max = [[i for i in range(1,42)]]
    SQD_wait_max =[[i for i in range(1,27)]]
    SQ_wait = [0]*26 
    SZ_wait = [0]*41 
    SQD_wait = [0]*26 
    A_list = [] 
    DQ_list = []
    DZ_list = []
    
    for t in range(20): #From December 7th to December 26th, we still use the real Inewt
        print(t)
        #Mild part step1-1：The determination of the initial state
        DQ_new = round(αQ * WU_I_new_list_all[t])
        DQ = DQ + DQ_new - RQ_new
        BQ = BQ_new_list_all[t] 
        #Mild part step1-2：Update of Mild to be treated
        DQ_new_wait= DQ_to_DQ_wait_times(DQ_new,times)
        #DQ_new_wait= DQ_to_DQ_wait(DQ_new)
        DQ_new_wait_max = numpy.row_stack((DQ_new_wait_max,DQ_new_wait))
        DQ_sum  = DQ_new_wait_max[1:]
        DQ = DQ_sum.sum()
        
        
        #Mild part step1-3:To calculate the number of new mild cases admitted
        SQ_new = min(BQ,DQ) 
        DQ_new_wait_max = Rule_new_FF(DQ_new_wait_max,SQ_new)
        #Mild part step1-4 update DQ
        DQ = DQ - SQ_new
        BQ = BQ - SQ_new
        #Critical part step2-1：The determination of the initial state
        DZ_new = round(αZ * WU_I_new_list_all[t]) 
        DZ = DZ +  DZ_new - RZ_new 
        BZ = BZ_new_list_all[t]   
        
        #Critical part step2-2：    Update critically ill patients awaiting admission
        DZ_wait[13] = DZ_wait[13] + DZ_new
        
        #Critical part step2-3:To calculate the number of newly hospitalized patients with severe diseases
        SZ_new = min(DZ,BZ)
        SDZ_new_wait = SDZ_New_wait(SZ_new,DZ_wait)
        DZ_wait=[DZ_wait[i]-SDZ_new_wait[i] for i in range(len(DZ_wait))]
        #Critical partstep2-4 update DZ
        DZ = DZ - SZ_new
        BZ = BZ - SZ_new
        #Special mild part step3
        if BZ*DQ > 0:       
            SQD_new = min (BZ,DQ)
        else:
            SQD_new = 0
        
        DQ_new_wait_max = Rule_new_FF(DQ_new_wait_max,SQD_new)
        #step3-2：Update patients to be admitted
        DQ = DQ - SQD_new
        BZ = BZ - SQD_new

        ht = SQ_new + SZ_new + SQD_new 
        ht_list.append(ht) 
        #DQ_wait State transition
        RQ_new = sum(DQ_new_wait_max[1:,0]) 
        DQ_new_wait_max = DQ_new_wait_max[:,1:] 
    
        back_0 = numpy.ones((DQ_new_wait_max.shape[0],1))*0
        DQ_new_wait_max = numpy.column_stack((DQ_new_wait_max,back_0))
        #DZ_wait State transition
        RZ_new = DZ_wait[0]
        del DZ_wait[0]
        DZ_wait = DZ_wait +[0]
        
        
        SQ_new_list.append(SQ_new) 
        SZ_new_list.append(SZ_new) 
        SQD_new_list.append(SQD_new) 
        RQ_new_list.append(RQ_new) 
        RZ_new_list.append(RZ_new) 
        
        I = I +   WU_I_new_list_all[t] - ht_list[t] - RQ_new_list[t] - RZ_new_list[t]
        
        I_list.append(I)
        
    
        #A mild case was admitted to the fangcang hospital
        SQ_new_wait= SQ_to_SQ_wait_times(SQ_new,times) 
        SQ_wait=[SQ_wait[i]+SQ_new_wait[i] for i in range(len(SQ_wait))]
        #SQ_wait_max_1 = numpy.row_stack((SQ_wait_max_1,SQ_wait)) 
        #Severe cases were admitted to designated hospitals
        SZ_new_wait= SZ_to_SZ_wait_times(SZ_new,times) 
        SZ_wait=[SZ_wait[i]+SZ_new_wait[i] for i in range(len(SZ_wait))]
        #SZ_wait_max = numpy.row_stack((SZ_wait_max ,SZ_wait))
        #Designated hospitals for mild cases
        SQD_new_wait = SQ_to_SQ_wait_times(SQD_new,times) 
        SQD_wait=[SQD_wait[i]+SQD_new_wait[i] for i in range(len(SQ_wait))]
        SQD_wait_max = numpy.row_stack((SQD_wait_max,SQD_wait)) 
        
       
        #SQ_waitState transition
        LQ_new = SQ_wait[0]
        del SQ_wait[0]
        SQ_wait = SQ_wait +[0]
          #SQD_wait State transition
        LQD_new = SQD_wait[0]
        del SQD_wait[0]
        SQD_wait = SQD_wait +[0]
        #SZ_wait State transition
        LZ_new = SZ_wait[0]
        del SZ_wait[0]
        SZ_wait = SZ_wait +[0]
        
        #save data
        LQ_new_list.append(LQ_new)
        LQD_new_list.append(LQD_new)
        LZ_new_list.append(LZ_new)
        
        DQ_list.append(DQ)
        DZ_list.append(DZ)


    #Update the bed information#
    BQ_new_list = BQ_new_list_all[20:]
    BZ_new_list = BZ_new_list_all[20:]
    I_new_list = [28] 
    population_in = population_in_all[20:]
    population_out = population_out_all[20:]
    BZ_back_list = BZ_back_list_all[20:]
    BQ_back_list = BQ_back_list_all[20:]
    ######Let's start with updates from December 27 onwards
    αQ=0.81
    αZ=0.19
    times = 100  
    
    #########The variable _temp that needs to be stored
    I_new_list_100 = [i for i in range(times)]
    DQ_list_100 = [i for i in range(times)]
    
    RQ_new_list_100 = [i for i in range(times)] 
    SQ_new_list_100 = [i for i in range(times)]
    
    DZ_list_100 = [i for i in range(times)] 
    RZ_new_list_100 = [i for i in range(times)] 
    SZ_new_list_100 = [i for i in range(times)] 
    SQD_new_list_100 = [i for i in range(times)] 
    ht_list_100 = [i for i in range(times)] 
    I_list_100 = [i for i in range(times)] 
    
    LQ_new_list_100 = [i for i in range(times)] 
    LZ_new_list_100 = [i for i in range(times)] 
    LQD_new_list_100 = [i for i in range(times)] 
    

    
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
