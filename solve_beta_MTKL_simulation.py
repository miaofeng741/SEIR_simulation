# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
# This program is used to realize the solution of beta, using monte Carlo simulation principle
# It is worth noting that when using different initial values, you need to change the program at line 456 to give 10 different initial values for the experiment
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

def Rule_new_FF(DQ_new_wait_max,SQ_new): #Update the queue for admission
    SQ_new_remained = SQ_new #rest number 
    for t in range(1,DQ_new_wait_max.shape[0]):
        if SQ_new_remained>0:
            DQ_new_wait_max[t,:],SQ_new_remained = SDQ_New_wait(SQ_new_remained,DQ_new_wait_max[t,:]) 
        else:
            break
    return DQ_new_wait_max

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

def SZ_to_SZ_wait_times(SZ_new,times): #Time to cure critically ill patients_MTKL
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

def ht_model(WU_I_new_list_all,BZ_new_list_all,BQ_new_list_all):
    WU_I_new_list = WU_I_new_list_all
    BZ_new_list = BZ_new_list_all
    BQ_new_list = BQ_new_list_all

    #Define constants
    aerfaQ = 0.81 #Proportion of patients with mild symptoms
    aerfaZ = 0.19 #Proportion of critically ill patients
    
    #Assign an initial value to a state variable
    DQ_new = 0 #New initial value for patients with mild disease
    DQ = 0 #The total number of new mild cases
    RQ_new = 0
    BQ = 0 
    SQ_new = 0
    
    DZ_new = 0 #The initial value of the new critically ill patients
    DZ = 0 #The total number of new mild cases
    RZ_new = 0 
    BZ = 0
    
    #Initial state matrix assignment
    
    DZ_wait = [0]*14 #Critically ill patients awaiting admission
    
    SQ_wait = [0]*26 
    SZ_wait = [0]*41 
    SQD_wait = [0]*26 
    SQ_wait_max_1 = [[i for i in range(1,27)]] 
    SZ_wait_max = [[i for i in range(1,42)]] 
    SQD_wait_max =[[i for i in range(1,27)]]
    DQ_new_wait_max = [i for i in range(28) ]
    #The amount that needs to be stored
    SQ_new_list = [] 
    SZ_new_list = []  
    SQD_new_list = []  
    RQ_new_list = [0] 
    RZ_new_list = [0]
    ht_list = []
    LQ_new_list = [0]
    LQD_new_list = [0]
    LZ_new_list = [0]
        
    for t in range(len(WU_I_new_list)):#We start at time T
        #mild step1-1：Determine the initial state
        DQ_new = round(aerfaQ * WU_I_new_list[t]) #new number of mild cases
        DQ = DQ + DQ_new - RQ_new #Total number of mild cases
        BQ = BQ_new_list[t]   #Number of empty bed
        #step1-2：
        DQ_new_wait= DQ_to_DQ_wait(DQ_new) 
        DQ_new_wait_max = numpy.row_stack((DQ_new_wait_max,DQ_new_wait))
        #DQ_wait=[DQ_wait[i]+DQ_new_wait[i] for i in range(len(DQ_new_wait))]#updata DQ_wait

        #mild step1-3:To calculate the number of new mild cases admitted
        SQ_new = min(BQ,DQ) #The number of people admitted to the hospital at time T
        # #SDQ_new_wait = SDQ_New_wait(SQ_new,DQ_wait) 
        # SDQ_new_wait = SDZ_New_wait(SQ_new,DQ_wait)  
        # DQ_wait=[DQ_wait[i]-SDQ_new_wait[i] for i in range(len(DQ_wait))]
        DQ_new_wait_max = Rule_new_FF(DQ_new_wait_max,SQ_new) 
        #mild step1-4 up dateDQ
        DQ = DQ - SQ_new
       #v step2-1：The determination of the initial state
        DZ_new = round(aerfaZ * WU_I_new_list[t]) #Calculate the number of new severe cases 
        DZ = DZ +   DZ_new  - RZ_new #Update the number of critically ill patients to be admitted
        BZ = BZ_new_list[t]   
        
        #Critical partstep2-2：Update critically ill patients awaiting admission
        DZ_wait[13] = DZ_wait[13] + DZ_new
        
        #Critical partstep step2-3:To calculate the number of newly hospitalized patients with severe diseases
        SZ_new = min(DZ,BZ) #Number of patients admitted to designated hospitals at time t
        SDZ_new_wait = SDZ_New_wait(SZ_new,DZ_wait)
        DZ_wait=[DZ_wait[i]-SDZ_new_wait[i] for i in range(len(DZ_wait))]
        
        #Critical part step2-4 up date DZ
        DZ = DZ - SZ_new
        BZ = BZ - SZ_new
        # step3
        #Determine if it can be treated
        if BZ*DQ > 0:
            #step3-1：up date DQ
            SQD_new = min (BZ,DQ)
        else:
            SQD_new = 0
            
        #SQD_reduce_wait = SDQ_New_wait(SQD_new,DQ_wait)  
        #SQD_reduce_wait = SDZ_New_wait(SQD_new,DQ_wait)  
        #DQ_wait=[DQ_wait[i]-SQD_reduce_wait[i] for i in range(len(DQ_wait))] 
        DQ_new_wait_max = Rule_new_FF(DQ_new_wait_max,SQD_new)
        #step3-2：Update patients to be admitted
        DQ = DQ - SQD_new
         #ht
        ht = SQ_new + SZ_new + SQD_new 
        
        #DQ_wait
        RQ_new = sum(DQ_new_wait_max[1:,0]) 
        DQ_new_wait_max = DQ_new_wait_max[:,1:] 
    
        back_0 = numpy.ones((DQ_new_wait_max.shape[0],1))*0
        DQ_new_wait_max = numpy.column_stack((DQ_new_wait_max,back_0))
        #DZ_wait State transition
        RZ_new = DZ_wait[0]
        del DZ_wait[0]
        DZ_wait = DZ_wait +[0]
        
        SQ_new_list.append(SQ_new) #save data
        SZ_new_list.append(SZ_new) #save data
        SQD_new_list.append(SQD_new) #save data
        RQ_new_list.append(RQ_new) #save data
        RZ_new_list.append(RZ_new) #save data
        ht_list.append(ht) #save data
        
        
        #Treat patients with mild symptoms in fangcang 
        SQ_new_wait= SQ_to_SQ_wait(SQ_new) 
        SQ_wait=[SQ_wait[i]+SQ_new_wait[i] for i in range(len(SQ_wait))]
        SQ_wait_max_1 = numpy.row_stack((SQ_wait_max_1,SQ_wait))  #save data
        #Designated hospitals admit patients with serious illness
        SZ_new_wait= SZ_to_SZ_wait(SZ_new) 
        SZ_wait=[SZ_wait[i]+SZ_new_wait[i] for i in range(len(SZ_wait))]
        SZ_wait_max = numpy.row_stack((SZ_wait_max ,SZ_wait))  
        #Designated hospitals treat mild patients
        SQD_new_wait = SQ_to_SQ_wait(SQD_new) 
        SQD_wait=[SQD_wait[i]+SQD_new_wait[i] for i in range(len(SQ_wait))]
        SQD_wait_max = numpy.row_stack((SQD_wait_max,SQD_wait))  ##save data
        
        ########Calculation of discharge
        #SQ_wait state transition
        LQ_new = SQ_wait[0]
        del SQ_wait[0]
        SQ_wait = SQ_wait +[0]
          #SQD_waitstate transition
        LQD_new = SQD_wait[0]
        del SQD_wait[0]
        SQD_wait = SQD_wait +[0]
        #SZ_waitstate transition
        LZ_new = SZ_wait[0]
        del SZ_wait[0]
        SZ_wait = SZ_wait +[0]
        
        #save data 
        LQ_new_list.append(LQ_new)
        LQD_new_list.append(LQD_new)
        LZ_new_list.append(LZ_new)
        
        
    return ht_list,RQ_new_list,RZ_new_list,LQ_new_list,LQD_new_list,LZ_new_list

def solve_beta(t,I,WU_I_new,lognormal,Et):#Solve quadratic equations with one variable
    #########def  abc 
    a = ((lognormal[0]+lognormal[1]+lognormal[2]+lognormal[3])*(lognormal[1]+lognormal[2])*(I[t]+Et[0]+Et[1])+
    (lognormal[0]+lognormal[1]+lognormal[2])*((lognormal[2]+lognormal[3])*(I[t]+Et[0]+Et[1])+(lognormal[1]+lognormal[2])*(I[t+1]+Et[1]+Et[2]))+
    (lognormal[0]+lognormal[1])*((lognormal[3]+lognormal[4])*(I[t]+Et[0]+Et[1])+(lognormal[2]+lognormal[3])*(I[t+1]+Et[1]+Et[2])+(lognormal[1]+lognormal[2])*(I[t+2]+Et[2]+Et[3]))+
    lognormal[0]*((lognormal[4]+lognormal[5])*(I[t]+Et[0]+Et[1])+(lognormal[3]+lognormal[4])*(I[t+1]+Et[1]+Et[2])+(lognormal[2]+lognormal[3])*(I[t+2]+Et[2]+Et[3])+(lognormal[1]+lognormal[2])*(I[t+3]+Et[3]+Et[4]))
    )
    
    b = ((lognormal[0]+lognormal[1]+lognormal[2]+lognormal[3]+lognormal[4])*(I[t]+Et[0]+Et[1])+
    (lognormal[0]+lognormal[1]+lognormal[2]+lognormal[3])*(I[t+1]+Et[1]+Et[2])+
    (lognormal[0]+lognormal[1]+lognormal[2])*(I[t+2]+Et[2]+Et[3])+
    (lognormal[0]+lognormal[1])*(I[t+3]+Et[3]+Et[4])+
    (lognormal[0])*(I[t+4]+Et[4]+Et[5]))
    
    c =(Et[0]+Et[1]+Et[2]+Et[3]+Et[4]
        -(WU_I_new[t+5]+WU_I_new[t+1]+WU_I_new[t+2]+WU_I_new[t+3]+WU_I_new[t+4]))
    
    d = (b**2)-(4*a*c)
    
    if d < 0:
        beta = 0.00001
    if d >= 0:
        solve_1 = (-b+math.sqrt(d))/(2*a) 
        solve_2 = (-b-math.sqrt(d))/(2*a) 
        if solve_1>0:
            beta=solve_1
        else:
            beta= 0
            
    return beta

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
def calculate_lognormal(Enew,times):#estimate lognormal distribution
    Et_max = [i for i in range(21)]
    for i in range(times):
        E_new_wait= E_to_E_wait(Enew)
        Et_max = numpy.row_stack((Et_max,E_new_wait))
    Et_max = Et_max[1:,:]
    Et_sum = [i for i in range(21)] 
    lognormal_estimation =  [i for i in range(21)] 
    for t in range(len(Et_sum)):
         Et_sum[t] = sum(Et_max[:,t])
    for t in range(len(lognormal_estimation)):
         lognormal_estimation[t] = Et_sum[t]/sum(Et_sum)
    return lognormal_estimation
def E_to_E_wait_10time_mean(E_new,times):
    E_new_wait_max = [i for i in range(21)]
    E_new_wait_mean = [i for i in range(21)]
    for i in range(times):
        E_new_wait= E_to_E_wait(E_new)
        E_new_wait_max = numpy.row_stack((E_new_wait_max,E_new_wait))
    E_new_wait_max = E_new_wait_max[1:,:]
    for i in range(len(E_new_wait_mean)):
        E_new_wait_mean[i] = numpy.mean(E_new_wait_max[:,i])
    return E_new_wait_mean
    
    
if __name__ == '__main__':
    ###Read known data### 
    bed_number_data = numpy.genfromtxt("bed_number.csv",  dtype = str , 
                                delimiter=",", skip_header = 1) #Read known Inewt and bed number
    WU_I_new_list_all = bed_number_data[:,1] #known Inewt
    BZ_new_list_all = bed_number_data[:,4] # Number of empty beds in designated hospitals
    BQ_new_list_all = bed_number_data[:,5] # Number of empty beds in fangcang  hospitals
    #data type conversion 
    WU_I_new_list_all = numpy.array(WU_I_new_list_all, dtype=numpy.float32)
    BZ_new_list_all = numpy.array(BZ_new_list_all, dtype=numpy.float32)
    BQ_new_list_all = numpy.array(BQ_new_list_all, dtype=numpy.float32)
    
    ##########Read logNormal  distribution  data#############
    lognormal_data = numpy.genfromtxt("Incubation_period_lognormal_list.csv", 
                         dtype=str,  delimiter=",", skip_header = 1)
    lognormal =  numpy.array(lognormal_data[:,1], dtype=numpy.float32)
    ##########Access to population data################
    population_data = numpy.genfromtxt("movement of population.csv", 
                         dtype=str,  delimiter=",", skip_header = 1)
    population_in_all = population_data[:,1]  
    population_out_all = population_data[:,2] 
    population_in_all =  numpy.array(population_in_all, dtype=numpy.float32)
    population_out_all =  numpy.array(population_out_all, dtype=numpy.float32)
    
    ########## read estimated Enew data###############
    Enew_estimation_data = numpy.genfromtxt("ENEW_ estimation.csv", 
                         dtype=str,  delimiter=",", skip_header = 1)
    Enew_estimation_list = Enew_estimation_data[:,1]
    Enew_estimation_list = numpy.array(Enew_estimation_list, dtype=numpy.float32)
    #Store the initial value
    ##########Read the initial value data of Et
    Et_random_data = numpy.genfromtxt("Et_random.csv", 
                         dtype=str,  delimiter=",", skip_header = 1)
    Et_random =  numpy.array(Et_random_data[:,1:], dtype=numpy.float32)
    #n= int(numpy.random.uniform(0,50-0.00001,1))
    #Et = Et_random[n]
    #
    times = 50
    Et = Et_random[6]
    N = 11212000 # total population
    S = 11211861  
    #To store data
    use_I_new_t_all = WU_I_new_list_all
    WU_I_new_list = WU_I_new_list_all[20:]
    Et_Max=[i for i in range(1,22)] 
    E_new_t=[]  
    E_new_wait_MAX=[i for i in range(1,22)]  
    beta_list_max=[i for i in range(times)]  
    #recorded data
    record_Eout = []
    record_Sout = []
    record_S = [11211861] 
    record_N = [11212000]
    beta_mean_list = []
    #step1：
    #step1-1:estimate ht  RQ RZ 
    ht_list_all, RQ_new_list_all,RZ_new_list_all,LQ_new_list,LQD_new_list,LZ_new_list  = ht_model(use_I_new_t_all,BZ_new_list_all,BQ_new_list_all)
    #step1-2:Intercept the appropriate date
    population_in = population_in_all[20:]
    population_out = population_out_all[20:]
    ###Main program start####
    for t in range(68):
        use_I_new_t = use_I_new_t_all[20:]
        #up date ht
        ht_list_all, RQ_new_list_all,RZ_new_list_all,LQ_new_list,LQD_new_list,LZ_new_list = ht_model(use_I_new_t_all,BZ_new_list_all,BQ_new_list_all)
        #step1-2:ht
        ht_list = ht_list_all[20:]
        RQ_new_list =  RQ_new_list_all[20:]
        RZ_new_list =  RZ_new_list_all[20:]
        use_I_new_t = use_I_new_t_all[20:]
        #step2-1:dispose ht
        
        for a in range(45):
            ht_list_all[45] = ht_list_all[45] + ht_list_all[a] 
            ht_list_all[a]=0
        ht_list_all[45] =ht_list_all[45] -28
        ht_initial =[0]*20+ [0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,4,1,0,0,6,3,4,6,1,0]
        for b in range(45):
            ht_list_all[b] = ht_list_all[b] + ht_initial[b]
        #step2-2:update It
        
        I = 0 
        I_list_all = [0]
        for c in range(len(ht_list_all)): 
            I = I + use_I_new_t_all[c] - ht_list_all[c] - RQ_new_list_all[c] - RZ_new_list_all[c]
            I_list_all.append(I)
        I_list = I_list_all[20:]
        #step2-3：calculate beta 
        print(t)
        beta_list =[] 
        for i in range(times):
            #step3-1:get lognormal_estimation
            lognormal_estimation = calculate_lognormal(Enew_estimation_list[t],times)
            beta = solve_beta(t,I_list,use_I_new_t,lognormal_estimation,Et)
            if beta<0:
                beta=0
            beta_list.append(beta)
            print(t,i)
        beta_list_max =numpy.row_stack((beta_list_max,beta_list)) 
        beta_mean = numpy.mean(beta_list)
        beta_mean_list.append(beta_mean)
        
        #step2-4:calculate E_new(t) 
        E_new = round((S/N)*beta_mean*(Et[0]+Et[1]+I_list[t]))
        E_new_t=E_new_t+[E_new] 
        #step2-5:
        #E_new_wait= E_to_E_wait(E_new) 
        E_new_wait = E_to_E_wait_10time_mean(E_new,times) #average
        E_new_wait=numpy.rint(E_new_wait)
        E_new_wait_MAX = numpy.row_stack((E_new_wait_MAX,E_new_wait))
        #step2-6:
        our_I_new=Et[0]+E_new_wait[0]
        use_I_new_t_all[t+21] = our_I_new
        #step2-7update Et and store:
        Et_Max=numpy.row_stack((Et_Max,Et))
        for i in range(20):
            Et[i]=Et[i+1]+E_new_wait[i+1]
            Et[20]=0
       #Consider population mobility
       #only consider E and S 
        #Step 1: Calculate E S DQ DZ at the end of the day
        E = sum(Et) #The total number of people in incubation
        S = S-E_new # the total number of susceptible people
        #Step 2: Calculate Eout Sout at the end of the day
        Eout = round((E/N)*population_out[t])
        Sout = round((S/N)*population_out[t])
        #Step 3: Update the status of E, S, DQ and DZ at the end of the day
        #1.Update Et 
        Et = Et_out(Et,Eout)  
        


        S=S+population_in[t]-Sout
        N=N+population_in[t]-population_out[t]
        record_Eout.append(Eout)
        record_Sout.append(Sout)
        record_S.append(S)
        record_N.append(N)
        
        
        #####End of main program#########
