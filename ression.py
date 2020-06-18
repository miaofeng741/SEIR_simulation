# -*- coding: utf-8 
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

#The regression model fits the value of beta
import numpy as np  
from sklearn import linear_model
from sklearn.linear_model import LinearRegression 

def get_lr_stats(x, y, model):
    from scipy import stats
    n     = len(x)
    y_prd = model.predict(x)
    Regression = sum((y_prd - np.mean(y))**2) 
    Residual   = sum((y - y_prd)**2)          
    R_square   = Regression / (Regression + Residual) 
    F          = (Regression / 1) / (Residual / ( n - 2 )) 
    pf         = stats.f.sf(F, 1, n-2)
    ## T
    L_xx  =  n * np.var(x)
    sigma =  np.sqrt(Residual / n) 
    t     =  model.coef_ * np.sqrt(L_xx) / sigma
    pt    =  stats.t.sf(t, n-2)
    message2 = '           t ： ' + str(t[0][0])+ '；' +  '\t' + 'pt ： '  + str(pt[0][0])
    return 
if __name__ == '__main__':
    plan_1= np.loadtxt('FUPAN_DATA.csv',delimiter=',',
                     skiprows=1,usecols=(3),unpack = True)
    beta_list = plan_1[24:38]
    x = [i for i in range(1,15)]
    plt.scatter(x, beta_list)
    plt.xlabel('date')
    plt.ylabel("beta")
    plt.show()
    x_in = np.array(x).reshape(-1,1)
    y_in = np.array(beta_list).reshape(-1,1)
    lreg = LinearRegression()
    lreg.fit(x_in, y_in)
    y_prd = lreg.predict(x_in)
    get_lr_stats(x_in, y_in, lreg)
    
