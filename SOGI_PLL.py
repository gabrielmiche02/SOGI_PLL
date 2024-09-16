import pyComtrade
import comtrade
import matplotlib.pyplot as plt
import pandas as pd
import ClarkePark
import numpy as np
import control

rec = comtrade.load('COMTRADE/TesteASCII_Caso_0001_S.cfg','COMTRADE/TesteASCII_Caso_0001_S.dat')
df = rec.to_dataframe()
print("Trigger time = {}s".format(rec.trigger_time))

wdw_length = 256

# Va = df['CA00 - VA00'].iloc[:256]
# Vb = df['CA01 - VB00'].iloc[:256]
# Vc = df['CA02 - VC00'].iloc[:256]
Va = df['CA00 - VA00']
Vb = df['CA01 - VB00']
Vc = df['CA02 - VC00']

ke = 1


Valpha = np.array([])
Vbeta = np.array([])
Va_l = np.array([])
Va_ql = np.array([])
Va_l1 = np.array([])
Va_ql1 = np.array([])
Vb_l = np.array([])
Vb_ql = np.array([])
Valpha_pos = np.array([])
Vbeta_pos = np.array([])
Vd = np.array([])
Vq = np.array([])
X = np.array([])
theta = np.array([])
w = np.array([])
erro = np.array([])
freq = 2*np.pi*50

Valpha = np.append(Valpha,[0]*5500)
Vbeta = np.append(Vbeta, [0]*5500)
Va_l = np.append(Va_l,[0]*5500)
Va_ql = np.append(Va_ql,[0]*5500)
Va_l1 = np.append(Va_l1,[0]*5500)
Va_ql1 = np.append(Va_ql1,[0]*5500)
Vb_l = np.append(Vb_l,[0]*5500)
Vb_ql = np.append(Vb_ql,[0]*5500)
Valpha_pos = np.append(Valpha_pos,[0]*5500)
Vbeta_pos = np.append(Vbeta_pos,[0]*5500)
Vd = np.append(Vd,[0]*5500)
Vq = np.append(Vq,[0]*5500)
X = np.append(X,[0]*5500)
theta = np.append(theta,[0.89]*5500)
w = np.append(w,[314]*5500)
erro = np.append(erro,[0.001]*5500)

#Variables to PI
kp = 0.040
ki = 0.00005
tau = 0.0000015
ang = 0.89

for i in range(3,len(df['CA00 - VA00'])):
    # w = w[i-1]
    #calculate Valpha and Vbeta from Vabc
    Valpha[i] = 2/3 * (Va.iloc[i] - 0.5*Vb.iloc[i] - 0.5*Vc.iloc[i])    #calculate a Valpha from Vabc
    Vbeta[i] = 2/3 * (0 + 0.5*np.sqrt(3)*Vb.iloc[i] - 0.5*np.sqrt(3)*Vc.iloc[i])    #Calculate Vbeta from Vabc
   
    #calculate the sample time (não sei se precisa - talvez deixar fixo)
    T = Va.index[i] - Va.index[i-1]   

    # calculate Valpha', qValpha', Vbeta' e qVbeta' from Valpha and Vbeta
    # Va_l[i] = (2*ke*w[i-1]*T*(Valpha[i] - Valpha[i-2]) + (8 - 2*(pow(w[i-1]*T,2)))*Va_l[i-1] - (pow(w[i-1]*T,2)-2*ke*w[i-1]*T+4)*Va_l[i-2]) / (pow(w[i-1]*T,2)+2*ke*w[i-1]*T+4)
    # Va_ql [i] = ((8-2*pow(w[i-1]*T,2))*Va_ql[i-1] + (2*ke*w[i-1]*T - 4 - pow(w[i-1]*T,2))*Va_ql[i-2] + ke*pow(w[i-1]*T,2)*(Va.iloc[i] - 2*Va.iloc[i-1] + Va.iloc[i-2]))/(2*ke*w[i-1]*T + 4 + pow(w[i-1]*T,2))    
    # Vb_l [i] = (2*ke*w[i-1]*T*(Vbeta[i] - Vbeta[i-2]) + (8 - 2*(pow(w[i-1]*T,2)))*Vb_l[i-1] - (pow(w[i-1]*T,2)-2*ke*w[i-1]*T+4)*Vb_l[i-2]) / (pow(w[i-1]*T,2)+2*ke*w[i-1]*T+4)
    # Vb_ql [i] = ((8-2*pow(w[i-1]*T,2))*Vb_ql[i-1] + (2*ke*w[i-1]*T - 4 - pow(w[i-1]*T,2))*Vb_ql[i-2] + ke*pow(w[i-1]*T,2)*(Vb.iloc[i] - 2*Vb.iloc[i-1] + Vb.iloc[i-2]))/(2*ke*w[i-1]*T + 4 + pow(w[i-1]*T,2))

    Va_l[i] = ((8/pow(T,2) - 2*pow(w[i-1],2))*Va_l[i-1] + (2*ke*w[i-1]/T - 4/pow(T,2) - pow(w[i-1],2))*Va_l[i-2] + 2*ke*w[i-1]/T*(Valpha[i]-Valpha[i-2]))/(2*ke*w[i-1]/T + 4/pow(T,2) + pow(w[i-1],2))
    Va_ql[i] = ((8/pow(T,2) - 2*pow(w[i-1],2))*Va_ql[i-1] + (2*ke*w[i-1]/T - 4/pow(T,2) - pow(w[i-1],2))*Va_ql[i-2] + 2*ke*w[i-1]/T*(Valpha[i]-2*Valpha[i-1]+Valpha[i-2]))/(2*ke*w[i-1]/T + 4/pow(T,2) + pow(w[i-1],2))
    Vb_l[i] = ((8/pow(T,2) - 2*pow(w[i-1],2))*Vb_l[i-1] + (2*ke*w[i-1]/T - 4/pow(T,2) - pow(w[i-1],2))*Vb_l[i-2] + 2*ke*w[i-1]/T*(Vbeta[i]-Vbeta[i-2]))/(2*ke*w[i-1]/T + 4/pow(T,2) + pow(w[i-1],2))
    Vb_ql[i] = ((8/pow(T,2) - 2*pow(w[i-1],2))*Vb_ql[i-1] + (2*ke*w[i-1]/T - 4/pow(T,2) - pow(w[i-1],2))*Vb_ql[i-2] + 2*ke*w[i-1]/T*(Vbeta[i]-2*Vbeta[i-1]+Vbeta[i-2]))/(2*ke*w[i-1]/T + 4/pow(T,2) + pow(w[i-1],2))
  
    # Va_l[i] = (2*ke*freq*T*(Valpha[i] - Valpha[i-2]) + (8 - 2*(pow(freq*T,2)))*Va_l[i-1] - (pow(freq*T,2)-2*ke*freq*T+4)*Va_l[i-2]) / (pow(freq*T,2)+2*ke*freq*T+4)
    # Va_ql [i] = ((8-2*pow(freq*T,2))*Va_ql[i-1] + (2*ke*freq*T - 4 - pow(freq*T,2))*Va_ql[i-2] + ke*pow(freq*T,2)*(Va.iloc[i] - 2*Va.iloc[i-1] + Va.iloc[i-2]))/(2*ke*freq*T + 4 + pow(freq*T,2))
    # Vb_l [i] = (2*ke*freq*T*(Vbeta[i] - Vbeta[i-2]) + (8 - 2*(pow(freq*T,2)))*Vb_l[i-1] - (pow(freq*T,2)-2*ke*freq*T+4)*Vb_l[i-2]) / (pow(freq*T,2)+2*ke*freq*T+4)
    # Vb_ql [i] = ((8-2*pow(freq*T,2))*Vb_ql[i-1] + (2*ke*freq*T - 4 - pow(freq*T,2))*Vb_ql[i-2] + ke*pow(freq*T,2)*(Vb.iloc[i] - 2*Vb.iloc[i-1] + Vb.iloc[i-2]))/(2*ke*freq*T + 4 + pow(freq*T,2))


    #calculate the positive sequence of Valpha and Vbeta
    Valpha_pos[i] = Va_l[i] - Vb_ql[i]
    Vbeta_pos[i] = Va_ql[i] + Vb_l[i]



    #plark transformation to get Vd+ and Vq+ from Valpha+ and Vbeta+
    # Vd[i] = np.cos(freq*Va.index[i]+ang)*Valpha_pos[i] + np.sin(freq*Va.index[i]+ang)*Vbeta_pos[i]
    # Vq[i] = -np.sin(freq*Va.index[i]+ang)*Valpha_pos[i] + np.cos(freq*Va.index[i]+ang)*Vbeta_pos[i]
    
    # Vd[i] = np.cos(w[i-1]*Va.index[i-1]+theta[i-1])*Valpha_pos[i] + np.sin(w[i-1]*Va.index[i-1]+theta[i-1])*Vbeta_pos[i]
    # Vq[i] = -np.sin(w[i-1]*Va.index[i-1]+theta[i-1])*Valpha_pos[i] + np.cos(w[i-1]*Va.index[i-1]+theta[i-1])*Vbeta_pos[i]

    Vd[i] = np.cos(w[i-1]*Va.index[i-1]+ang)*Valpha_pos[i] + np.sin(w[i-1]*Va.index[i-1]+ang)*Vbeta_pos[i]
    Vq[i] = -np.sin(w[i-1]*Va.index[i-1]+ang)*Valpha_pos[i] + np.cos(w[i-1]*Va.index[i-1]+ang)*Vbeta_pos[i]

    # X[i] = ((kp*T + 2*kp*tau)*Vq[i] + (kp*T - 2*kp*tau)*Vq[i-1] + 2*tau*X[i-2])/(2*tau)
    
    erro[i] = kp*Vq[i] + ki*Vq[i-1]
    w[i] = erro[i] + freq
    theta[i] = T/2*(w[i] + w[i-1]) + theta[i-1]



# Valpha = ((2*ke*w*T*(Valpha[i] - Valpha[i-2])-((pow(w*T,2)-8)*Valpha[i-1]))-(pow(w*T,2)-2*w*T+4)*Valpha[i-2])/(pow(w*T,2)+2*k*w*T+4)
tam = 5200
plt.figure
plt.plot(Va.index[:tam],Valpha[:tam], label='Valpha')
plt.plot(Va.index[:tam],Vbeta[:tam], label='Vbeta')
plt.plot(Va.index[:tam],Va_l[:tam], label='Va_l')
plt.plot(Va.index[:tam],Va_ql[:tam], label='Va_ql')
# plt.plot(Va.index[:tam],Va_l1[:tam], label='Va_ql1')
# plt.plot(Va.index[:tam],Va_ql1[:tam], label='Va_ql1')
# plt.plot(Va.index[:tam],Vb_l[:tam], label='Vb_l')
# plt.plot(Va.index[:tam],Vb_ql[:tam], label='Vb_ql')
# plt.plot(Va.index[:tam],Valpha_pos[:tam], label='Valpha_pos')
# plt.plot(Va.index[:tam],Vbeta_pos[:tam], label='Vbeta_pos')
# plt.plot(Va.index[:tam],Vd[:tam], label='Vd')
# plt.plot(Va.index[:tam],Vq[:tam], label='Vq')
# plt.plot(Va.index[:tam],X[:tam], label='X')
# plt.plot(Va.index[:tam],w[:tam], label='w')
# plt.plot(Va.index[:tam],theta[:tam], label='theta')
# plt.plot(Va.index[:tam],erro[:tam], label='erro')
# plt.plot(Va.index,Va, label='Va')
# plt.plot(Va.index,Vb, label='Vb')
# plt.plot(Va.index,Vc, label='Vc')
plt.legend()
plt.grid()
plt.show()
