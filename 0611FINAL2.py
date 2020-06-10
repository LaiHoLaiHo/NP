#for all
import numpy as np
import scipy as sp
import pylab as plt
from scipy.integrate import odeint
from scipy import fftpack


r = []
jf = []

for i in range(5):
    R = 22+0.15*i
    FT = 10 #sample total time  
    SP = 0.1*0.001 #sample period
    
    class HodgkinHuxley():
        C_m  =   1.0
        g_Na = 120.0
        g_K  =  36.0
        g_L  =   0.3
        E_Na =  50.0
        E_K  = -77.0
        E_L  = -54.387
        #FT = 5 #sample total time  
        #SP = 0.001 #sample period
        t = sp.arange(0.0, FT*1000, SP*1000)

        def alpha_m(self, V):
            return 0.1*(V+40.0)/(1.0-sp.exp(-(V+40.0)/10.0))
        def beta_m(self, V):
            return 4.0*sp.exp(-(V+65.0)/18.0)
        def alpha_h(self, V):
            return 0.07*sp.exp(-(V+65.0)/20.0)
        def beta_h(self, V):
            return 1.0/(1.0 + sp.exp(-(V+35.0)/10.0))
        def alpha_n(self, V):
            return 0.01*(V+55.0)/(1.0-sp.exp(-(V+55.0)/10.0))
        def beta_n(self, V):
            return 0.125*sp.exp(-(V+65)/80.0)
        def I_Na(self, V, m, h):
            return self.g_Na*m**3*h*(V-self.E_Na)
        def I_K(self, V, n):
            return self.g_K*n**4*(V-self.E_K)
        def I_K2(self, V, n):
            return self.g_K*n**6*(V-self.E_K)
        def I_L(self, V):
            return self.g_L*(V-self.E_L)

        @staticmethod
        def Couple(X, t, self):
            V1, m1, h1, n1, V2, m2, h2, n2 = X

            if V1 < V2:
                dV1dt=((V2-V1)/R - self.I_Na(V1, m1, h1) - self.I_K2(V1, n1) - self.I_L(V1)) / self.C_m
                dV2dt=(- (V2-V1)/R - self.I_Na(V2, m2, h2) - self.I_K(V2, n2) - self.I_L(V2)) / self.C_m
            elif V1 > V2:
                dV1dt=(- (V1-V2)/R - self.I_Na(V1, m1, h1) - self.I_K2(V1, n1) - self.I_L(V1)) / self.C_m
                dV2dt=((V1-V2)/R - self.I_Na(V2, m2, h2) - self.I_K(V2, n2) - self.I_L(V2)) / self.C_m
            else:
                dV1dt=(-self.I_Na(V1, m1, h1) - self.I_K2(V1, n1) - self.I_L(V1)) / self.C_m
                dV2dt=(-self.I_Na(V2, m2, h2) - self.I_K(V2, n2) - self.I_L(V2)) / self.C_m
            dm1dt = self.alpha_m(V1) * (1.0-m1) - self.beta_m(V1) * m1
            dh1dt = self.alpha_h(V1) * (1.0-h1) - self.beta_h(V1) * h1
            dn1dt = self.alpha_n(V1) * (1.0-n1) - self.beta_n(V1) * n1
            dm2dt = self.alpha_m(V2) * (1.0-m2) - self.beta_m(V2) * m2
            dh2dt = self.alpha_h(V2) * (1.0-h2) - self.beta_h(V2) * h2
            dn2dt = self.alpha_n(V2) * (1.0-n2) - self.beta_n(V2) * n2
            return dV1dt, dm1dt, dh1dt, dn1dt,  dV2dt, dm2dt, dh2dt, dn2dt
        def cri(self):
            X = odeint(self.Couple, [-41, 0.05, 0.6, 0.32, -39, 0.05, 0.6, 0.32], self.t, args=(self,))
            V1 = X[:,0]
            m1 = X[:,1]
            h1 = X[:,2]
            n1 = X[:,3]
            V2 = X[:,4]
            m2 = X[:,5]
            h2 = X[:,6]
            n2 = X[:,7]


            SN = FT/SP  #Sample number
            Fs = 1/SP  #sample freq
            fr = Fs/SN

            fft1=fftpack.fft(V1)
            yf1=abs(fft1)/len(V1)
            xf1=sp.arange(len(V1))*fr

            fft2=fftpack.fft(V2)
            yf2=abs(fft2)/len(V2)
            xf2=sp.arange(len(V2))*fr
            #ave = np.mean(yf2)
            
            A = []
            B = []
            for h in range(len(xf1)):
                if xf1[h] > 10:
                    if xf1[h] <100:
                        A.append(xf1[h])
                        B.append(yf1[h])
            pos = B.index(max(B))
            ffp = A[pos]
            #yf1[yf1>25] = 0
            #mmax = np.argmax(yf1)
            #ffp = xf1[mmax]
            A = []   #j th index in xf2
            B = []  # amp in yf2
            for j in range(len(xf2)):
                if xf2[j] > 10:
                    if xf2[j] <ffp-1:
                        A.append(xf2[j])
                        B.append(yf2[j])
            pos = B.index(max(B))
            #maxx = max(B)
            fp = A[pos]  
            if B[pos] < 0.2:
                C = []
                D = []
                for j in range(len(xf2)):
                    if xf2[j] > ffp-5:
                        if xf2[j] <ffp+5:
                            C.append(xf2[j])
                            D.append(yf2[j])
                pos = D.index(max(D))
                fp = C[pos]
            print("B[pos] = ")
            print(B[pos])
            return fp/ffp
            
            ####################################################
    r.append(R)
    jf.append(HodgkinHuxley().cri())
    print("fn = "+str(i))
print("go")        
plt.scatter(r,jf)
plt.show()
