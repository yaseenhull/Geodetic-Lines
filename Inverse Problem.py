# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 22:08:09 2016

@author: Yaseen Hull

Inverse Problem: Determine ellipsoidal distance, forward and reverse azimuth
"""
import math 
#import numpy as np

db = {}
dc = {}
values=[]
d={}
f = open("anaconda.txt","r")
data = f.read() 
spl = data.splitlines()
i = 0

while i < len(spl): 
    for i in spl:
        values = i.split('\t')
        sta =  values.pop(0)
        values = map(float,values)
        d[sta]=[values.pop(0),values.pop(0),values.pop(0)]
        d2 = {sta:(values.pop(0))-(values.pop(0)/60)-(values.pop(0)/3600)}
        update1 = db.update(d2)
        d3 = {sta:(values.pop(0))+(values.pop(0)/60)+(values.pop(0)/3600)}
        update2 = dc.update(d3)
        
f.close()



f2 =open('boa.txt','w')
f3 =open('viper.txt','w')
f2.write("staTo -StaFrm"+"\t"+"DD"+"\t"+"MM"+'\t'+"SS"+"\t"+"DD"+"\t"+"MM"+'\t'+"SS"+'\t'+'S'+"\t"+"\t"+"DD"+"\t"+"MM"+'\t'+"SS"+"\t"+"DD"+"\t"+"MM"+'\t'+"SS"+'\n') 

a = 6378137
b = 6356752.314 
f = 1/298.257224 
e2 = ((a**2)-(b**2))/(a**2)
e2x = ((a**2)-(b**2))/(b**2)
n = 0.0001679220

v = db.items()
w =dc.items()


for i in range(len(v)):
    for k in range(len(v)):
        dlat = math.radians(v[i][1]-v[k][1])
        dlong = math.radians(w[i][1]-w[k][1])
        
        dD = int(v[k][1])
        mM = int((v[k][1] - dD)*60)
        sS = round((((v[k][1] - dD)*60)-mM)*60,4)
        latTo = str(dD)+'\t'+str(-1*mM)+'\t'+str(-1*sS)
        latto =v[k][1]
        
        dl = int(w[k][1])
        ml = int((w[k][1] - dl)*60)
        sl = round((((w[k][1] - dl)*60)-ml)*60,4)
        longTo = str(dl)+'\t'+str(ml)+'\t'+str(sl)
        longto = w[k][1]
        
        nameLat = v[i][0]+'-'+v[k][0]
        nameLong = w[i][0]+'-'+w[k][0]
        nameFrm = v[k][0]
        nameTo = v[i][0]
        
        mLat = math.radians((v[i][1]+v[k][1]))/2.0
        mLong = math.radians((w[i][1]+w[k][1]))/2.0
        
        if nameLat ==nameLong and dlat!=0 and dlong!=0:
            M = (a*(1-e2))/math.pow((1-(e2*math.pow(math.sin(mLat),2))),3.0/2.)
            N = (a)/math.sqrt(1-e2*(math.sin(mLat))**2)
            Tau = math.tan(mLat)
            eta = e2x*(math.cos(mLat))**2
            
            ScosA = M*dlat+(M/24.)*dlat*(((3*eta-3*eta*math.pow(Tau,2)+3*math.pow(eta,2))*(math.pow(dlat,2)/(1+eta)))-((2+3*math.pow(Tau,2)+2*eta)*math.pow(dlong,2)*math.pow(math.cos(mLat),2)))
            SsinA = N*math.cos(mLat)*dlong+(((N*math.cos(mLat))/24.)*dlong)*(((1+eta-9*eta*math.pow(Tau,2))*((math.pow(dlat,2))/(1+eta)))-math.pow(Tau,2)*math.pow(dlong,2)*math.pow(math.cos(mLat),2))
            S = math.sqrt(ScosA**2+SsinA**2)
            Deg = math.atan(SsinA/ScosA)
               
            dA = (S/N)*(Tau*math.sin(Deg))+((S**3)/(24*N**3))*Tau*((2+(Tau**2)+2*eta)*(math.sin(Deg)**3)+(2+7*eta+9*eta*(Tau**2)+5*(eta**2)*math.sin(Deg)*(math.cos(Deg)**2)))
            DegF = (Deg-(0.5*dA))
            DegB = (Deg+(0.5*dA))
            
           
            x = math.pi
            if SsinA>0 and ScosA<0:
                DegF = math.degrees(DegF +2*x)
                DegB = math.degrees(DegB +x)
            elif SsinA<0 and ScosA<0:
                DegF =math.degrees(DegF)
                DegB =math.degrees(DegB+x)
            elif SsinA<0 and ScosA>0:
                DegF = math.degrees(DegF +x)
                DegB = math.degrees(DegB +2*x)            
            else:
                DegF =math.degrees(DegF+x)
                DegB = math.degrees(DegB)
        
            DD = int(DegF)
            MM = int((DegF - DD)*60)
            SS = round((((DegF - DD)*60)-MM)*60,2)
            
            Dd = int(DegB)
            Mm = int((DegB - Dd)*60)
            Ss = round((((DegB - Dd)*60)-Mm)*60,2)
            Azi = str(DD)+'\t'+str(MM)+'\t'+str(SS)+'\t'+str(Dd)+'\t'+str(Mm)+'\t'+str(Ss)
            
            f2.write(nameLat+'\t'+latTo+'\t'+longTo+'\t'+str(S)+'\t'+Azi+'\n')
            f3.write(nameLat+'\t'+str(math.radians(latto))+'\t'+str(math.radians(longto))+'\t'+str(S)+'\t'+str(math.radians(DegF))+'\t'+str(math.radians(DegB))+'\n')
f2.close           