# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 08:15:32 2016

@author: Yaseen Hull

Inverse Problem: Determine ellipsoidal distance, forward and reverse azimuth
"""
import math

d1= {}
dDis = {}
values=[]
f = open("viper.txt","r")
data = f.read() 
spl = data.splitlines()
i = 0

while i < len(spl): 
    for i in spl:
        values = i.split('\t')
        staFrm =  values.pop(0)
        values = map(float,values)
        glat =  values.pop(0)
        glong = values.pop(0)
        eS = values.pop(0)
        fAzi = values.pop(0)
        bAzi = values.pop(0)        
        d1[staFrm] = [glat,glong,eS,fAzi,bAzi]
        update3 = dDis.update(d1)
        
a = 6378137
b = 6356752.314 
f = 1/298.257224 
e2 = ((a**2)-(b**2))/(a**2)
e2x = ((a**2)-(b**2))/(b**2)
n = 0.0001679220 


x = dDis.items()

f2 = open("red.txt","w")

for i in range(len(x)):
    lat1= (x[i][1][0]) #obtaining lat
    fAzi1 = (x[i][1][3]) #forward azimuth to point 2   
    name = x[i][0] #name of point 2'''
    S  = x[i][1][2]

    phi2 = lat1
    alpha2 = fAzi1    


    k = 20
    while k>0:
        k -= 1
        mLat =((phi2 + lat1))/2 #mid lat of rbay and hnus
        mfAzi =((fAzi1 + alpha2))/2
                   
       # mid azi of rbay and hnus
              
        M = (a*(1-e2))/math.pow((1-(e2*math.pow(math.sin(mLat),2))),3.0/2.)
        N = (a)/math.sqrt(1-e2*(math.sin(mLat))**2)
        Tau = math.tan(mLat)
        eta = e2x*(math.cos(mLat))**2        
                    
        dphi = (S/M) * math.cos(mfAzi) + (S**3/ (24 * (M**3)) )*1/( (1+eta)**2 )*((2 + 3*(Tau**2) + 2*eta)*(math.sin(mfAzi)**2)*math.cos(mfAzi)+3*(-eta + 3*eta*(Tau**2) - eta**2)*(math.cos(mfAzi)**3 ))    
        dA = (S/N)*(Tau*math.sin(mfAzi))+((S**3)/(24*N**3))*Tau*((2+(Tau**2)+2*eta)*(math.sin(mfAzi)**3)+(2+7*eta+9*eta*(Tau**2)+5*(eta**2)*math.sin(mfAzi)*(math.cos(mfAzi)**2)))
            
            
        phi2 = lat1 + dphi
        alpha2 = mfAzi + 0.5*(dA)
        Deg =math.degrees(phi2)
        
        DD = int(Deg)
        MM = int((Deg - DD)*60)
        SS = round((((Deg - DD)*60)-MM),2)
        Deg2 = math.degrees(lat1)
        
    dlon_1 = (S/(N*math.cos(mfAzi)))*math.sin(mfAzi) +( (S**3)/(24*(N**3)*math.cos(mfAzi)) )*((Tau**2)*(math.sin(mfAzi)**3)+(-1 - eta + 9*eta*(Tau**2))*math.sin(mfAzi)*(math.cos(mfAzi)**2))    
    print(name+' '+str(DD)+'"'+str(MM)+"'"+str(SS)+'"')
    f2.write(name+' '+str(DD)+'"'+str(MM)+"'"+str(SS)+'"'+'\n')
#print(str(Deg))

f2.close()