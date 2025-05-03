import numpy as np
from scipy.optimize import fsolve 
import numpy.linalg as la 
from scipy.optimize import least_squares 
#lewismatheson
def get_activity_coeff(x,T):
	xi=x
	R_k=np.array([0.9011,0.6744,1.528,1.4311,0.9200,1.380])
	nk=np.array([[1,16,1,0,0,0],[0,0,0,1,0,0],[2,16,0,0,0,1],[0,0,0,0,1,0]])
	Q_k=np.array([0.8480,0.5400,1.5320,1.4320,1.400,1.2])
	q_ki=np.zeros(4)
	r_ki=np.zeros(4)
	for i in range (4):
		for j in range(6):
			q_ki[i]+=nk[i][j]*Q_k[j]
			r_ki[i]+=nk[i][j]*R_k[j]
	e_ki=np.zeros((6,4))
	for i in range(6):
		for j in range(4):
			e_ki[i][j]=nk[j][i]*Q_k[i]/q_ki[j]

	a_ij=np.array([[0,86.02,663.5,697.2,1318,387.1],[-35.36,0,318.9,787.6,270.6,48.33],[315,1264.00,0,339.8,-66.17,-337],[16.51,-12.52,-202.00,0,-180.95,165.7],[300.00,496.1,-14.09,289.6,0,-197.5],[529.00,1397,1179.00,171,284.4,0]])
	tau_ij=np.zeros((6,6))
	for i in range(6):
		for j in range(6):
			tau_ij[i][j]=np.exp(-a_ij[i][j]/T)
	beta_ik=np.zeros((4,6))
	for i in range(4):
		for j in range(6):
			for k in range(6):
				beta_ik[i][j] += e_ki[k][i]*tau_ij[k][j]
	theta_k=np.zeros(6)
	for i in range(6):
		for j in range(4):
			sum=0
			for k in range(4):
				sum+=xi[k]*q_ki[k]
			theta_k[i]+=xi[j]*q_ki[j]*e_ki[i][j]/sum
	s_k=np.zeros(6)
	for i in range(6):
		for j in range(6):
			s_k[i]+=theta_k[j]*tau_ij[j][i]
	J_i=np.zeros(4)
	for i in range(4):
		sum=0
		for j in range(4):  
			sum+=r_ki[j]*xi[j]
		J_i[i]=r_ki[i]/sum
	l_i=np.zeros(4)
	for i in range(4):
		sum=0
		for j in range(4):  
			sum+=q_ki[j]*xi[j]
		l_i[i]=q_ki[i]/sum
	lngamc=np.zeros(4)
	lngamr=np.zeros(4)
	lngam=np.zeros(4)
	for i in range(4):
		lngamc[i]=(1-J_i[i]+np.log(J_i[i])-(5*q_ki[i]*(1-(J_i[i]/l_i[i]))+np.log(J_i[i]/l_i[i])))
	for i in range(4):
		sum=0
		for j in range(6):
			sum+=(theta_k[j]*beta_ik[i][j]/s_k[j])-(e_ki[j][i]*np.log(beta_ik[i][j]/s_k[j]))    
		lngamr[i]=q_ki[i]*(1-sum)
	for i in range(4):
		lngam[i]=lngamc[i]+lngamr[i]
	gam=np.zeros(4)
	for i in range(4):
		gam[i]=np.exp(lngam[i])
	return gam
T=590
print(get_activity_coeff(np.array([0.01,0.0,0.989,0.001]),T))
ant=np.zeros((4,3))
ant[0][:]= np.array([766,97500,624])
ant[1][:]= np.array([766,35200,337])
ant[2][:]=np.array([766,96000,602])
ant[3][:]=np.array([766,40650,373])

def sat_pressure(A,B,C,T):
	satp=A*np.exp((B/8.314)*((1/C)-(1/T)))
	return satp
Pressure=766 #mmHg 
def y_i(soln,psat,xc,act):
	return np.array([((soln[0]*soln[4])-(psat[0]*xc[0]*act[0])),((soln[1]*soln[4])-(psat[1]*xc[1]*act[1])),((soln[2]*soln[4])-(psat[2]*xc[2]*act[2])),((soln[3]*soln[4])-(psat[3]*xc[3]*act[3])),(soln[0]+soln[1]+soln[2]+soln[3]-1)])
stage=0
sol=np.zeros(5)
def yvalu(xvals,T):
	psat=np.zeros(4)
	for i in range(0,4):
		psat[i]= sat_pressure(ant[i,0],ant[i,1],ant[i,2],T)
	gamb=get_activity_coeff(xvals,T)
	lower_bounds=np.array([0,0,0,0,0])
	upper_bounds=np.array([0.9999,0.9999,0.9999,0.9999,np.inf])
	solu=least_squares(y_i,x0=np.array([0.25,0.25,0.25,0.25,766]),bounds=(lower_bounds,upper_bounds),args=(psat,xvals,gamb))
	return solu
enth=np.zeros((4,2))
enth[0][0]=6
enth[1][0]=2.436
enth[2][0]=1.62
enth[3][0]=5.4
enth[0][1]=6+97.5
enth[1][1]=35.2+2.436
enth[2][1]=96+1.62
enth[3][1]=40.65+5.4
c_p=np.zeros((4,2))
c_p[0][0]=0.4
c_p[1][0]=0.812
c_p[2][0]=0.54
c_p[3][0]=0.075
c_p[0][1]=0.463
c_p[1][1]= 0.43
c_p[2][1]=0.27
c_p[3][1]=0.0365
Ljp1=909
Ljp1i=100
vflowr=32
def iter(xjp1,xb,yval,yval0,xjp0,Ta,flow,tval):
	return np.array([xjp1[0]-(((yval[0][0]*(vflowr))+(xb[0]*flow)-(xjp0*yval0[0][0]))/(vflowr+flow-xjp0)),vflowr-(((((enth[0][0]+(c_p[0][0]*(xjp1[4]-303.0)))*xjp1[0])+((enth[1][0]+(c_p[1][0]*(xjp1[4]-303.0)))*xjp1[1])+((enth[2][0]+(c_p[2][0]*(xjp1[4]-303.0)))*(xjp1[2]))+((enth[3][0]+(c_p[3][0]*(xjp1[4]-303.0)))*xjp1[3]))*(vflowr+flow-xjp0))+(((enth[0][0]+(c_p[0][0]*(tval-303.0)))*xb[0])+((enth[1][0]+(c_p[1][0]*(tval-303.0)))*xb[1])+((enth[2][0]+(c_p[2][0]*(tval-303.0)))*xb[2])+((enth[3][0]+(c_p[3][0]*(tval-303.0)))*xb[3]))-(xjp0*(((enth[0][1]+c_p[0][1]*(Ta-303))*yval0[0][0])+((enth[1][1]+c_p[1][1]*(Ta-303))*yval0[0][1])+((enth[2][1]+c_p[2][1]*(Ta-303))*yval0[0][2])+((enth[3][1]+c_p[3][1]*(Ta-303))*yval0[0][3]))))/((((enth[0][1])+(c_p[0][1]*(tval-303.0)))*yval[0][0])+((enth[1][1]+(c_p[1][1]*(tval-303.0)))*yval[0][1])+((enth[2][1]+(c_p[2][1]*(tval-303.0)))*yval[0][2])+((enth[3][1]+(c_p[3][1]*(tval-303.0)))*yval[0][3])),xjp1[1]-(((yval[0][1]*(vflowr))+(xb[1]*flow)-(xjp0*yval0[0][1]))/(vflowr+flow-xjp0)),xjp1[2]-(((yval[0][2]*(vflowr))+(xb[2]*flow)-(xjp0*yval0[0][2]))/(vflowr+flow-xjp0)),xjp1[3]-((yval[0][3]*(vflowr)+(xb[3]*flow)-(xjp0*yval0[0][3]))/(vflowr+flow-xjp0)),1-xjp1[0]-xjp1[1]-xjp1[2]-xjp1[3]])
def flowrates(xb,yval0,xjp0,tval,flow,Ta):
	yval=np.zeros((1,4))
	yval[0][:]=yvalu(xb,tval).x[:4]
	soln=least_squares(iter,x0=np.array([0.25,0.5,0.25,0.25,320]),bounds=(np.array([0.0,0.0,0.0,0.0,0.0]),np.array([0.9999,0.99999,0.99999,0.9999,1000])),args=(xb,yval,yval0,xjp0,Ta,flow,tval))
	return soln
Ta=590
xc=np.array([0.01,0,0.989,0.001])
yflow=np.zeros((1,4))
print(flowrates(xc,yflow,0,590.0,909.0,590.0).x)
stages=1
flowval=909
yval=np.zeros((1,4))
for i in range(4):
	yval[0][i]=yvalu(xc,590).x[i]
tvalu=590
xjp0=0
yval1=np.zeros((1,4))
yval1[0][:]=yvalu(xc,tvalu).x[0:4]
yval2=yval1
xcp=xc
xc=flowrates(xc,yflow,0,590,909,590).x[0:4]
Tapp=Ta
Tap=tvalu
tvalu=flowrates(xcp,yflow,0,Tap,909,Tapp).x[4]
print(tvalu,"tvalu")
Ta=Tap
while xc[2]>0.05:
	if xc[2]>0.95:
		yval2=np.zeros((1,4))
		yval2[0][:]=yval1[0][:]
		yval1[0][:]=yvalu(xc,tvalu).x[0:4]
		xcp=xc
		vflowr=32
		xc=flowrates(xc,yval2,0,tvalu,941,Ta).x[0:4]
		Tapp=Ta
		Tap=tvalu
		tvalu=flowrates(xcp,yval2,0,tvalu,941,Tapp).x[4]
		Ta=Tap
		print(Tap)
		print(tvalu)
		print("stripping")
	elif xc[2]<0.95:
		yval2=np.zeros((1,4))
		yval2[0][:]=yval1[0][:]
		yval1[0][:]=yvalu(xc,tvalu).x[0:4]
		xcp=xc
		vflowr=941
		xc=flowrates(xc,yval2,vflowr,tvalu,flowval,Ta).x[0:4]
		Tapp=Ta
		Tap=tvalu
		tvalu=flowrates(xcp,yval2,vflowr,tvalu,flowval,Tapp).x[4]
		Ta=Tap
		print(Tap)
		print(tvalu)
		print("rectifying")
	stages+=1
	print(stages)
print(flowrates(xc,yval2,0,Tap,flowval,Tapp).x,stages,flowval,"flowval")
y00=yvalu(xcp,Tap)
print(y00.x,"y00")
#entropybalanceifneeded#(vflowr+flow-xjp0)*((xjp1[0]*((0.452+(c_p[0][0]*np.log(xjp1[5]/303.0)))))+(xjp1[1])*(0.126+c_p[1][0]*np.log(xjp1[5]/303.0))+(xjp1[2])*(0.495+(c_p[2][0]*np.log(xjp1[5]/303.0)))+(xjp1[3])*(0.069+(c_p[3][0]*np.log(xjp1[5]/303.0))))-((flow)*((xb[0]*((0.452+(c_p[0][0]*np.log(Ta/303.0)))))+(xb[1])*(0.126+c_p[1][0]*np.log(Ta/303.0))+(xb[2])*(0.495+(c_p[2][0]*np.log(Ta/303.0)))+(xb[3])*(0.069+(c_p[3][0]*np.log(Ta/303.0)))))-((xjp1[4])*(yval[0][0]*((0.452+c_p[0][1]*(np.log(Ta/303.0))))+(yval[0][1]*(0.239+c_p[1][1]*np.log(Ta/303.0)))+(yval[0][2]*(0.288+c_p[2][1]*np.log(Ta/303.0)))+(yval[0][3]*(0.188+c_p[3][1]*np.log(Ta/303.0)))))+((xjp0)*(yval0[0][0]*((0.452+c_p[0][1]*(np.log(tval/303.0))))+(yval0[0][1]*(0.239+c_p[1][1]*np.log(tval/303.0)))+(yval0[0][2]*(0.288+c_p[2][1]*np.log(tval/303.0)))+(yval0[0][3]*(0.188+c_p[3][1]*np.log(tval/303.0)))
