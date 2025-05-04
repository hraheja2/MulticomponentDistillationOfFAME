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
	lower_bounds=np.array([1e-4,1e-4,1e-4,1e-4,1e-4])
	upper_bounds=np.array([0.9999,0.9999,0.9999,0.9999,4500])
	solu=least_squares(y_i,x0=np.array([0.25,0.25,0.25,0.25,766]),bounds=(lower_bounds,upper_bounds),args=(psat,xvals,gamb))
	return solu
enth=np.zeros((4,2))
enth[0][0]=0.06
enth[1][0]=0.02436
enth[2][0]=0.0162
enth[3][0]=0.054
enth[0][1]=0.06+97.5
enth[1][1]=35.2+0.02436
enth[2][1]=96+0.0162
enth[3][1]=40.65+0.054
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
pi0=np.zeros((1,4))
pi0[0][0]=sat_pressure(ant[0,0],ant[0,1],ant[0,2],303)
pi0[0][1]=sat_pressure(ant[1,0],ant[1,1],ant[1,2],303)
pi0[0][2]=sat_pressure(ant[2,0],ant[2,1],ant[2,2],303)
pi0[0][3]=sat_pressure(ant[3,0],ant[3,1],ant[3,2],303)
def iter(xjp1,xb,yval,yval0,pi0,xjp0,Ta,flow,tval,xf,feed,q,ptd,ptu):
	sump=0
	sumd=0
	sumgb=0
	sumgf=0
	sump=0
	for i in range(0,4):
		sump+= yval[0][i]*np.log(yval[0][i]*ptu/pi0[0][i])
		sumd+= yval0[0][i]*np.log(yval0[0][i]*ptd/pi0[0][i])
		sumgb+=xb[i]*np.log(xb[i]*1)#get_activity_coeff(xb,tval)[i]
		sumgf+=xf[i]*np.log(xf[i]*1)#get_activity_coeff(xf,400)[i]
		sump+=xjp1[i]*np.log(xjp1[i])
	sump=0.00831*sump
	sumd=0.00831*sumd
	sumgf=0.00831*sumgf
	sumgb=0.00831*sumgb
	sump=0.00831*sump
	return np.array([xjp1[0]-(((yval[0][0]*(xjp1[5]))+(xb[0]*flow)-(xjp0*yval0[0][0])-(feed*(xf[0])))/(xjp1[6])),xjp1[5]-(((q)+((((enth[0][0]+(c_p[0][0]*(xjp1[4]-303.0)))*xjp1[0])+((enth[1][0]+(c_p[1][0]*(xjp1[4]-303.0)))*xjp1[1])+((enth[2][0]+(c_p[2][0]*(xjp1[4]-303.0)))*(xjp1[2]))+((enth[3][0]+(c_p[3][0]*(xjp1[4]-303.0)))*xjp1[3]))*(xjp1[6]))+(((enth[0][0]+(c_p[0][0]*(tval-303.0)))*xb[0])+((enth[1][0]+(c_p[1][0]*(tval-303.0)))*xb[1])+((enth[2][0]+(c_p[2][0]*(tval-303.0)))*xb[2])+((enth[3][0]+(c_p[3][0]*(tval-303.0)))*xb[3]))-(xjp0*(((enth[0][1]+c_p[0][1]*(Ta-303))*yval0[0][0])+((enth[1][1]+c_p[1][1]*(Ta-303))*yval0[0][1])+((enth[2][1]+c_p[2][1]*(Ta-303))*yval0[0][2])+((enth[3][1]+c_p[3][1]*(Ta-303))*yval0[0][3])))-((((enth[0][0]+(c_p[0][0]*(400.0-303.0)))*xf[0])+((enth[1][0]+(c_p[1][0]*(400.0-303.0)))*xf[1])+((enth[2][0]+(c_p[2][0]*(400.0-303.0)))*(xf[2]))+((enth[3][0]+(c_p[3][0]*(400.0-303.0)))*xf[3]))*(feed)))/((((enth[0][1])+(c_p[0][1]*(tval-303.0)))*yval[0][0])+((enth[1][1]+(c_p[1][1]*(tval-303.0)))*yval[0][1])+((enth[2][1]+(c_p[2][1]*(tval-303.0)))*yval[0][2])+((enth[3][1]+(c_p[3][1]*(tval-303.0)))*yval[0][3]))),xjp1[1]-(((yval[0][1]*(xjp1[5]))+(xb[1]*flow)-(xjp0*yval0[0][1])-(feed*xf[1]))/(xjp1[6])),xjp1[2]-(((yval[0][2]*(xjp1[5]))+(xb[2]*flow)-(xjp0*yval0[0][2])-(feed*xf[2]))/(xjp1[6])),xjp1[3]-((((yval[0][3]*xjp1[5])+(xb[3]*flow)-(xjp0*yval0[0][3])-(feed*xf[3])))/(xjp1[6])),1-xjp1[0]-xjp1[1]-xjp1[2]-xjp1[3],((xjp1[6])*((xjp1[0]*((0.452+(c_p[0][0]*np.log(xjp1[4]/303.0)))))+(xjp1[1])*(0.126+c_p[1][0]*np.log(xjp1[4]/303.0))+(xjp1[2])*(0.495+(c_p[2][0]*np.log(xjp1[4]/303.0)))+(xjp1[3])*(0.069+(c_p[3][0]*np.log(xjp1[4]/303.0)))-sump))-((flow)*((xb[0]*((0.452+(c_p[0][0]*np.log(Ta/303.0)))))+(xb[1])*(0.126+c_p[1][0]*np.log(Ta/303.0))+(xb[2])*(0.495+(c_p[2][0]*np.log(Ta/303.0)))+(xb[3])*(0.069+(c_p[3][0]*np.log(Ta/303.0)))-sumgb))-((xjp1[5])*((yval[0][0]*((0.452+c_p[0][1]*(np.log(Ta/303.0))+(0.00831*np.log(pi0[0][0]/ptu))))+(yval[0][1]*(0.239+c_p[1][1]*np.log(Ta/303.0))+(0.00831*np.log(pi0[0][1]/ptu)))+(yval[0][2]*(0.288+c_p[2][1]*np.log(Ta/303.0)+(0.00831*np.log(pi0[0][2]/ptu))))+(yval[0][3]*(0.188+c_p[3][1]*np.log(Ta/303.0)+(0.00831*np.log(pi0[0][3]/ptu)))))-sump))+((xjp0)*((yval0[0][0]*((0.452+c_p[0][1]*(np.log(tval/303.0)))+(0.00831*np.log(pi0[0][0]/ptd)))+(yval0[0][1]*(0.239+c_p[1][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][1]/ptd))))+(yval0[0][2]*(0.288+c_p[2][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][2]/ptd))))+(yval0[0][3]*(0.188+c_p[3][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][3]/ptd)))))-sumd))+((feed)*((xf[0]*((0.452+(c_p[0][0]*np.log(400/303.0)))))+(xf[1])*(0.126+c_p[1][0]*np.log(400/303.0))+(xf[2])*(0.495+(c_p[2][0]*np.log(400/303.0)))+(xf[3])*(0.069+(c_p[3][0]*np.log(400/303.0)))-sumgf))])
def flowrates(xb,yval0,xjp0,tval,flow,Ta,xf,feed,q,pi0,ptd):
	yval=np.zeros((1,4))
	yval[0][:]=yvalu(xb,tval).x[:4]
	ptu=yvalu(xb,tval).x[4]
	soln=least_squares(iter,x0=np.array([0.3,0.2,0.2,0.3,590,50,10]),bounds=(np.array([1e-4,1e-4,1e-4,1e-4,303,24,1]),np.array([0.9999,0.99999,0.99999,0.9999,2000,10000,10000])),args=(xb,yval,yval0,pi0,xjp0,Ta,flow,tval,xf,feed,q,ptd,ptu))
	return soln
Ta=590.0
feed=941
xf=np.array([0.01,0.0133,0.956,0.02])
xc=np.array([0.01,0.0000001,0.989,0.001])
xd=np.array([0,18.82/32,0.004,12.87/32])
yflow=np.zeros((1,4))
for i in range(0,4):
	yflow[0][i]=5e-10
ptd=5e-10
q0=0
sumd=0
sumf=0
sumb=0
for i in range(0,4):
	sumd+=xd[i]*(enth[i][0]+(c_p[i][0]*(320-303)))
	sumb+=xc[i]*(enth[i][0]+(c_p[i][0]*(590-303)))
	sumf+=xf[i]*(enth[i][0]+(c_p[i][0]*(400-303)))
q0=sumd+sumb-sumf
print(q0)
print(flowrates(xc,yflow,0,590.0,909.0,590.0,xf,0,q0,pi0,ptd).x)
stages=1
flowval=909
yval=np.zeros((1,4))
for i in range(4):
	yval[0][i]=yvalu(xc,590).x[i]
tvalu=590.0
xjp0=0
yval1=np.zeros((1,4))
yval1[0][:]=yvalu(xc,tvalu).x[0:4]
yval2=yval1.copy()
xcp=xc.copy()
xc=flowrates(xc,yflow,0,590,909,600,xf,0,q0,pi0,ptd).x[0:4]
print(xc,"xc")
xjp00=xjp0
xjp0=flowrates(xcp,yflow,xjp00,590,909,600,xf,0,q0,pi0,ptd).x[5]
flowval1=flowval
flowval=flowrates(xcp,yflow,xjp00,590,909,600,xf,0,q0,pi0,ptd).x[6]
Tapp=Ta
Tap=tvalu
tvalu=flowrates(xcp,yflow,0,Tap,909,600,xf,0,q0,pi0,ptd).x[4]
print(tvalu,"tvalu")
Ta=Tap
Purity=0.01
rectifying_stages=0
while xc[2]>Purity:
	if xc[2]>xf[2]:
		yval2[0][:]=yval1.copy()
		yval1[0][:]=yvalu(xc,tvalu).x[0:4]
		print(yval2,'yval2')
		print(yval1,"yval1")
		ptd0=ptd
		ptd=yvalu(xcp,Ta).x[4]
		xcp=xc.copy()
		xc=flowrates(xc,yval2,xjp0,tvalu,flowval,Ta,xf,0,0,pi0,ptd).x[0:4]
		xjp00=xjp0
		xjp0=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,0,0,pi0,ptd0).x[5]
		print(xc,"xc")
		print(xcp,"xcp")
		print(xjp0,"xjp0")
		flowval1=flowval
		flowval=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,0,0,pi0,ptd0).x[6]
		Tapp=Ta
		Tap=tvalu
		tvalu=flowrates(xcp,yval2,xjp00,tvalu,flowval1,Tapp,xf,0,0,pi0,ptd0).x[4]
		Ta=Tap
		print(Tap)
		print(tvalu)
		print("stripping")
	elif xc[2]<xf[2]:
		if rectifying_stages <=0:
			yval2[0][:]=yval1.copy()
			yval1[0][:]=yvalu(xc,tvalu).x[0:4]
			print(yval2,'yval2')
			print(yval1,"yval1")
			ptd0=ptd
			ptd=yvalu(xcp,Ta).x[4]
			xcp=xc.copy()
			xc=flowrates(xc,yval2,xjp0,tvalu,flowval,Ta,xf,feed,0,pi0,ptd).x[0:4]
			xjp00=xjp0
			xjp0=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,feed,0,pi0,ptd).x[5]
			print(xc,"xc")
			print(xcp,"xcp")
			print(xjp0,"xjp0")
			flowval1=flowval
			flowval=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,feed,0,pi0,ptd0).x[6]
			Tapp=Ta
			Tap=tvalu
			tvalu=flowrates(xcp,yval2,xjp00,tvalu,flowval1,Tapp,xf,feed,0,pi0,ptd0).x[4]
			Ta=Tap
			print(Tap)
			print(tvalu)
			print("rectifying")
			rectifying_stages+=1
		else:
			yval2[0][:]=yval1.copy()
			yval1[0][:]=yvalu(xc,tvalu).x[0:4]
			print(yval2,'yval2')
			print(yval1,"yval1")
			ptd0=ptd
			ptd=yvalu(xcp,Ta).x[4]
			xcp=xc.copy()
			xc=flowrates(xc,yval2,xjp0,tvalu,flowval,Ta,xf,0,0,pi0,ptd).x[0:4]
			xjp00=xjp0
			xjp0=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,0,0,pi0,ptd0).x[5]
			print(xc,"xc")
			print(xcp,"xcp")
			print(xjp0,"xjp0")
			flowval1=flowval
			flowval=flowrates(xcp,yval2,xjp00,tvalu,flowval,Ta,xf,0,0,pi0,ptd0).x[6]
			Tapp=Ta
			Tap=tvalu
			tvalu=flowrates(xcp,yval2,xjp00,tvalu,flowval1,Tapp,xf,0,0,pi0,ptd0).x[4]
			Ta=Tap
			print(Tap)
			print(tvalu)
			print("rectifying")
			rectifying_stages+=1
	stages+=1
	print(stages)
print(flowrates(xc,yval1,xjp0,tvalu,flowval,Ta,xf,0,0,pi0,ptd).x,stages,flowval,"flowval",xjp0,"xjp0")
y00=yvalu(xcp,Tapp)
print(y00.x,"y00")
