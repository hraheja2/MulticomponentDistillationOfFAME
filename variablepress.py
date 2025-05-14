import numpy as np

from scipy.optimize import fsolve 
import numpy.linalg as la 
from scipy.optimize import least_squares 
#lewismatheson
ptu=50
def get_activity_coeff(x,T):
	xi=x
	R_k=np.array([0.9011,0.6744,1.4457,1.4311,0.9200,1.3013])
	nk=np.array([[1,16,0,0,0,1],[0,0,0,1,0,0],[2,16,1,0,0,0],[0,0,0,0,1,0]])
	Q_k=np.array([0.8480,0.5400,1.180,1.4320,1.400,1.2240])
	q_ki=np.zeros(4)
	r_ki=np.zeros(4)
	for i in range (4):
		for j in range(6):
			q_ki[i]+=nk[i][j]*Q_k[j]
			r_ki[i]+=nk[i][j]*R_k[j]
	e_ki=np.zeros((6,4))
	for i in range(6):
		for j in range(4):
			if nk[j][i]==0:
				e_ki[i][j]=0
			else:
				e_ki[i][j]=nk[j][i]*Q_k[i]/q_ki[j]

	a_ij=np.array([[0,0,476.4,697.2,1318,663.5],[0,0,476.4,697.2,1318,663.5],[267.6,267.6,0,108.65,472.5,669.4],[16.51,16.51,23.39,0,-180.95,-202],[300.00,300,-195.4,289.6,0,-14.09],[315.3,315.3,-297.8,339.8,-66.17,0]])
	tau_ij=np.zeros((6,6))
	for i in range(6):
		for j in range(6):
			tau_ij[j][i]=np.exp(-a_ij[j][i]/T)
	beta_ik=np.zeros((4,6))
	for i in range(4):
		for j in range(6):
			for k in range(6):
				beta_ik[i][j] += e_ki[k][i]*tau_ij[k][j]
	theta_k=np.zeros(6)
	sumsum=0
	for i in range(6):
		sumsum=0
		for j in range(4):
			sum=0
			for k in range(4):
				sum+=xi[k]*q_ki[k]
			sumsum+=xi[j]*q_ki[j]*e_ki[i][j]/sum
		theta_k[i]=sumsum
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
		lngamc[i]=(1-J_i[i]+np.log(J_i[i])-(5*q_ki[i]*((1-(J_i[i]/l_i[i]))+np.log(J_i[i]/l_i[i]))))
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
print(get_activity_coeff(np.array([0.25,0.25,0.25,0.25]),750))
s=0.96
ant=np.zeros((4,3))
ant[0][:]= np.array([766,97500,624])
ant[1][:]= np.array([766,35200,337])
ant[2][:]=np.array([766,96000,602])
ant[3][:]=np.array([766,40650,373])
def temp(xt,xvv,act,anto,P):
	satu_press=np.zeros(5)
	for i in range(0,4):
		satu_press[i]=xt[i]-(ant[i][0]*np.exp((anto[i][1]/8.314)*((1/anto[i][2])-(1/xt[4])))*xvv[i]*act[i]/P)
	satu_press[4]=1-xt[0]-xt[1]-xt[2]-xt[3]
	return satu_press
def solverT(xvv,anto,P):
	act=np.zeros(4)
	act=get_activity_coeff(xvv,450)
	Ti=450
	Tip=0
	while (abs(Ti-Tip)/Ti)>0.005:
		Tip=Ti
		Ti=least_squares(temp,x0=np.array([0.01,0.0001,0.989,0.0009,Ti]),bounds=(np.array([1e-4,1e-4,1e-4,1e-4,340]),np.array([0.9999,0.9999,0.9999,0.9999,2000])),args=(xvv,act,anto,P)).x[4]
	return(least_squares(temp,x0=np.array([0.01,0.0001,0.989,0.0009,Ti]),bounds=(np.array([1e-4,1e-4,1e-4,1e-4,340]),np.array([0.9999,0.9999,0.9999,0.9999,2000])),args=(xvv,act,anto,P)))	
def sat_pressure(A,B,C,T):
	satp=A*np.exp((B/8.314)*((1/C)-(1/T)))
	return satp
def y_i(soln,psat,xc,act):
	return np.array([((soln[0])-(psat[0]*xc[0]*act[0]/soln[4])),((soln[1])-(psat[1]*xc[1]*act[1]/soln[4])),((soln[2])-(psat[2]*xc[2]*act[2]/soln[4])), ((soln[3])-(psat[3]*xc[3]*act[3]/soln[4])),(soln[0]+soln[1]+soln[2]+soln[3]-1)])
x00=np.zeros(5)
def yvalu(xvals,T,yvalt):
	psat=np.zeros(4)
	for i in range(0,4):
		psat[i]= sat_pressure(ant[i,0],ant[i,1],ant[i,2],T)
	gamb=get_activity_coeff(xvals,T)
	xp=np.zeros(4)
	xp[:]=xvals.copy()
	x00=np.append(xp,((ptd+0)/2))
	lower_bounds=np.append(np.array([0.0001,0.0001,0.0001,0.0001]),(0))
	upper_bounds=np.append(np.array([0.9999,0.9999,0.9999,0.9999]),ptd)
	solu=least_squares(y_i,x0=x00,bounds=(lower_bounds,upper_bounds),args=(psat,xvals,gamb))
	return solu
enth=np.zeros((4,2))
enth[0][0]=0
enth[1][0]=0
enth[2][0]=0
enth[3][0]=0
enth[0][1]=97.5+(0.8*(624-303))
enth[1][1]=35.2+(0.068*(337-303))
enth[2][1]=96+(0.54*(602-303))
enth[3][1]=40.65+(0.075*(373-303))
c_p=np.zeros((4,2))
c_p[0][0]=0.8
c_p[1][0]=0.0636
c_p[2][0]=0.54
c_p[3][0]=0.075
c_p[0][1]=0.4
c_p[1][1]=0.0318
c_p[2][1]=0.27
c_p[3][1]=0.0365
def ssr(a,b):
	sum=0
	for i in range(0,4):
		sum+=((a[0][i]-b[i])**2)
	return sum**0.5
pi0=np.zeros((1,4))
pi0[0][0]=sat_pressure(ant[0,0],ant[0,1],ant[0,2],303)
pi0[0][1]=sat_pressure(ant[1,0],ant[1,1],ant[1,2],303)
pi0[0][2]=sat_pressure(ant[2,0],ant[2,1],ant[2,2],303)
pi0[0][3]=sat_pressure(ant[3,0],ant[3,1],ant[3,2],303)
def reboiler(xjp1,xb,yval,yval0,pi0,xjp0,Ta,flow,tval1,xf,feed,q,ptd,ptu,xjp10,b1):
	sum=0
	for i in range(0,4):
		sum += yimpv[i]*(enth[i][1]+(c_p[i][1]*(feed_temp-ant[i][2])))
	sumpsg=sum*(feed)*(1-qratio)
	return np.array([xjp10-(((xjp1[5])+((((enth[0][0]+(c_p[0][0]*(xjp1[4]-303.0)))*xjp1[0])+((enth[1][0]+(c_p[1][0]*(xjp1[4]-303.0)))*xjp1[1])+((enth[2][0]+(c_p[2][0]*(xjp1[4]-303.0)))*(xjp1[2]))+((enth[3][0]+(c_p[3][0]*(xjp1[4]-303.0)))*xjp1[3]))*(b1))-((((enth[0][0]+(c_p[0][0]*(tval1-303.0)))*xb[0])+((enth[1][0]+(c_p[1][0]*(tval1-303.0)))*xb[1])+((enth[2][0]+(c_p[2][0]*(tval1-303.0)))*xb[2])+((enth[3][0]+(c_p[3][0]*(tval1-303.0)))*xb[3]))*(flow))+(xjp0*(((enth[0][1]+c_p[0][1]*(Ta-624))*yval0[0][0])+((enth[1][1]+c_p[1][1]*(Ta-337))*yval0[0][1])+((enth[2][1]+c_p[2][1]*(Ta-602))*yval0[0][2])+((enth[3][1]+c_p[3][1]*(Ta-373))*yval0[0][3])))+(sumpsg+((((enth[0][0]+(c_p[0][0]*(feed_temp-303.0)))*xf[0])+((enth[1][0]+(c_p[1][0]*(feed_temp-303.0)))*xf[1])+((enth[2][0]+(c_p[2][0]*(feed_temp-303.0)))*(xf[2]))+((enth[3][0]+(c_p[3][0]*(feed_temp-303.0)))*xf[3]))*(feed*(qratio)))))/((((enth[0][1])+(c_p[0][1]*(tval1-624.0)))*yval[0][0])+((enth[1][1]+(c_p[1][1]*(tval1-337.0)))*yval[0][1])+((enth[2][1]+(c_p[2][1]*(tval1-602.0)))*yval[0][2])+((enth[3][1]+(c_p[3][1]*(tval1-373.0)))*yval[0][3]))),xjp1[1]-(((yval[0][1]*(xjp10))+(xb[1]*flow)-(xjp0*yval0[0][1])-(feed*xf[1]))/(b1)),xjp1[0]-(((yval[0][0]*(xjp10))+(xb[0]*flow)-(xjp0*yval0[0][0])-(feed*xf[0]))/(b1)),xjp1[2]-(((yval[0][2]*(xjp10))+(xb[2]*flow)-(xjp0*yval0[0][2])-(feed*xf[2]))/(b1)),xjp1[3]-((((yval[0][3]*xjp10)+(xb[3]*flow)-(xjp0*yval0[0][3])-(feed*xf[3])))/(b1)),xjp1[0]+xjp1[1]+xjp1[2]+xjp1[3]-1])
def reboilsolve(xb,yval0,xjp0,tval1,flow,Ta,xf,feed,q,qratio,pi0,ptd,s):
	yval=np.zeros((1,4))
	ytest=np.zeros((1,4))
	ytest[0][:]=xb.copy()
	yval[0][:]=yvalu(xb,tval1,ytest).x[:4]
	xjp10=(s/(1-s))*(bottoms+feed)
	if feed<=0:
		b1=xjp10+bottoms
	elif(feed>0):
		s=1-(bottoms/(((1-qratio)*feed)+(bottoms)))
		print(s,"s")
		xjp10=(1-qratio)*feed
		b1=xjp10+bottoms-feed
	xtest=xb.copy()
	t=np.append(xtest,tval1)
	t=np.append(t,1e5)
	soln=least_squares(reboiler,x0=t,bounds=(np.array([1e-4,1e-4,1e-4,1e-4,340,0]),np.array([0.9999,0.9999,0.9999,0.9999,2000,1e10])),args=(xb,yval,yval0,pi0,xjp0,Ta,flow,tval1,xf,feed,q,ptd,ptu,xjp10,b1))
	jj=np.append(soln.x,xjp10)
	jj=np.append(jj,b1)	
	return jj
def iter(xjp1,xb,yval,yval0,pi0,xjp0,Ta,flow,tval1,xf,feed,q,ptd,ptu,upflow):
	return np.array([xjp1[5]-(((q)+((((enth[0][0]+(c_p[0][0]*(xjp1[4]-303.0)))*xjp1[0])+((enth[1][0]+(c_p[1][0]*(xjp1[4]-303.0)))*xjp1[1])+((enth[2][0]+(c_p[2][0]*(xjp1[4]-303.0)))*(xjp1[2]))+((enth[3][0]+(c_p[3][0]*(xjp1[4]-303.0)))*xjp1[3]))*(upflow))-((((enth[0][0]+(c_p[0][0]*(tval1-303.0)))*xb[0])+((enth[1][0]+(c_p[1][0]*(tval1-303.0)))*xb[1])+((enth[2][0]+(c_p[2][0]*(tval1-303.0)))*xb[2])+((enth[3][0]+(c_p[3][0]*(tval1-303.0)))*xb[3]))*(flow))+(xjp0*(((enth[0][1]+c_p[0][1]*(Ta-624))*yval0[0][0])+((enth[1][1]+c_p[1][1]*(Ta-337))*yval0[0][1])+((enth[2][1]+c_p[2][1]*(Ta-602))*yval0[0][2])+((enth[3][1]+c_p[3][1]*(Ta-373))*yval0[0][3])))+((((enth[0][0]+(c_p[0][0]*(feed_temp-303.0)))*xf[0])+((enth[1][0]+(c_p[1][0]*(feed_temp-303.0)))*xf[1])+((enth[2][0]+(c_p[2][0]*(feed_temp-303.0)))*(xf[2]))+((enth[3][0]+(c_p[3][0]*(feed_temp-303.0)))*xf[3]))*(feed)))/((((enth[0][1])+(c_p[0][1]*(tval1-624.0)))*yval[0][0])+((enth[1][1]+(c_p[1][1]*(tval1-337.0)))*yval[0][1])+((enth[2][1]+(c_p[2][1]*(tval1-602.0)))*yval[0][2])+((enth[3][1]+(c_p[3][1]*(tval1-373.0)))*yval[0][3]))),xjp1[1]-(((yval[0][1]*(xjp1[5]))+(xb[1]*flow)-(xjp0*yval0[0][1])-(feed*xf[1]))/(upflow)),xjp1[2]-(((yval[0][2]*(xjp1[5]))+(xb[2]*flow)-(xjp0*yval0[0][2])-(feed*xf[2]))/(upflow)),xjp1[3]-((((yval[0][3]*xjp1[5])+(xb[3]*flow)-(xjp0*yval0[0][3])-(feed*xf[3])))/(upflow)),xjp1[0]-((((yval[0][0]*xjp1[5])+(xb[0]*flow)-(xjp0*yval0[0][0])-(feed*xf[0])))/(upflow)),1-xjp1[0]-xjp1[1]-xjp1[2]-xjp1[3]])
def flowrates(xb,yval0,xjp123,tval1,flow,Ta,xf,feed,q,pi0,ptd):
	yval=np.zeros((1,4))
	if xjp123==0:
		yval[0][i]=yvalu(xb,tval1,yval0).x[0:4]
	else:
		yval[0][:]=yvalu(xb,tval1,yval0).x[0:4]
	ptu=1e-10
	upflow=flow
	upflowp=0
	xtest=xb.copy()
	Temp_max=solverT(xb,ant,ptd)
	t=np.append(xtest,((tval1)))
	t=np.append(t,((min+max)/2))
	soln=np.zeros(6)
	while abs((upflow-upflowp)/upflow)>0.01:
		soln[:]=least_squares(iter,x0=t,bounds=(np.array([1e-4,1e-4,1e-4,1e-4,300,min]),np.array([0.9999,0.9999,0.9999,0.9999,2000,max])),args=(xb,yval,yval0,pi0,xjp123,Ta,flow,tval1,xf,feed,q,ptd,ptu,upflow)).x[:]
		upflowp=upflow
		upflow=soln[5]+flow-feed-xjp0

	lalap=np.append(soln,upflow)
	return lalap
summm=0
def xvaluuuuu(xarr,yvalt,gamm,T):
	nana=np.zeros(5)
	for i in range(0,4):
		nana[i]=(xarr[i]*gamm[i])-(yvalt[0][i]*xarr[4])/sat_pressure(ant[i][0],ant[i][1],ant[i][2],T)
	nana[4]=1-xarr[0]-xarr[1]-xarr[2]-xarr[3]
	return nana
def solvex(yvalt,T):
	ttt=np.array([0.01,0.01,0.01,0.01])
	pressure=(ptd)/2
	ttt=np.append(ttt,pressure)
	xttup=yvalt[0][:]
	sum=100
	haha=np.zeros(5)
	while abs(sum)>=1e-3:
		sum=0
		gamm=np.ones(4)#get_activity_coeff(xttup,T)
		haha=least_squares(xvaluuuuu,x0=ttt,bounds=(np.array([1e-4,1e-4,1e-4,1e-4,0]),np.array([0.9999,0.9999,0.9999,0.9999,ptd])),args=(yvalt,gamm,T)).x
		for i in range(0,4):
			sum+=xttup[i]-haha[i]
		xttup=haha[:4]
	return haha
def iterfeed(xjp1,xb,xf,yval,yval0,pi0,xjp0,Ta,flow,tval1,feed,q,ptd,ptu,upflow,xjp10):
	sum=0
	for i in range(0,4):
		sum+=yimpv[i]*(enth[i][1]+(c_p[i][1]*(feed_temp-ant[i][2])))
	sumpsg=sum*(feed)*(1-qratio)
	return np.array([xjp10-((((((enth[0][0]+(c_p[0][0]*(xjp1[4]-303.0)))*xjp1[0])+((enth[1][0]+(c_p[1][0]*(xjp1[4]-303.0)))*xjp1[1])+((enth[2][0]+(c_p[2][0]*(xjp1[4]-303.0)))*(xjp1[2]))+((enth[3][0]+(c_p[3][0]*(xjp1[4]-303.0)))*xjp1[3]))*(upflow))-((((enth[0][0]+(c_p[0][0]*(tval1-303.0)))*xb[0])+((enth[1][0]+(c_p[1][0]*(tval1-303.0)))*xb[1])+((enth[2][0]+(c_p[2][0]*(tval1-303.0)))*xb[2])+((enth[3][0]+(c_p[3][0]*(tval1-303.0)))*xb[3]))*(flow))+(xjp0*(((enth[0][1]+c_p[0][1]*(Ta-624))*yval0[0][0])+((enth[1][1]+c_p[1][1]*(Ta-337))*yval0[0][1])+((enth[2][1]+c_p[2][1]*(Ta-602))*yval0[0][2])+((enth[3][1]+c_p[3][1]*(Ta-373))*yval0[0][3])))+(sumpsg+((((enth[0][0]+(c_p[0][0]*(feed_temp-303.0)))*xf[0])+((enth[1][0]+(c_p[1][0]*(feed_temp-303.0)))*xf[1])+((enth[2][0]+(c_p[2][0]*(feed_temp-303.0)))*(xf[2]))+((enth[3][0]+(c_p[3][0]*(feed_temp-303.0)))*xf[3]))*(feed*(q)))))/((((enth[0][1])+(c_p[0][1]*(tval1-624.0)))*yval[0][0])+((enth[1][1]+(c_p[1][1]*(tval1-337.0)))*yval[0][1])+((enth[2][1]+(c_p[2][1]*(tval1-602.0)))*yval[0][2])+((enth[3][1]+(c_p[3][1]*(tval1-373.0)))*yval[0][3]))),xjp1[1]-(((yval[0][1]*(xjp10))+(xb[1]*flow)-(xjp0*yval0[0][1])-(feed*xf[1]))/(upflow)),xjp1[2]-(((yval[0][2]*(xjp10))+(xb[2]*flow)-(xjp0*yval0[0][2])-(feed*xf[2]))/(upflow)),xjp1[3]-((((yval[0][3]*xjp10)+(xb[3]*flow)-(xjp0*yval0[0][3])-(feed*xf[3])))/(upflow)),xjp1[0]+xjp1[1]+xjp1[2]+xjp1[3]-1])
def flowratesfeed(xb,yval0,xjp123,tval1,flow,Ta,xf,feed,qratio,pi0,ptd):
	yval=np.zeros((1,4))
	if xjp123==0:
		yval[0][i]=solverT(xb,ant,ptd).x[0:4]
	else:
		yval[0][:]=yvalu(xb,tval1,yval0).x[0:4]
	upflow=flow-((qratio)*feed)
	print(upflow,"upflow")
	xjp1a=xjp123+((1-qratio)*feed)
	ptu=1e-10
	xtest=xb.copy()
	t=np.append(xtest,tval1)
	solt=least_squares(iterfeed,x0=t,bounds=(np.array([1e-4,1e-4,1e-4,1e-4,340]),np.array([0.9999,0.9999,0.9999,0.9999,2000])),args=(xb,xf,yval,yval0,pi0,xjp123,Ta,flow,tval1,feed,qratio,ptd,ptu,upflow,xjp1a))
	lalap=np.append(solt.x,xjp1a)
	lalap=np.append(lalap,upflow)
	return lalap
yval22=np.zeros((1,4))
print(yval22)
ptd=2341
Ta=600
Pressure=25 #mmHg  
feed=941
feed_temp=458
cc=18.82/32
dd=12.87/32
xf=np.array([0.02,0.02,0.956,0.004])
xd=np.array([0.0001,0.8572,0.004,0.1387])
xbot=np.array([0.0205,0.0001,0.9786,0.0008])
yval22[0][:]=xbot
yimpv=yvalu(xbot,solverT(xbot,ant,2500).x[4],xbot).x[:4]
xc=xbot.copy()
Ti=solverT(xc,ant,2500).x[4]
print(Ti,"Ti")
Ta=Ti
R=1.5
maxp=5000
distillate=21.846
bottoms=feed-distillate
xw=xbot.copy()
yflow=np.zeros((1,4))
print(yvalu(xbot,600,yval22).x[4])
for i in range(0,4):
	yflow[0][i]=1e-4
ptd=yvalu(xc,Ta,yval22).x[4]
print(ptd,"initialpress")
q0=0
sumd=0
suml=0
sumb=0
sumv=0
sumf=0
def vaporTsolver(xt,yval,act,P):
	satu_press=np.zeros(5)
	for i in range(0,4):
		satu_press[i]=yval[i]-(ant[i][0]*np.exp((ant[i][1]/8.314)*((1/ant[i][2])-(1/xt[4])))*xt[i]*act[i]/P)
	satu_press[4]=1-xt[0]-xt[1]-xt[2]-xt[3]
	return satu_press
def solverTv(yval,P):
	yoyo=np.zeros(4)
	xass=xbot.copy()
	xassp=np.zeros(4)
	Tass=450
	act=get_activity_coeff(xass,Tass)
	Tassu=0
	sum=100
	while abs(sum)>1e-2:
		yoyo=least_squares(vaporTsolver,x0=np.array([0.59,0.001,0.39,0.02,450]),bounds=(np.array([1e-4,1e-4,1e-4,1e-4,337]),np.array([0.9999,0.9999,0.9999,0.9999,2000])),args=(yval,act,766)).x
		Tassu=Tass
		Tass=yoyo[4]
		xass[:]=yoyo[:4]
		for i in range(0,4):
			sum+=xass[i]-xassp[i]
		return yoyo
temptemp=solverTv(yimpv,766)[4]
for i in range(0,4):
	sumd+=xd[i]*(enth[i][0]+(c_p[i][0]*(345-303)))
	sumb+=xc[i]*(enth[i][0]+(c_p[i][0]*(Ti-303)))
	sumf+=xf[i]*(enth[i][0]+(c_p[i][0]*(feed_temp-303)))+(yimpv[i]*(enth[i][1]+(c_p[i][1]*(feed_temp-ant[i][2]))))
	suml+=xf[i]*(enth[i][0]+(c_p[i][0])*(solverT(xf,ant,766).x[4]-303))
	sumv+=yimpv[i]*(enth[i][1]+((c_p[i][1])*(temptemp-ant[i][2])))
qratio=(sumv-sumf)/(sumv-suml)
at=qratio
q0 = (sumd*distillate)+(sumb*bottoms)-(sumf*feed)
print(q0)
max=feed
min=bottoms*(s/(1-s))
yaaa=(reboilsolve(xc,yflow,0,Ta,bottoms,Ta,xf,0,q0,qratio,pi0,ptd,s))
xcp=xc.copy()
yval2=yflow.copy()
stages=0
flowval1=bottoms
xjp00=0
xjp0=xjp00
xcps=xc.copy()
Tapps=Ti
i=0
qc0=0
qc=0
q0p=0
q0h=1
qh=q0
#while abs((qh-q0h)/q0h)>0.3:  #abs((xw[1]-xbot[1])/xbot[1])>0.001 :
Rmin=0
alphad=np.zeros(4)
alphabot=np.zeros(4)
alphaavg=np.zeros(4)
alphaf=np.zeros(4)
for i in range(0,4):
	Rmin+=(xd[i]-yvalu(xf,feed_temp,yflow).x[i])/(yvalu(xf,feed_temp,yflow).x[i]-xf[i])
	alphabot[i]=(sat_pressure(ant[0][0],ant[0][1],ant[0][2],Ti)*get_activity_coeff(xbot,Ti)[0])/(sat_pressure(ant[i][0],ant[i][1],ant[i][2],Ti)*get_activity_coeff(xbot,Ti)[i])
	alphaf[i]=(sat_pressure(ant[0][0],ant[0][1],ant[0][2],feed_temp)*get_activity_coeff(xf,feed_temp)[0])/(sat_pressure(ant[i][0],ant[i][1],ant[i][2],feed_temp)*get_activity_coeff(xf,feed_temp)[i])
	alphad[i]=(sat_pressure(ant[0][0],ant[0][1],ant[0][2],345)*get_activity_coeff(xd,345)[0])/(sat_pressure(ant[i][0],ant[i][1],ant[i][2],345)*get_activity_coeff(xd,345)[i])
	alphaavg[i]=((alphad[i]*alphabot[i])**0.5)
print(alphaavg,"alphaavg")
def theta(the,alphaf):
	sum=0
	for i in range(0,4):
		sum+=(alphaf[i]*xf[i])/(alphaf[i]-the)
	return sum-(1-qratio)
def thetasolver(alphaavg):
	solbay=fsolve(theta,1.3,args=(alphaavg))
	return solbay[0]
def rminsolver(xd,the,alphaavg):
	sum=0
	for i in range(0,4):
		sum+=(alphaavg[i]*xd[i]/(alphaavg[i]-the))
	return sum-1
print(thetasolver(alphaavg))
print(rminsolver(xd,thetasolver(alphaavg),alphaavg),"rmin")
Rmin=Rmin**(0.5)
min=(1-qratio)*feed
max=bottoms*(s/(1-s))
Ta=Ti
stages=1
flowval=bottoms
tvalu=Ti
xjp0=0
ptd=2500
ptdp=ptd
yval1=np.zeros((1,4))
yval22=yflow.copy()
yval2=yval22.copy()
yval1[0][:]=yvalu(xbot,tvalu,yval22).x[0:4]
ptd=2500
pressure=Pressure 
print(pressure,"press")
print(yval1,"yval1",yval2,"yval2")
yay=reboilsolve(xbot,yval22,xjp0,tvalu,bottoms,Ta,xf,0,qh,qratio,pi0,ptd,s)
yval22=yval2.copy()
xcps=xbot.copy()
xcp=xbot.copy()
xc=yay[0:4]
print(xc,"xc",xcp,"xcp")
xjp00=xjp0
xjp0=yay[6]
flowval1=flowval
flowval=yay[7]
Taps=Ta
Tapp=Ta
Tap=tvalu
tvalu=yay[4]
qh=yay[5]
print(tvalu,"tvalu",xjp0,"xjp0",flowval,"flowval")
Ta=Tap
ptd0=ptd
pflw=32
Purity=0.05
solsol=yay
i=0
rectifying_stages=0
stripping_stages=1
print(qratio,"qratio")
sumt=0
while sumt<1e-5:
	min=0
	max=bottoms*(s)/(1-s)
	print(yval2,'yval2')
	print(yval1,"yval1")
	ptdp=ptd0
	ptd0=ptd
	ptd=yvalu(xc,tvalu,yval1).x[4]
	solsol=flowrates(xc,yval1,xjp0,tvalu,flowval,Ta,xf,0,0,pi0,ptd)
	yval22=yval2.copy()
	yval2=yval1.copy()
	yval1[0][:]=yvalu(xc,tvalu,yval1).x[0:4]
	xcps=xcp.copy()
	xcp=xc.copy()
	xc=solsol[0:4]
	xjp00=xjp0
	xjp0=solsol[5]
	print(xc,"xc")
	print(xcp,"xcp")
	print(xjp0,"xjp0")
	flowval1=flowval
	flowval=solsol[6]
	print(flowval1,"flowval1",flowval,"flowval")
	Tapps=Taps
	Taps=Ta
	Tap=tvalu
	tvalu=solsol[4]
	Ta=Tap
	print(Tap)
	print(tvalu)
	print("stripping")
	for i in range(0,4):
		sumt+=(yval1[0][i]-yval2[0][i])+(xc[i]-xd[i])
	sumt=(sumt+((tvalu-Ta)/485))/stripping_stages
	stages+=1
	stripping_stages+=1
	i+=1
ptd=ptdp
stripping_stages=stripping_stages-1
sum=0
while rectifying_stages <=0:
	print(yval2,'yval2')
	print(yval22,"yval22")
	if i<=0:
		solsol=reboilsolve(xcp,yval2,0,Ta,bottoms,Taps,xf,feed,0,qratio,pi0,ptd,s)
		qh=solsol[5]
		flowval2=flowval1
		flowval1=solsol[7]
		xjp000=xjp00
		xjp00=solsol[6]
	else:
		solsol=flowratesfeed(xcp,yval2,xjp00,Ta,flowval1,Taps,xf,feed,qratio,pi0,ptd)
		xjp000=xjp00
		xjp00=solsol[5]
		flowval2=flowval1
		flowval1=solsol[6]
	yval22=yval2.copy()
	ptdp=ptd
	ptd=yvalu(xcp,Ta,yval2).x[4]
	yval2[0][:]=yvalu(xcp,Ta,yval2).x[0:4]
	xcps=xcp.copy()
	xcpp=np.vstack([xcp,np.ones(len(xcp))])
	xcp=solsol[0:4]
	print(xc,"xc")
	print(xcp,"xcp")
	print(xcps,"xcps")
	print(xjp00,"xjp0")
	Tapps=Taps
	Taps=Ta
	Ta=solsol[4]
	print(Taps)
	print(Ta)
	print(yval2,"yval new")
	print("rectifyingwfeed")
	rectifying_stages+=1
	stages+=1
	for i in range(0,4):
			sum += (yval2[0][i]-yval22[0][i]) + (xcp[i]-xcps[i])
			print(sum)
	sum=(sum+((Ta-Taps)/345))/stages
	print(sum,"sum difference")
sum=0
print(stages,"stages")
print(xjp0,"xjp0")
print(xjp00,"xjp00")
while abs(sum)<=1e-1 and ((xjp000-flowval2)/xjp000)>0.001: #and xcp[1]<xd[1] and xcp[3]<xd[3]:# and yval2[0][3]>yval22[0][3]: #xcp[3]<xd[3] and xcp[1]<xd[1]:
	min=distillate
	max=xjp00
	print(yvalu(xd,345,xbot).x,"345")
	print(yval2,'yval2')
	print(yval22,"yval22")
	solsol=flowrates(xcp,yval2,xjp00,Ta,flowval1,Taps,xf,0,0,pi0,ptdp)
	ptdp=ptd
	ss=yvalu(xcp,Ta,xcp)
	yy=yvalu(xd,solsol[4],xd)
	ptd=ss.x[4]
	yval22=yval2.copy()
	yval2[0][:]=ss.x[0:4]
	xcps=xcp
	xcp=np.array(solsol[0:4])
	xcpp=np.vstack([xcp,np.ones(len(xcp))])
	xjp000=xjp00
	xjp00=solsol[5]
	print(xc,"xc")
	print(xcp,"xcp")
	print(xcps,"xcps")
	print(xjp000,"xjp000",xjp00,"xjp0")
	flowval2=flowval1
	flowval1=solsol[6]
	print(flowval2,"flowval2",flowval1,"flowval1")
	Tapps=Taps
	Taps=Ta
	Ta=solsol[4]
	print(Taps)
	print(Ta)
	print("rectifying")
	rectifying_stages+=1
	stages+=1
	for i in range(0,4):
			sum+=(yval2[0][i]-yval22[0][i])+(xcp[i]-xcps[i])#solvex(yval2,Ta)[i]-xd[i]
			print(sum)
	sum=(sum+((Ta-Taps)/345))/stages
	print(sum,"sumsum")
print(Rmin,"Rmin")
print(s,"s")
print(stripping_stages,"ss")
print(rectifying_stages,"rs")
print(xjp000,"xjp000")
print(flowval2,"flowval2")
print(solvex(yval2,Ta),"SOLVEX")
sumvap=0
sumdis=0
for i in range(0,4):
	sumvap+=yvalu(xcp,Ta,yval2).x[i]*(enth[i][1]+(c_p[i][1]*(Taps-ant[i][2])))
	sumdis+=xcp[i]*(enth[i][0]+(c_p[i][0]*(Ta-303)))
qc0=qc
qc=qh-q0
print(qc)
q0p=q0
stages= stripping_stages + rectifying_stages
def xvals(xbott,xfeed,xtops):
	return np.array([1-xbott[0]-xbott[1]-xbott[2]-xbott[3],xbott[0]-(((xfeed[0]*feed)-(xtops[0]*distillate))/bottoms),xbott[1]-(((xfeed[1]*feed)-xtops[1]*distillate)/bottoms),xbott[2]-(((xfeed[2]*feed)-(xtops[2]*distillate))/bottoms)])
def xsolve(xfeed,xtops):
	tt=xbot
	soln=least_squares(xvals,x0=xbot,bounds=(np.array([1e-4,1e-4,1e-4,1e-4]),np.array([0.999,0.999,0.999,0.999])),args=(xfeed,xtops))
	return soln
sumbb=0
print(solvex(yval2,Ta),"xvals")
print(xjp0,"xjp0")
#for i in range(0,4):
	#sumbb+=xbot[i]*(c_p[i][0]*(Ti-303))
#q0=(sumd*distillate)-(sumf*feed)+(sumbb*bottoms)
#q0h=qh
#qh=q0+qc
#xw=xsolve(xf,xcp).x[0:4]
print(xcp,stages, " " + "+" + " " + "condensor",flowval1,"flowval1",xjp00,"xjp00")
#print(xvalu(yval2,Ta,ptd),"xvalu")
y00=yvalu(xcps,Taps,yval2)
print(yval2,"yval2")
print(yval22,"yval22")
print(y00.x,"yval1")
print(yvalu(xcp,Ta,yval2).x[4],"presstops")
ptd=2500
print(yvalu(xbot,Ti,yval2).x[4],"press")
print(xcps,"xcps")
print(xbot,"xbot")
print(xw,"xw")
print(qratio,"qratio")
print(qh,"qh",qc,"qc")
print(yvalu(xcp,Ta,xbot).x,"yvalnext")
#while xjp000 >0 and abs((xjp000-(2.5*flowval2))/xjp000)>0.0001:
#(q/Ta)+((xjp1[6])*((xjp1[0]*((0.452+(c_p[0][0]*np.log(xjp1[4]/303.0)))))+(xjp1[1])*(0.126+c_p[1][0]*np.log(xjp1[4]/303.0))+(xjp1[2])*(0.495+(c_p[2][0]*np.log(xjp1[4]/303.0)))+(xjp1[3])*(0.069+(c_p[3][0]*np.log(xjp1[4]/303.0)))-sump))-((flow)*((xb[0]*((0.452+(c_p[0][0]*np.log(Ta/303.0)))))+(xb[1])*(0.126+c_p[1][0]*np.log(Ta/303.0))+(xb[2])*(0.495+(c_p[2][0]*np.log(Ta/303.0)))+(xb[3])*(0.069+(c_p[3][0]*np.log(Ta/303.0)))-sumgb))-((xjp1[5])*((yval[0][0]*((0.452+c_p[0][1]*(np.log(Ta/303.0))+(0.00831*np.log(pi0[0][0]/ptu))))+(yval[0][1]*(0.239+c_p[1][1]*np.log(Ta/303.0))+(0.00831*np.log(pi0[0][1]/ptu)))+(yval[0][2]*(0.288+c_p[2][1]*np.log(Ta/303.0)+(0.00831*np.log(pi0[0][2]/ptu))))+(yval[0][3]*(0.188+c_p[3][1]*np.log(Ta/303.0)+(0.00831*np.log(pi0[0][3]/ptu)))))-sump))+((xjp0)*((yval0[0][0]*((0.452+c_p[0][1]*(np.log(tval/303.0)))+(0.00831*np.log(pi0[0][0]/ptd)))+(yval0[0][1]*(0.239+c_p[1][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][1]/ptd))))+(yval0[0][2]*(0.288+c_p[2][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][2]/ptd))))+(yval0[0][3]*(0.188+c_p[3][1]*np.log(tval/303.0)+(0.00831*np.log(pi0[0][3]/ptd)))))-sumd))+((feed)*((xf[0]*((0.452+(c_p[0][0]*np.log(400/303.0)))))+(xf[1])*(0.126+c_p[1][0]*np.log(400/303.0))+(xf[2])*(0.495+(c_p[2][0]*np.log(400/303.0)))+(xf[3])*(0.069+(c_p[3][0]*np.log(400/303.0)))-sumgf))
print(solvex(yval22,345),"xval")
print(sumf,"sumf")
print(sumv,"sumv")