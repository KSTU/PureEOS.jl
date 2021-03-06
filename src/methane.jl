"""
Calculate methane critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
"""
function methaneCritical()
	sub=newCritical()
	sub.T=190.564
	sub.p=4.5922*10.0	#MPa to bar
	sub.ro=162.66/16.04	#mol/l #322 kg/m^3
	return sub
end

"""
Calculate water properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function methaneTriple()
	sub=newTriple()
	sub.T=90.6941
	sub.p=116.96*100*10.0^-5	#Pa to bar
	sub.lRo=-1.0	#kg/m^3 to mol/l
	sub.vRo=-1.0	#kg/m^3s to mol/l
	return sub
end

"""
calculate methane properties at saturation line
source
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function methaneSaturate(T)
a=[-6.036219, 1.409353, -0.4945199, -1.443048]
b=[1.9906389, -0.78756197, 0.036976723]
c=[-1.8802840, -2.8526531, -3.0006480, -5.2511690, -13.191859, -37.553961]

sub=newSaturate()
crit=methaneCritical()
if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=crit.p*exp(crit.T/T*(a[1]*v+a[2]*v^1.5+a[3]*v^2.0+a[4]*v^4.5))
	sub.lRo=crit.ro*exp(b[1]*v^0.354+b[2]*v^0.5+b[3]*v^(5.0/2.0))
	sub.vRo=crit.ro*exp(c[1]*v^0.354+c[2]*v^(5.0/6.0)+c[3]*v^(3.0/2.0)+c[4]*v^(5.0/2.0)+c[5]*v^(25.0/6.0)+c[6]*v^(47.0/6.0))
	sub.T=T
else
	sub.p=-1
	sub.lRo=-1
	sub.vRo=-1
	sub.T=-1
end
return sub
end


"""
calculate methane properties
source
T --- temperature, Kelvin
p --- pressure, bar
ro --- density, mol/l
z ---
s --- entropy,
s0 --- ideal part of entropy,
ds --- nonideal part derived k_B, dimensionless
g --- Gibbs energy
g0 --- ideal part of Gibbs energy
dg --- nonideal part of Gibbs energy
f ---  Helmholtz free energy
f0 --- ideal part
df --- nonideal part, dimensionless
u --- internal energy
du --- nonideal part
h --- entalpy
dh --- nonideal part
w --- speed of sound,
cp ---
cv ---
"""
function methaneSinglephase(T,ro)

sub=newSingle()
crit=methaneCritical()
delta=ro/crit.ro
tau=crit.T/T
R=8.3144598
mm=0.01604

nOne=13
nTwo=23
nThree=4
nFour=0


c=[0.00, 0.00, 0.00, 0.00, 0.00,	#corrected
   0.00, 0.00, 0.00, 0.00, 0.00,	#10
   0.00, 0.00, 0.00, 1.00, 1.00,
   1.00, 1.00, 1.00, 1.00, 1.00,	#20
   2.00, 2.00, 2.00, 2.00, 2.00,
   3.00, 3.00, 3.00, 3.00, 4.00,	#30
   4.00, 4.00, 4.00, 4.00, 4.00,
   4.00 ]

d=[1.0, 1.0, 1.0, 2.0, 2.0,	#corrected
   2.0, 2.0, 3.0, 4.0, 4.0,	#10
   8.0, 9.0, 10.0, 1.0, 1.0,
   1.0, 2.0, 4.0, 5.0, 6.0,	#20
   1.0, 2.0, 3.0, 4.0, 4.0,
   3.0, 5.0, 5.0, 8.0, 2.0,	#30 
   3.0, 4.0, 4.0, 4.0, 5.0,
   6.0, 2.0, 0.0, 0.0, 0.0 ]	#40

t=[-0.5, 0.5, 1.0, 0.5, 1.0,
1.5, 4.5, 0.0, 1.0, 3.0,	#10
1.0, 3.0, 3.0, 0.0, 1.0,
2.0, 0.0, 0.0, 2.0, 2.0,	#20
5.0, 5.0, 5.0, 2.0, 4.0,
12.0, 8.0, 10.0, 10.0, 10.0,	#30
14.0, 12.0, 18.0, 22.0, 18.0,
14.0, 2.0, 0.0, 1.0, 2.0]

n=[0.4367901028e-01, 0.6709236199,	#corrected
-0.1765577859e01, 0.8582330241,
-0.1206513052e01, 0.5120467220,
-0.4000010791e-03, -0.1247842423e-01,
0.3100269701e-01, 0.1754748522e-02,
-0.3171921605e-05, -0.2240346840e-05,
0.2947056156e-06, 0.1830487909,
0.1511883679, -0.4289363877,
0.6894002446e-01, -0.1408313996e-01,
-0.3063054830e-01, -0.2969906708e-01,
-0.1932040831e-01, -0.1105739959,
0.9952548995e-01, 0.8548437825e-02,
-0.6150555662e-01, -0.4291792423e-01,
-0.1813207290e-01, 0.3445904760e-01,
-0.2385919450e-02, -0.1159094939e-01,
0.6641693602e-01, -0.2371549590e-01,
-0.3961624905e-01, -0.1387292044e-01,
0.3389489599e-01, -0.2927378753e-02,
0.9324799946e-04, -0.6287171518e01,
0.1271069467e02, -0.6423953466e01]

al=[0.0, 0.0, 0.0, 0.0, 0.0,	#corrected
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 20.0, 40.0, 40.0, 40.0]	#40

bt=[0.0, 0.0, 0.0, 0.0, 0.0,	#corrected
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 200.0, 250.0, 250.0, 250.0]	#40

gm=[0.0, 0.0, 0.0, 0.0, 0.0,	#corrected
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 1.07, 1.11, 1.11, 1.11]	#40

ep=[0.0, 0.0, 0.0, 0.0, 0.0,	#corrected delta in article
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 1.0, 1.0, 1.0, 1.0]	#40


phidelta=0.0
phir=0.0
phitau=0.0
phidelta2=0.0
phitau2=0.0
dphideltadtau=0.0

i=0
for ni=1:nOne
	i=i+1
	phidelta=phidelta+n[i]*d[i]*delta^(d[i]-1.0)*tau^t[i]
	phir=phir+n[i]*delta^d[i]*tau^t[i]
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1)
	phidelta2=phidelta2+n[i]*d[i]*(d[i]-1.0)*delta^(d[i]-2)*tau^t[i]
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2.0)
	dphideltadtau=dphideltadtau+n[i]*d[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)
end
for ni=1:nTwo
	i=i+1
	phidelta=phidelta+n[i]*exp(-delta^c[i])*(delta^(d[i]-1.0)*tau^t[i]*(d[i]-c[i]*delta^c[i]))
	phir=phir+n[i]*delta^(d[i])*tau^t[i]*exp(-delta^c[i])
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1.0)*exp(-delta^c[i]);
	phidelta2=phidelta2+n[i]*exp(-delta^c[i])*(delta^(d[i]-2)*tau^(t[i])*((d[i]-c[i]*delta^c[i])*(d[i]-1.0-c[i]*delta^c[i])-c[i]^2*delta^c[i]))
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2)*exp(-delta^c[i])
	dphideltadtau=dphideltadtau+n[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)*(d[i]-c[i]*delta^c[i])*exp(-delta^c[i])
end
for pn=1:nThree
	i=i+1
	phidelta=phidelta+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta-2.0*al[i]*(delta-ep[i]))
	phir=phir+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)
	phitau=phitau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(t[i]/tau-2*bt[i]*(tau-gm[i]))
	phidelta2=phidelta2+n[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(-2*al[i]*delta^d[i]+4*al[i]^2*delta^d[i]*(delta-ep[i])^2-4*d[i]*al[i]*delta^(d[i]-1)*(delta-ep[i])+d[i]*(d[i]-1)*delta^(d[i]-2))
	phitau2=phitau2+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*((t[i]/tau-2*bt[i]*(tau-gm[i]))^2-t[i]/tau^2-2*bt[i])
	dphideltadtau=dphideltadtau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta -2*al[i]*(delta-ep[i]))*(t[i]/tau-2*bt[i]*(tau-gm[i]))
end

for pn=1:nFour
	i=i+1
	theta=(1.0-tau)+A[i]*((delta-1.0)^2)^(1.0/(2.0*bt[i]))
	delb=(theta)^2+B[i]*((delta-1.0)^2)^a[i]
	psif=exp(-C[i]*(delta-1.0)^2-D[i]*(tau-1.0)^2)

	dpsidd=-2.0*C[i]*(delta-1.0)*psif
	dpsidd2=(2*C[i]*(delta-1)^2-1)*2*C[i]*psif
	dpsidtau=-2*D[i]*(tau-1)*psif
	dpsidtau2=(2*D[i]*(tau-1)^2-1)*2*D[i]*psif
	dpsidtaudd=4*C[i]*D[i]*(delta-1)*(tau-1)*psif

	dddd=(delta-1.0)*(A[i]*theta*2/bt[i]*((delta-1)^2)^(1/2/bt[i]-1)+2*B[i]*a[i]*((delta-1)^2)^(a[i]-1))

	if ((delta==1) && (tau==1))
		ddbdd=0
		ddbdtau=0
		ddbdd2=0
		ddbdtau2=0
		ddbdddtau=0
	else
		dddd2=1/(delta-1)*dddd+(delta-1)^2*(4*B[i]*a[i]*(a[i]-1)*((delta-1)^2)^(a[i]-2)+2*A[i]^2*(1/bt[i])^2*(((delta-1)^2)^(1/2/bt[i]-1))^2+A[i]*theta*4/bt[i]*(1/2/bt[i]-1)*((delta-1)^2)^(1/2/bt[i]-2))

		ddbdd=b[i]*delb^(b[i]-1)*dddd
		ddbdd2=b[i]*(delb^(b[i]-1)*dddd2+(b[i]-1)*delb^(b[i]-2)*(dddd)^2)
		ddbdtau=-2*theta*b[i]*delb^(b[i]-1)
		ddbdtau2=2*b[i]*delb^(b[i]-1)+4*theta^2*b[i]*(b[i]-1)*delb^(b[i]-2)
		ddbdddtau=-A[i]*b[i]*2/bt[i]*delb^(b[i]-1)*(delta-1)*((delta-1)^2)^(1/2/bt[i]-1)-2*theta*b[i]*(b[i]-1)*delb^(b[i]-2)*dddd
	end
	#
	phidelta=phidelta+n[i]*(delb^b[i]*(psif+delta*dpsidd)+ddbdd*delta*psif);
	phir=phir+n[i]*delb^b[i]*delta*psif
	phitau=phitau+n[i]*delta*(ddbdtau*psif+delb^b[i]*dpsidtau)
	phidelta2=phidelta2+n[i]*(delb^b[i]*(2*dpsidd+delta*dpsidd2)+2*ddbdd*(psif+delta*dpsidd)+ddbdd2*delta*psif)
	phitau2=phitau2+n[i]*delta*(ddbdtau2*psif+2*ddbdtau*dpsidtau+delb^b[i]*dpsidtau2)
	dphideltadtau=dphideltadtau+n[i]*(delb^b[i]*(dpsidtau+delta*dpsidtaudd)+delta*ddbdd*dpsidtau+ddbdtau*(psif+delta*dpsidd)+ddbdddtau*delta*psif)
end

sub.T=T
sub.ro=ro
sub.p=(1.0+phidelta*delta)*ro*R*T/100.0	#bar
sub.z=1.0+phidelta*delta	#
#ideal part
n0=[9.91243972, -6.33270087, 3.0016, 0.008449, 4.6942, 3.4865, 1.6572, 1.4115]
gamma0=[0.0, 0.0, 0.0, 3.40043240, 10.26951575, 20.43932747, 29.93744884, 79.13351945]

#entropy
phi0=log(delta)+n0[1]+n0[2]*tau+n0[3]*log(tau)+n0[4]*log(1.0-exp(-gamma0[4]*tau))+n0[5]*log(1.0-exp(-gamma0[5]*tau))+n0[6]*log(1.0-exp(-gamma0[6]*tau))+n0[7]*log(1.0-exp(-gamma0[7]*tau))+n0[8]*log(1.0-exp(-gamma0[8]*tau))
phi0tau=n0[2]+n0[3]/tau+n0[4]*gamma0[4]*(1.0/(1.0-exp(-gamma0[4]*tau))-1.0)+n0[5]*gamma0[5]*(1.0/(1.0-exp(-gamma0[5]*tau))-1.0)+n0[6]*gamma0[6]*(1.0/(1.0-exp(-gamma0[6]*tau))-1.0)+n0[7]*gamma0[7]*(1.0/(1.0-exp(-gamma0[7]*tau))-1.0)+n0[8]*gamma0[8]*(1.0/(1.0-exp(-gamma0[8]*tau))-1.0)
phi0delta=1/delta
phi0delta2=-1/delta^2
phi0tau2=-n0[3]/tau^2-n0[4]*gamma0[4]^2*exp(-gamma0[4]*tau)*(1.0-exp(-gamma0[4]*tau))^(-2)-n0[5]*gamma0[5]^2*exp(-gamma0[5]*tau)*(1.0-exp(-gamma0[5]*tau))^(-2)-n0[6]*gamma0[6]^2*exp(-gamma0[6]*tau)*(1.0-exp(-gamma0[6]*tau))^(-2)-n0[7]*gamma0[7]^2*exp(-gamma0[7]*tau)*(1.0-exp(-gamma0[7]*tau))^(-2)-n0[8]*gamma0[8]^2*exp(-gamma0[8]*tau)*(1.0-exp(-gamma0[8]*tau))^(-2)
phi0taudelta=0

#Entropy
sub.s=(tau*(phi0tau+phitau)-phi0-phir)*R
sub.s0=(tau*(phi0tau)-phi0)*R
sub.ds=(tau*(phitau)-phir)	#devidet to R

#out.entrfrac=(tau*(phi0tau+phitau)-phi0-phir)/(tau*(phi0tau)-phi0)

dtaudtau=-tau/T
dphidt=dtaudtau*dphideltadtau

#sub.dpdt=(delta*cur_temp*dphidt+delta*phidelta+1)*R*T/100;
#phi0
#phir
#phidelta
#phidelta*delta
#(1+phi0+phir+phidelta*delta)
sub.g=(1+phi0+phir+phidelta*delta)*R*T	#R*cur_temp*
sub.g0=(1+phi0)*R*T
sub.dg=(phir+phidelta*delta)	#deriver RT
#
sub.f=(phi0+phir)*R*T
sub.f0=phi0*R*T
sub.df=phir

#Energy
sub.u=tau*(phi0tau+phitau)*R*T
sub.du=tau*(phitau)
#
sub.h=(1+tau*(phi0tau+phitau)+delta*phidelta)*R*T
sub.dh=1+tau*(phitau)+delta*phidelta;
#speed of sound
sub.w=sqrt(R/mm*T*(1+2*delta*phidelta+delta^2*phidelta2-(1+delta*phidelta-delta*tau*dphideltadtau)^2/(tau^2*(phi0tau2+phitau2))))	#0.018 перевод в килограмы
#
sub.cv=-tau^2*(phi0tau2+phitau2)*R
sub.cp=(-tau^2*(phi0tau2+phitau2)+(1+delta*phidelta-delta*tau*dphideltadtau)^2/(1+2*delta*phidelta+delta^2*phidelta2))*R
sub.dpdt=R/100.0*T*(1.0+2*delta*phidelta+delta^2*phidelta2)
#debug
# println("phi0 $(phi0)")
# println("phi0delta $(phi0delta)")
# println("phi0delta2 $(phi0delta2)")
# println("phi0tau $(phi0tau)")
# println("phi0tau2 $(phi0tau2)")
# println("phi0taudelta $(phi0taudelta)")
#
# println("phir $(phir)")
# println("phidelta $(phidelta)")
# println("phidelta2 $(phidelta2)")
# println("phitau $(phitau)")
# println("phitau2 $(phitau2)")
# println("dphideltadtau $(dphideltadtau)")
#enddebug
return sub
end


"""
calculate water viscosity
Huber M.L. et al. New International Formulation for the Viscosity of H[sub 2]O // J. Phys. Chem. Ref. Data. 2009. Vol. 38, № 2. P. 101.
10.1063/1.3088050


"""
function methaneVisc(T,n)
mu_ref=1	#*10^(-6)	#\mu Pa*s
Tb=T/waterCritical().T;
nb=n/waterCritical().ro;

H0=[1.67752, 2.20462, 0.6366564, -0.241605]

mu_0=100*sqrt(Tb)
sum_mu_0=0;
for ii=1:4
	sum_mu_0=sum_mu_0+H0[ii]/Tb^(ii-1);
end
mu_0=mu_0/sum_mu_0;

H1=[[5.20094e-1, 8.50895e-2, -1.08374, -2.89555e-1 , 0, 0] [2.22531e-1, 9.99115e-1, 1.88797, 1.26613, 0, 1.20573e-1] [-2.81378e-1, -9.06851e-1, -7.72479e-1, -4.89837e-1, -2.57040e-1,0] [1.61913e-1, 2.57399e-1,0, 0, 0, 0] [-3.25372e-2, 0, 0, 6.98452e-2, 0, 0] [0, 0, 0, 0, 8.72102e-3, 0] [0, 0, 0, -4.35673e-3, 0, -5.93264e-4]]

mu_1=0
sum_ii=0
for ii=1:6
	sum_jj=0
	for jj=1:7
		sum_jj=sum_jj+H1[ii,jj]*(nb-1)^(jj-1)
	end
	sum_ii=sum_ii+(1/Tb-1)^(ii-1)*sum_jj
end
mu_1=exp(nb*sum_ii)

ksim=0.068
qC=1.9
qD=1.1
vv=0.63
gam=1.239
dz0=0.13
G0=0.06
TR=1.5
mu_2=1
return mu_ref*mu_0*mu_1*mu_2
end


function methaneViscMax(T)
	sub=newViscMax()
	if(T<=373.0)
		sub.p=10000.0
		init=65.0
	elseif((T>373.0)&&(T<=433.0))
		sub.p=5000.0
		init=60.0
	elseif((T>433.0)&&(T<=873.0))
		sub.p=3500.0
		init=57.0
	else
		sub.p=3000.0
		init=52.0
	end
	sub.ro=fzero(x -> waterSinglephase(T,x).p-sub.p, init)
	return sub
end

"""
thermalconductivity of water
*warning* due using IAPWS-95 wrong prediction in near critical region

"""
function methaneTherm(T,n)
	crit=waterCritical()
	Tb=T/crit.T
	nb=n/crit.ro
	lab=1	#0.001
	R=8.3144598

	L0=[2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]

	la0=sqrt(Tb)/(L0[1]+L0[2]/Tb+L0[3]/Tb^2+L0[4]/Tb^3+L0[5]/Tb^4)

	L1=[[1.60397357, 2.33771842, 2.19650529, -1.21051378, -2.7203370] [-0.646013523, -2.78843778, -4.54580785, 1.60812989, 4.57586331] [0.111443906, 1.53616167, 3.55777244, -0.621178141, -3.18369245] [0.102997357, -0.463045512, -1.40944978, 0.0716373224, 1.1168348] [-0.0504123634, 0.0832827019, 0.275418278, 0, -0.19268305] [0.00609859258, -0.00719201245, -0.0205938816, 0, 0.012913842]]

	temp1=0.0

	for ii=1:5
		temp2=0;
		for jj=1:6
			temp2=temp2+L1[ii,jj]*(nb-1)^(jj-1);
		end
		temp1=temp1+temp2*(1/Tb-1)^(ii-1);
	end
	la1=exp(nb*temp1);

	LA=177.8514
	qD=1.0/0.4
	vv=0.63
	gam=1.239
	dz0=0.13
	g0=0.06
	Tr=1.5
#
	# if(nb<=0.310559006)\
	# 	A=[6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709, 1.97815050331519]
	# elseif((nb>=0.310559006) && (nb<0.776397516))
	# 	A=[6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395, 5.54349664571295]
	# elseif((nb>=0.776397516) && (nb<1.242236025))
	# 	A=[5.35500529896124, -3.96415689925446, 8.91990208918795, -12.0338729505790, 9.19494865194302, -2.16866274479712]
	# elseif((nb>=1.242236025) && (nb<1.863354037))
	# 	A=[1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.16780999933360, -0.965458722086812]
	# else
	# 	A=[1.1199926419994, 0.595748562571649, 9.88952565078920, -10.3255051147040, 4.66861294457414, -0.503243546373828 ]
	# end
    #
	# tempsum=0
	# for ii=1:6
	# 	tempsum+=A[ii]*nb^(ii-1)
	# end
#
	cur=waterSinglephase(T,n)

	grt=1/cur.dpdt*crit.p/crit.ro
	grtr=1/waterSinglephase(Tr*waterCritical().T,n).dpdt*crit.p/crit.ro
	dksi=nb*(grt-grtr*Tr/Tb)
	k=cur.cp/cur.cv
	if (dksi>0)
		dz=dz0*(dksi/g0)^(vv/gam)
	else
		dz=dz0
	end
	y=qD*dz
	if (y<1.7e-7)
		Z=0
	else
		Z=2/pi/y*(((1-1/k)*atan(y)+y/k)-(1.0-exp(-1/(1/y+y^2/3/nb^2))))
	end
	mub=waterVisc(T,n)

	la2=LA*nb*waterSinglephase(T,n).cp/R*Tb/mub*Z
	# println("$(la1) - $(la2)  ")
	return lab*(la0*la1+la2)
end
