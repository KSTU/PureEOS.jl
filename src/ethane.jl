"""
Calculate ethane critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
"""
function ethaneCritical()
	sub=newCritical()
	sub.T=305.322
	sub.p=4.8722*10.0	#MPa to bar
	sub.ro=206.18/30.06904	#mol/l #322 kg/m^3
	return sub
end

"""
Calculate ethane properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function ethaneTriple()
	sub=newTriple()
	sub.T=90.368
	sub.p=1.14*10^-5	#Pa to bar
#	sub.lRo=999.793/18.015268	#kg/m^3 to mol/l
#	sub.vRo=0.00485458/18.015268	#kg/m^3s to mol/l
	return sub
end

"""
calculate water properties at saturation line
Wagner W. The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use // J. Phys. Chem. Ref. Data. 1999. Vol. 31, № 2. P. 387.
10.1063/1.1461829
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function ethaneSaturate(T)
a=[-6.48647577, 1.47010078, -1.66261122, 3.57898378, -4.79105705]
b=[1.56138026, -0.381552776, 0.0785372040, 0.0370315089]
c=[-1.89879145, -3.65459262, 0.850562745, 0.363965487, -1.50005943, -2.26690389]
sub=newSaturate()
crit=ethaneCritical()
if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=crit.p*exp(crit.T/T*(a[1]*v+a[2]*v^1.5+a[3]*v^2.5+a[4]*v^3.5+a[5]*v^4.0))
	sub.lRo=crit.ro*exp(b[1]*v^0.329+b[2]*v^(4.0/6.0)+b[3]*v^(8.0/6.0)+b[4]*v^(19.0/6.0))
	sub.vRo=crit.ro*exp(crit.T/T*(c[1]*v^(0.346)+c[2]*v^(5.0/6.0)+c[3]*v+c[4]*v^2.0+c[5]*v^3.0+c[6]*v^5))
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
calculate ethane properties
A Reference Equation of State for the Thermodynamic Properties
of Ethane for Temperatures from the Melting Line to 675 K
and Pressures up to 900 MPa
D. Bücker  and W. Wagner 
10.1063/1.1859286
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
function ethaneSinglephase(T,ro)

sub=newSingle()
crit=ethaneCritical()
delta=ro/crit.ro
tau=crit.T/T
R=8.3144598
mm=0.03006904

nOne=5
nTwo=34
nThree=5


c=[0.00, 0.00, 0.00, 0.00, 0.00,
   1.00, 1.00, 1.00, 1.00, 1.00,	#10
   1.00, 1.00, 1.00, 1.00, 1.00,
   2.00, 2.00, 2.00, 2.00, 2.00,	#20
   2.00, 2.00, 2.00, 3.00, 3.00,
   3.00, 3.00, 3.00, 3.00, 3.00,	#30
   3.00, 3.00, 3.00, 3.00, 4.00,
   4.00, 4.00, 4.00, 4.00, 0.00,	#40
   0.00, 0.00, 0.00, 0.00]

d=[1.00, 1.00, 2.00, 2.00, 4.00,
   1.00, 1.00, 2.00, 2.00, 3.00,	#10
   6.00, 6.00, 7.00, 9.00, 10.00,
   2.00, 4.00, 4.00, 5.00, 5.00,	#20
   6.00, 8.00, 9.00, 2.00, 3.00,
   3.00, 3.00, 4.00, 4.00, 5.00,	#30
   5.00, 6.00, 11.00, 14.00, 3.00,
   3.00, 4.00, 8.00, 10.00, 1.00,	#40
   1.00, 3.00, 3.00, 2.00 ]

t=[ 0.25, 1.00, 0.25, 0.75, 0.75,
	2.00, 4.25, 0.75, 2.25, 3.00,	#10
	1.00, 1.25, 2.75, 1.00, 2.00,
	2.50, 5.50, 7.00, 0.50, 5.50,	#20
	2.50, 4.00, 2.00, 10.0, 16.0,
	18.0, 20.0, 14.0, 18.0, 12.0,	#30
	19.0, 7.00, 15.0, 9.00, 26.0,
	28.0, 28.0, 22.0, 13.0, 0.00,	#40
	3.00, 3.00, 0.00, 3.00]

n=[0.83440745735241, -0.14287360607171e1,
0.34430242210927, -0.42096677920265,
0.12094500886549e-1, -0.57976201597341,
-0.33127037870838e-1, -0.11751654894130,
-0.11160957833067, 0.62181592654406e-1,
0.98481795434443e-1, -0.98268582682358e-1,
-0.23977831007049e-3, 0.69885663328821e-3,
0.19665987803305e-4, -0.14586152207928e-1,
0.46354100536781e-1, 0.60764622180645e-2,
-0.26447330147828e-2, -0.42931872689904e-1,
0.29987786517263e-2, 0.52919335175010e-2,
-0.10383897798198e-2, -0.54260348214694e-1,
-0.21959362918493, 0.35362456650354,
-0.12477390173714, 0.18425693591517,
-0.16192256436754, -0.82770876149064e-1,
0.50160758096437e-1, 0.93614326336655e-2,
-0.27839186242864e-3, 0.23560274071481e-4,
0.39238329738527e-2, -0.76488325813618e-3,
-0.49944304440730e-2, 0.18593386407186e-2,
-0.61404353331199e-3,-0.23312179367924e-2,
0.29301047908760e-2, -0.26912472842883e-3,
0.18413834111814e3, -0.10397127984854e2]

al=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 15.0,	#40
	15.0, 15.0, 20.0, 20.0]

bt=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 150.0,	#40
	150.0, 150.0, 275.0, 400.0]

gm=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 1.05,	#40
	1.05, 1.05, 1.22, 1.16]

ep=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#20
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#30
	0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 1.0,	#40
	1.0, 1.0, 1.0, 1.0]

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


sub.T=T
sub.ro=ro
sub.p=(1.0+phidelta*delta)*ro*R*T/100.0	#bar
sub.z=1.0+phidelta*delta	#
#ideal part
n0=[9.212802589, 
-4.682248550,
3.003039265,
1.117433359,
3.467773215,
6.941944640,
5.970850948,
0.0]
gamma0=[0.0, 0.0, 0.0, 1.4091052332, 4.0099170712, 6.5967098342, 13.9798102659,1.0]

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


