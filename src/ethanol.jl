"""
Calculate ethanol critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
"""
function ethanolCritical()
	sub=newCritical()
	sub.T=514.71
	sub.p=6.268*10	#MPa to bar
	sub.ro=5.93		#mol/l #322 kg/m^3
	return sub
end

"""
Calculate ethanol properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function ethanolTriple()
	sub=newTriple()
	sub.T=159.0
	sub.p=7.185*10^(-9)	#MPa to bar
	sub.lRo=19.731	#mol/l
	sub.vRo=5.43*10*(-10)	#mol/l
	return sub
end


"""
calculate ethanol properties at saturation line
A Fundamental Equation of State for Ethanol
J. A. Schroeder, S. G. Penoncello, and J. S. Schroeder
Citation: Journal of Physical and Chemical Reference Data 43, 043102 (2014); doi: 10.1063/1.4895394
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function ethanolSaturate(T)
N=[9.00921, -23.1668, 30.9092, -16.5459, 3.64294]
k=[0.5, 0.8, 1.1, 1.5, 3.3]
Nv=[-1.75362, -10.5323, -37.6407, -129.762]
kv=[0.21, 1.1, 3.4, 10.0]
Np=[-8.94161, 1.61761, -51.1428, 53.1360]
kp=[1.0, 1.5, 3.4, 3.7]

sub=newSaturate()
crit=ethanolCritical()
if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=crit.p*exp(crit.T/T*(Np[1]*v^kp[1]+Np[2]*v^kp[2]+Np[3]*v^kp[3]+Np[4]*v^kp[4]))
	sub.lRo=crit.ro*(1.0+N[1]*v^k[1]+N[2]*v^k[2]+N[3]*v^k[3]+N[4]*v^k[4]+N[5]*v^k[5])
	sub.vRo=crit.ro*exp(Nv[1]*v^kv[1]+Nv[2]*v^kv[2]+Nv[3]*v^kv[3]+Nv[4]*v^kv[4])
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
calculate ethanol properties at single phase
A Fundamental Equation of State for Ethanol
J. A. Schroeder, S. G. Penoncello, and J. S. Schroeder
Citation: Journal of Physical and Chemical Reference Data 43, 043102 (2014); doi: 10.1063/1.4895394
Valid for up to 2800 bar and fromn 160 to 200 K


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
function ethanolSinglephase(T,ro)

sub=newSingle()
crit=ethanolCritical()
delta=ro/crit.ro
tau=crit.T/T
R=8.3144598
mm=0.04606844

nOne=6
nTwo=10
nThree=9
nFour=0


c=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.000, 1.000, 2.000, 1.000,
2.000, 1.000, 2.000, 1.000, 1.000, 
1.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000]

d=[4.0, 1.0, 1.0, 2.0, 2.0,
3.0, 1.0, 1.0, 1.0, 3.0,
3.0, 2.0, 2.0, 6.0, 6.0,
8.0, 1.0, 1.0, 2.0, 3.0,
3.0, 2.0, 2.0, 2.0, 1.0]

t=[1.0, 1.04, 2.72, 1.174,
1.329, 0.195, 2.43, 1.274, 
4.16, 3.3, 4.177, 2.5, 0.81,
2.02, 1.606, 0.86, 2.5, 3.72,
1.19, 3.25, 3.0, 2.0, 2.0, 1.0, 1.0]

n=[0.058200796, 0.94391227, -0.80941908, 
0.55359038, -1.4269032, 0.13448717,
0.42671978, -1.1700261, -0.92405872,
0.34891808, -0.91327720, 0.022629481,
-0.15513423, 0.21055146, -0.21997690,
-0.0065857238, 0.75564749, 0.1069411,
-0.069533844, -0.24947395, 0.027177891,
-0.00090539530, -0.12310953, -0.089779710,
-0.39512601]

al=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 
0.000, 1.075, 0.463, 0.876, 1.108, 
0.741, 4.032, 2.453, 2.300, 3.143]

bt=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.207, 0.090, 0.581, 0.947,
2.356, 27.010, 4.542, 1.287, 3.090]

gm=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.194, 1.986, 1.583, 0.756,
0.495, 1.002, 1.077, 1.493, 1.542]

ep=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.779, 0.805, 1.869, 0.694,
1.312, 2.054, 0.441, 0.793, 0.313]

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
	i+=1
	phidelta=phidelta+n[i]*exp(-delta^c[i])*(delta^(d[i]-1.0)*tau^t[i]*(d[i]-c[i]*delta^c[i]))
	phir=phir+n[i]*delta^(d[i])*tau^t[i]*exp(-delta^c[i])
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1.0)*exp(-delta^c[i]);
	phidelta2=phidelta2+n[i]*exp(-delta^c[i])*(delta^(d[i]-2)*tau^(t[i])*((d[i]-c[i]*delta^c[i])*(d[i]-1.0-c[i]*delta^c[i])-c[i]^2*delta^c[i]))
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2)*exp(-delta^c[i])
	dphideltadtau=dphideltadtau+n[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)*(d[i]-c[i]*delta^c[i])*exp(-delta^c[i])
end
for pn=1:nThree
	i+=1
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
n0=[-12.7531, 9.39094, 3.43069, 2.14326, 5.09206, 6.60138, 5.70777]
gamma0=[0.0, 0.0, 0.0, 0.816771, 2.59175, 3.80408, 8.58736]

#entropy
phi0=log(delta)+n0[1]+n0[2]*tau+n0[3]*log(tau)+n0[4]*log(1.0-exp(-gamma0[4]*tau))+n0[5]*log(1.0-exp(-gamma0[5]*tau))+n0[6]*log(1.0-exp(-gamma0[6]*tau))+n0[7]*log(1.0-exp(-gamma0[7]*tau))
phi0tau=n0[2]+n0[3]/tau+n0[4]*gamma0[4]*(1.0/(1.0-exp(-gamma0[4]*tau))-1.0)+n0[5]*gamma0[5]*(1.0/(1.0-exp(-gamma0[5]*tau))-1.0)+n0[6]*gamma0[6]*(1.0/(1.0-exp(-gamma0[6]*tau))-1.0)+n0[7]*gamma0[7]*(1.0/(1.0-exp(-gamma0[7]*tau))-1.0)
phi0delta=1/delta
phi0delta2=-1/delta^2
phi0tau2=-n0[3]/tau^2-n0[4]*gamma0[4]^2*exp(-gamma0[4]*tau)*(1.0-exp(-gamma0[4]*tau))^(-2)-n0[5]*gamma0[5]^2*exp(-gamma0[5]*tau)*(1.0-exp(-gamma0[5]*tau))^(-2)-n0[6]*gamma0[6]^2*exp(-gamma0[6]*tau)*(1.0-exp(-gamma0[6]*tau))^(-2)-n0[7]*gamma0[7]^2*exp(-gamma0[7]*tau)*(1.0-exp(-gamma0[7]*tau))^(-2)
phi0taudelta=0.0

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
calculates maximum of eos works
"""
function ethanolMax(T)
	sub = newViscMax()
	if(T < ethanolCritical().T)
		sub.p = 2800
		init = ethanolSaturate(T).lRo*1.2
	else
		sub.p = 2800
		init = ethanolCritical().ro * 1.2
	end
	sub.ro = find_zero(x -> ethanolSinglephase(T,x).p-sub.p, init)
	return sub
end


"""
uPa*s
Generalized SAFT-DFT/DMT Model for the Thermodynamic,
Interfacial, and Transport Properties of Associating Fluids:
Application for n-Alkanols
S. B. Kiselev and J. F. Ely
I. M. Abdulagatov and M. L. Huber
10.1021/ie050010e

"""
function ethanolVisc(T,ro)
	nu0=-1.03116 + 3.48379e-2*T -6.50264e-6*T^2
	b=[-19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158]
	t=[0.0, -0.25, -0.5, -0.75, -1.0, -1.25, -1.5, -2.5, -5.5]
	bb=0.0
	Tb=T/362.6
	for i=1:9
		bb+=b[i]*Tb^t[i]
	end
	br=bb*0.602*0.453^3
	
	a=[[0.0, 0.0, 0.0], [0.131194057, -0.382240694, 0.0], [-0.0805700894, 0.153811778, -0.110578307]]
	c=[23.7222995, -3.38264465, 12.7568864]
	Tr=513.9
	ror=5.991
	del=ro/ror
	tau=T/Tr
	del0=c[2]+c[3]*sqrt(tau)
	
#	println(" i $(3) j $(2) aij $(a[3][3])")
	dh=0
	for j=2:3
		for k=1:3
			dh+=a[j][k]*del^j/tau^(k-1)
		end
	end
	dh=1000*(dh+c[1]*del*(1/(del0-del)-1/del0))
#	println("nu0 $(nu0) br $(br) dh $(dh)")
	return nu0*(1.0+br*ro)+dh
end

"""
ethanol viscosity maximum
"""
function ethanolViscMax(T)
	sub = newViscMax()
	if(T < 480)
		sub.p = 1000
		init = ethanolSaturate(T).lRo*1.2
	else
		sub.p = 1000
		if(T < ethanolCritical().T)
			init = ethanolSaturate(T).lRo*1.2
		else
			init = ethanolCritical().ro * 1.2
		end
	end
	sub.ro = find_zero(x -> ethanolSinglephase(T,x).p-sub.p, init)
	return sub
end

"""
ethanol thermal conductivity, [m W m^-1 K^-1]
Reference Correlation of the Thermal Conductivity of Ethanol from the Triple Point to 600 K and up to 245 MPa
M. J. Assael, E. A. Sykioti, M. L. Huber and R. A. Perkins
10.1063/1.4797368
(Empirical critical enhancement)

"""
function ethanolTherm(T,ro)
	crit=ethanolCritical()
	Tr=T/crit.T
	ror=ro/crit.ro
	
	la0=(-2.09575 + 19.9045*Tr - 53.964*Tr^2 + 82.1223*Tr^3 - 1.98864*Tr^4 - 0.495513*Tr^5)/(0.17223 - 0.078273*Tr + Tr^2)
	B1=[2.67222e-2, 1.48279e-1, -1.30429e-1, 3.46232e-2, -2.44293e-3]
	B2=[1.77166e-2, -8.93088e-2, 6.84664e-2, -1.45702e-2, 8.09189e-4]
#	println("la0 $(la0)")
	dla=0
	for i=1:5
		dla+=(B1[i]+B2[i]*Tr)*ror^i
#		println("i $(B1[i]) dla $(B2[i])")
	end
	dlac=1.7e-3/(7.0e-2+Tr-1.0)*exp(-(1.7*(ror-1.0))^2)
#	println("dlac $(dlac)")
#more accuracy
#Rd=1.02
#v=0.63
#gam=1.239
#Gam=0.05885
#dz0=1.64296*10.0^-10
#qd=5.3*10.0^-10
#Tref=772.06
#kB=1.38064852*10.0^-23

#eos_crit=ethanolCritical()
#eos_cur=ethanolSinglephase(T,ro)

#println("$(eos_cur.dpdt) --- $(ethanolSinglephase(Tref,ro).dpdt)")
#println("$(eos_cur.dpdt-Tref/T*ethanolSinglephase(Tref,ro).dpdt)")
#println("$((v/gam))")

#if(eos_cur.dpdt-Tref/T*ethanolSinglephase(Tref,ro).dpdt>0)
#	dz=dz0*(eos_crit.p*ro/Gam/eos_crit.ro^2)^(v/gam)*(eos_cur.dpdt-Tref/T*ethanolSinglephase(Tref,ro).dpdt)^(v/gam)
#else
#	dz=0.0
#end
#if(qd*dz>1.7e-7)
#	Om0=2.0/pi*(1-exp(-1/(qd*dz)^-1+(qd*dz*eos_crit.ro/ro)^2/3))
#	Om=2/pi*((eos_cur.cp-eos_cur.cv)/eos_crit.cp*atan(qd*dz)+eos_cur.cv/eos_cur.cp*qd*dz)
#else
#	Om0=0.0
#	Om=0.0
#end


#dlac2=ro*eos_cur.cp*Rd*kB*T/6/pi/ethanolVisc(T,ro)/dz*(Om-Om0)
#println("dlac $(dlac) dlac2 $(dlac2)")

#
	return la0+dla*1000.0+dlac
end

"""


"""
function ethanolThermMax(T)
	sub = newViscMax()
	if(T < 480)
		sub.p = 2500
		init = ethanolSaturate(T).lRo*1.2
	else
		sub.p = 2500
		if(T < ethanolCritical().T)
			init = ethanolSaturate(T).lRo*1.2
		else
			init = ethanolCritical().ro * 1.2
		end
	end
	sub.ro = find_zero(x -> ethanolSinglephase(T,x).p-sub.p, init)
	return sub
	
end

"""
ethanol meltion line

"""
function ethanolIceBound(T)
	sub = newCritical()	#one point
	sub.p = 4369 * ((T/158.37)^2.6432 -1)
	init = ethanolTriple().lRo*1.02
	sub.ro = find_zero(x -> ethanolSinglephase(T,x).p-sub.p, init)


end

