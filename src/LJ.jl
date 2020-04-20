"""
Lennard Jones equation of state

"""
function LJSinglephase(T,ro)
	x=[0.8623085097507421; 2.976218765822098; -8.402230115796038; 0.1054136629203555; -0.8564583828174598;
1.582759470107601; 0.7639421948305453; 1.753173414312048; 2.798291772190376e03; -4.8394220260857657e-02;
0.9963265197721935; -3.698000291272493e01; 2.084012299434647e01; 8.305402124717285e01; -9.574799715203068e02;
-1.477746229234994e02; 6.398607852471505e01; 1.603993673294834e01; 6.805916615864377e01; -2.791293578795945e03;
-6.245128304568454; -8.116836104958410e03; 1.488735559561229e01; -1.059346754655084e04; -1.131607632802822e02;
-8.867771540418822e03; -3.986982844450543e01; -4.689270299917261e03; 2.593535277438717e02; -2.6945235894349033e03;
-7.218487631550215e02; 1.721802063863269e02;]
	
	a=Array{Float64}(undef, 8)
	b=Array{Float64}(undef, 6)
	G=Array{Float64}(undef, 6)
	
	a[1]=x[1]*T+x[2]*sqrt(T)+x[3]+x[4]/T+x[5]/(T^2);
	a[2]=x[6]*T+x[7]+x[8]/T+x[9]/(T^2);
	a[3]=x[10]*T+x[11]+x[12]/T;
	a[4]=x[13];
	a[5]=x[14]/T+x[15]/T^2;
	a[6]=x[16]/T;
	a[7]=x[17]/T+x[18]/T^2;
	a[8]=x[19]/T^2;
	
	b[1]=x[20]/T^2+x[21]/T^3;
	b[2]=x[22]/T^2+x[23]/T^4;
	b[3]=x[24]/T^2+x[25]/T^3;
	b[4]=x[26]/T^2+x[27]/T^4;
	b[5]=x[28]/T^2+x[29]/T^3;
	b[6]=x[30]/T^2+x[31]/T^3+x[32]/T^4;

	gm=3.0
	F=exp(-gm*ro*ro)
	G[1]=(1.0-F)/(2*gm)
	G[2]=-(F*ro^2-2*G[1])/(2*gm)
	G[3]=-(F*ro^4-4*G[2])/(2*gm)
	G[4]=-(F*ro^6-6*G[3])/(2*gm)
	G[5]=-(F*ro^8-8*G[4])/(2*gm)
	G[6]=-(F*ro^10-10*G[5])/(2*gm)
	
	pres=ro*T+a[1]*ro^2+a[2]*ro^3+a[3]*ro^4+a[4]*ro^5+a[5]*ro^6+a[6]*ro^7+a[7]*ro^8+a[8]*ro^9+F*(b[1]*ro^3+b[2]*ro^5+b[3]*ro^7+b[4]*ro^9+b[5]*ro^11+b[6]*ro^13)
	Ar=a[1]*ro^1+a[2]*ro^2/2.0+a[3]*ro^3/3.0+a[4]*ro^4/4.0+a[5]*ro^5/5.0+a[6]*ro^6/6.0+a[7]*ro^7/7.0+a[8]*ro^8/8.0+b[1]*G[1]+b[2]*G[2]+b[3]*G[3]+b[4]*G[4]+b[5]*G[5]+b[6]*G[6]
	Ge=Ar+pres/ro-T
	dmudn=a[1]+a[2]*ro+a[3]*ro^2+a[4]*ro^3+a[5]*ro^4+a[6]*ro^5+a[7]*ro^6+a[8]*ro^7+F*(b[1]*ro+b[2]*ro^3+b[3]*ro^5+b[4]*ro^7+b[5]*ro^9+b[6]*ro^11)
	dmudn+=a[1]+a[2]*2.0*ro+a[3]*3.0*ro^2+a[4]*4.0*ro^3+a[5]*5.0*ro^4+a[6]*6.0*ro^5+a[7]*7.0*ro^6+a[8]*8.0*ro^7-F*2*gm*(b[1]*ro^3+b[2]*ro^5+b[3]*ro^7+b[4]*ro^9+b[5]*ro^11+b[6]*ro^13)+2*F*(b[1]*ro+b[2]*2.0*ro^3+b[3]*3.0*ro^5+b[4]*4.0*ro^7+b[5]*5.0*ro^9+b[6]*6.0*ro^11)

	return [pres; Ar; Ge; dmudn]
end


