"""
Stockmayer potential
T --- temperature [K]
n --- densty [mol/l]
del --- charge  []
sigma --- [nm]
eps --- epsilon/kB [K]
Mm -- molar mass  []

"""
function StECH(T::Number, n::Number, del::Number, sigma::Number, eps::Number, Mm::Number)
omdel11=[0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5];
omtemp11=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3, 3.5, 
4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 75, 100 ]

oms11=[ 4.0079  4.002  4.655  5.521  6.454  8.214  9.824  11.31;
 3.13  3.164  3.355  3.721  4.198  5.23  6.225  7.16;
 2.6494  2.657  2.77  3.002  3.319  4.054  4.785  5.483;
 2.3144  2.32  2.402  2.572  2.812  3.386  3.972  4.539;
 2.0661  2.073  2.14  2.278  2.472  2.946  3.437  3.918;
 1.8767  1.885  1.944  2.06  2.225  2.628  3.054  3.474;
 1.7293  1.738  1.791  1.893  2.036  2.388  2.763  3.137;
 1.6122  1.622  1.67  1.76  1.886  2.198  2.535  2.872;
 1.5175  1.527  1.572  1.653  1.765  2.044  2.349  2.657;
 1.4398  1.45  1.49  1.564  1.665  1.917  2.196  2.478;
1.3204  1.33  1.364  1.425  1.509  1.72  1.956  2.199;
 1.2336  1.242  1.272  1.324  1.394  1.573  1.777  1.99;
 1.1679  1.176  1.202  1.246  1.306  1.461  1.639  1.827;
1.1166  1.124  1.146  1.185  1.237  1.372  1.53  1.698;
1.0753  1.082  1.102  1.135  1.181  1.3  1.441  1.592;
 1.0006  1.005  1.02  1.044  1.08  1.17  1.278  1.397;
0.95003  0.9538  0.9656  0.9852  1.012  1.082  1.168  1.265;
0.91311  0.9162  0.9256  0.9413  0.9626  1.019  1.09  1.17;
 0.88453  0.8871  0.8948  0.9076  0.9252  0.9721  1.031  1.098;
0.84277  0.8446  0.8501  0.8592  0.8716  0.9053  0.9483  0.9984;
0.81287  0.8142  0.8183  0.8251  0.8344  0.8598  0.8927  0.9316;
0.78976  0.7908  0.794  0.7993  0.8066  0.8265  0.8526  0.8836;
0.77111  0.772  0.7745  0.7788  0.7846  0.8007  0.8219  0.8474;
0.75553  0.7562  0.7584  0.7619  0.7667  0.78  0.7976  0.8189;
0.7422  0.7428  0.7446  0.7475  0.7515  0.7627  0.7776  0.7957;
0.72022  0.7206  0.722  0.7241  0.7271  0.7354  0.7464  0.76;
0.70254  0.7029  0.7039  0.7055  0.7078  0.7142  0.7228  0.7334;
0.68776  0.688  0.6888  0.6901  0.6919  0.697  0.704  0.7125;
0.6751  0.6753  0.676  0.677  0.6785  0.6827  0.6884  0.6955;
0.66405  0.6642  0.6648  0.6657  0.6669  0.6704  0.6752  0.6811;
0.64136  0.6415  0.6418  0.6425  0.6433  0.6457  0.649  0.6531;
0.6235  0.6236  0.6239  0.6243  0.6249  0.6267  0.6291  0.6321;
0.60882  0.6089  0.6091  0.6094  0.6099  0.6112  0.6131  0.6154;
0.5964  0.5964  0.5966  0.5969  0.5972  0.5983  0.5998  0.6017;
0.57626  0.5763  0.5764  0.5766  0.5768  0.5775  0.5785  0.5798;
0.54146  0.5415  0.5416  0.5416  0.5418  0.5421  0.5424  0.5429;
0.51803  0.5181  0.5182  0.5184  0.5184  0.5185  0.5186  0.5187]

omdel22=[0, 0.25, 0.5, 0.75, 1.0, 1.5, 2, 2.5]
omtemp22=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4, 
5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 75, 100 ]

oms22=[4.1005  4.266  4.833  5.742  6.729  8.624  10.34  11.89;
3.2626  3.305  3.516  3.914  4.433  5.57  6.637  7.618;
2.8399  2.836  2.936  3.168  3.511  4.329  5.126  5.874;
2.531  2.522  2.586  2.749  3.004  3.64  4.282  4.895;
2.2837  2.277  2.329  2.46  2.665  3.187  3.727  4.249;
2.0838  2.081  2.13  2.243  2.417  2.862  3.329  3.786;
1.922  1.924  1.97  2.072  2.225  2.614  3.028  3.435;
1.7902  1.795  1.84  1.934  2.07  2.417  2.788  3.156;
1.6823  1.689  1.733  1.82  1.944  2.258  2.596  2.933;
1.5929  1.601  1.644  1.725  1.838  2.124  2.435  2.746;
1.4551  1.465  1.504  1.574  1.67  1.913  2.181  2.451;
1.3551  1.365  1.4  1.461  1.544  1.754  1.989  2.228;
1.28  1.289  1.321  1.374  1.447  1.63  1.838  2.053;
1.2219  1.231  1.259  1.306  1.37  1.532  1.718  1.912;
1.1757  1.184  1.209  1.251  1.307  1.451  1.618  1.795;
1.0933  1.1  1.119  1.15  1.193  1.304  1.435  1.578;
1.0388  1.044  1.059  1.083  1.117  1.204  1.31  1.428;
0.99963  1.004  1.016  1.035  1.062  1.133  1.22  1.319;
0.96988  0.9732  0.983  0.9991  1.021  1.079  1.153  1.236;
0.92676  0.9291  0.936  0.9473  0.9628  1.005  1.058  1.121;
0.89616  0.8979  0.903  0.9114  0.923  0.9545  0.9955  1.044;
0.87272  0.8741  0.878  0.8845  0.8935  0.9181  0.9505  0.9893;
0.85379  0.8549  0.858  0.8632  0.8703  0.8901  0.9164  0.9482;
0.83795  0.8388  0.8414  0.8456  0.8515  0.8678  0.8895  0.916;
0.82435  0.8251  0.8273  0.8308  0.8356  0.8493  0.8676  0.8901;
0.80184  0.8024  0.8039  0.8065  0.8101  0.8201  0.8337  0.8504;
0.78363  0.784  0.7852  0.7872  0.7899  0.7976  0.8081  0.8212;
0.76834  0.7687  0.7696  0.7712  0.7733  0.7794  0.7878  0.7983;
0.75518  0.7554  0.7562  0.7575  0.7592  0.7642  0.7711  0.7797;
0.74364  0.7438  0.7445  0.7455  0.747  0.7512  0.7569  0.7642;
0.71982  0.72  0.7204  0.1211  0.7221  0.725  0.7289  0.7339;
0.70097  0.7011  0.7014  0.7019  0.7026  0.7047  0.7076  0.7112;
0.68545  0.6855  0.6858  0.6861  0.6867  0.6883  0.6905  0.6932;
0.67232  0.6724  0.6726  0.6728  0.6733  0.6745  0.6762  0.6784;
0.65099  0.651  0.6512  0.6513  0.6516  0.6524  0.6534  0.6546;
0.61397  0.6141  0.6143  0.6145  0.6147  0.6148  0.6148  0.6147;
0.5887  0.5889  0.5894  0.59  0.5903  0.5901  0.5895  0.5885]

spline11 = Spline2D(omtemp11, omdel11, oms11)
spline22 = Spline2D(omtemp22, omdel22, oms22)

sigw::Float64=sigma*10^(-9)	#to [m]


#condstants
Na::Float64 = 6.02 * 10.0^23
kb::Float64 = 1.38*10^-23;

#dimensionless units
tb::Float64 = T/eps
nb::Float64 = n*Na*1000.0*sigw^3

#integrals
omega11::Float64 = evaluate(spline11, tb, nb)
omega22::Float64 = evaluate(spline22, tb, nb)

#println("$(omega11) Na $(Na) $(pi) sigma $(sigma) tb $(tb) nb $(nb)")

D0::Float64=3.0/(8.0*sqrt(pi)*nb)*sqrt(tb)/omega11;
#println("D0 $(D0)")
D::Float64=D0*sigw/sqrt(Mm/1000/Na/(eps*kb))*10^9;
#println("D $(D)")


return [D, 0.1]
end

