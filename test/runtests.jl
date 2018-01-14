using PureEOS
using Base.Test

# write your own tests here
println("Viscosity test")
@test waterVisc(298.15,998/18.015268) ≈ 889.735100 atol=0.001
@test waterVisc(298.15,1200/18.015268) ≈ 1437.649467 atol=0.001
@test waterVisc(373.15,1000/18.015268) ≈ 307.883622 atol=0.001
@test waterVisc(433.15,1/18.015268) ≈ 14.538324 atol=0.001
@test waterVisc(433.15,1000/18.015268) ≈ 217.685358 atol=0.001
@test waterVisc(873.15,1/18.015268) ≈ 32.619287 atol=0.001
@test waterVisc(873.15,100/18.015268) ≈ 35.802262 atol=0.001
@test waterVisc(873.15,600/18.015268) ≈ 77.430195 atol=0.001
@test waterVisc(1173.15,1/18.015268) ≈ 44.217245 atol=0.001
@test waterVisc(1173.15,100/18.015268) ≈ 47.640433 atol=0.001
@test waterVisc(1173.15,400/18.015268) ≈ 64.154608 atol=0.001
println("thermalconductivity test")
@test waterTherm(647.35,1/18.015268) ≈ 51.9298924 atol=0.001
@test waterTherm(647.35,122/18.015268) ≈ 130.922885 atol=0.1
@test waterTherm(647.35,222/18.015268) ≈ 367.787459 atol=20
@test waterTherm(647.35,272/18.015268) ≈ 757.959776 atol=40
@test waterTherm(647.35,322/18.015268) ≈ 1443.75556 atol=200
@test waterTherm(647.35,372/18.015268) ≈ 650.319402 atol=40
@test waterTherm(647.35,422/18.015268) ≈ 448.883487 atol=20
@test waterTherm(647.35,750/18.015268) ≈ 600.961346 atol=0.001

#Data from NIST
