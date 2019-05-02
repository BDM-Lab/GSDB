%SimulationDemo
noiserate=0;
n=100;
method_type=0;

%%Helix Curve
 [binAnno_Helix FreqMat_Helix XX_Helix]=helixdata(n,noiserate);
 ChromSDE(binAnno_Helix,FreqMat_Helix,method_type);
%%Brownian Curve
[binAnno_Brownian FreqMat_Brownian XX_Brownian]=randwalkdata(n,noiserate);
ChromSDE(binAnno_Brownian,FreqMat_Brownian,method_type);

%%Uniform Random Points
 [binAnno_Rand FreqMat_Rand XX]=randdata(n,noiserate);
 ChromSDE(binAnno_Rand,FreqMat_Rand,method_type);