(* Independent direct angular integrals for the declared analytic interpolation. *)
args=Rest[$ScriptCommandLine]; If[Length[args]!=1,Exit[2]];
m=8; wp=40; phi=Pi/m (Range[m]-1/2);
d=N[Table[If[n==m,1,2]/m Sin[n phi[[j]]],{n,m},{j,m}],50];
g[j_,z_?NumericQ]:=Sum[d[[n,j]](z^n-(-1)^n),{n,m}];
z[v_?NumericQ]:=(2-Sqrt[4-v])/(2+Sqrt[4-v]);
s=N[8/(1+Cos[phi[[5]]]),50]; zs=N[Exp[I phi[[5]]],50];
t[u_?NumericQ]:=-(s-4)/2 (1-u); q[u_?NumericQ]:=-(s-4)/2 (1+u);
sig=3/2 g[5,zs]+NIntegrate[g[5,z[t[u]]]+g[5,z[q[u]]],{u,-1,1},
 WorkingPrecision->wp,AccuracyGoal->25,PrecisionGoal->25,
 Method->{"GlobalAdaptive","SymbolicProcessing"->0}]/4;
rho=NIntegrate[u (g[2,z[t[u]]](g[5,zs]+g[5,z[q[u]]])-g[2,z[q[u]]](g[5,z[t[u]]]+g[5,zs])),{u,-1,1},
 WorkingPrecision->wp,AccuracyGoal->25,PrecisionGoal->25,
 Method->{"GlobalAdaptive","SymbolicProcessing"->0}]/4;
num[x_]:=Module[{r=RealDigits[Abs[N[x,30]],10,30]},
 If[x<0,"-",""]<>"0."<>StringJoin[ToString/@r[[1]]]<>"e"<>ToString[r[[2]]]];
Export[args[[1]],<|"M"->m,"node_zero_based"->4,"working_precision"->wp,
 "sigma1_node4_S0"->{num[Re[sig]],num[Im[sig]]},
 "rho1_1_4_P1"->{num[Re[rho]],num[Im[rho]]},"version"->$Version|>,"RawJSON"];
Print["Analytic-family direct angular controls written to ",args[[1]]];
