(* Independent 2309v3 equations only. No later-author notebook is loaded. *)
args = Rest[$ScriptCommandLine];
If[Length[args] != 2, Print["usage: audit_2309.wls M OUT.json"]; Exit[2]];
m = ToExpression[args[[1]]]; out = args[[2]];
If[!IntegerQ[m] || m < 2, Exit[2]];
phi = Table[Pi (i - 1/2)/m, {i, m}];
nodes = N[8/(1 + Cos[phi]), 100];
weights = N[Pi/m 8 Sin[phi]/(1 + Cos[phi])^2, 100];
kt[0] = 0;
kt[j_Integer] := (1 - (-1)^j)/(2 m) Cot[j Pi/(2 m)];
kernel = N[Table[kt[i+j-2 m-1] - kt[i-j], {i,m}, {j,m}],100];
kin0 = Sqrt[6 Pi]/(16 Pi^3) nodes^(-1/4) ((nodes-4)/4)^(1/4);
kin1 = Sqrt[4 Pi/3]/(8 Pi^3) nodes^(-1/4) ((nodes-4)/4)^(3/4);
s0 = 3600/49;
targets = {44187/10^11 s0^2, 282014/10^12 s0^3,
           575484/10^10 s0, 269948/10^10 s0^2};
(* Exact unitary change of basis for the original complex Gram. *)
u = {{1/Sqrt[2], I/Sqrt[2], 0}, {1/Sqrt[2], -I/Sqrt[2], 0}, {0,0,1}};
gram = {{1,sr+I si,fr+I fi},{sr-I si,1,fr-I fi},{fr-I fi,fr+I fi,rho}};
realGram = {{1+sr,si,Sqrt[2] fr},{si,1-sr,Sqrt[2] fi},{Sqrt[2] fr,Sqrt[2] fi,rho}};
gramIdentity = Simplify[ConjugateTranspose[u].gram.u-realGram] === ConstantArray[0,{3,3}];
doc = <|"source"->"2309.12402v3 (2.33), (2.56), (3.60-61), (3.67-68)",
 "version"->$Version, "working_precision"->100, "gram_identity"->gramIdentity,
 "nodes"->N[nodes,17], "weights"->N[weights,17], "hilbert"->N[kernel,17],
 "kinematic"->N[{kin0,kin1},17], "targets"->N[targets,17]|>;
Export[out, doc, "RawJSON"];
Print["2309 audit M=",m," Gram identity=",gramIdentity," -> ",out];
If[!gramIdentity, Exit[1]];
