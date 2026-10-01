(* Anchor tests for GTBKernels.wl (plan Task A3).  Run:
     scripts/mma/wolfram.sh mma/gtb2403/test_kernels.wls rw:/playpen1/shiqiu/gtb2403_runs/scratch_mma:/out
   Exit status = number of failed checks.  Hand formulas below use the built-in LegendreQ (type 3), while the package
   uses the hypergeometric representation, so the single-entry checks compare two independent routes. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "GTBKernels.wl"}]];
fails = 0;
check[name_, ok_] := (Print[If[TrueQ[ok], "PASS ", "FAIL "], name]; If[!TrueQ[ok], fails++]);
close[a_, b_, tol_] := Abs[a - b] <= tol Max[1, Abs[b]];

g = GTBGrid[50, -20, 120];
x = g["x"]; M = 50;
check["n0 = 39", g["n0"] == 39];
check["x1", close[x[[1]], 4.0059227368677782684`25, 10^-18]];
check["x39", close[x[[39]], 172.08533939065126099`25, 10^-18]];
check["x50", close[x[[50]], 97256.336691443898924`25, 10^-18]];
k0 = GTBK0[g];
check["K0 range", close[Min[k0], -0.0204103, 10^-6] && close[Max[k0], -0.00026181, 10^-5]];
svz = {3.11259180672*^-7, 2.07506120448*^-7, 1.55629590336*^-7, 5.61716677265*^-5, 2.80858338632*^-5,
       1.87238892422*^-5, 5.13358746336*^-5, 2.56679373168*^-5, 1.71119582112*^-5};
check["SVZSR 11 digits", And @@ MapThread[Abs[#1 - #2] <= 10^-11 Abs[#2] &, {GTBSVZ[], svz}]];
check["Lambda19(x1)^2", close[GTBLambda[19, x[[1]]]^2, 6.21643*^-66, 10^-5]];
check["Q5 hypergeometric", close[GTBQ[5, 1 + 2 x[[11]]/(x[[4]] - 4)], 9.952517107684723245837481*^-13, 10^-22]];
check["Q parity z<-1 vs LegendreQ type 3",
  And @@ Table[close[GTBQ[l, -1 - 3/7 l - 2], Re[LegendreQ[l, 0, 3, N[-1 - 3/7 l - 2, 60]]], 10^-40], {l, 0, 19}]];

(* discrete conjugate function: exact for sin(k phi), 1 <= k <= M; fails at M+1 *)
kk = GTBConjugateKernel[g]; phi = g["phi"];
conj[k_] := Max[Abs[kk . Sin[k phi] - Cos[k phi]]];
check["conjugate identity k=1,2,10,M", And @@ (conj[#] < 10^-90 & /@ {1, 2, 10, M})];
check["conjugate identity fails at k=M+1", conj[M + 1] > 10^-3];

(* B-splines: partition of unity on [0, x39] *)
bs = GTBBSplineBasis[g];
check["B-spline partition of unity",
  And @@ Table[Abs[Total[bs[s]] - 1] < 10^-80, {s, {0, 1/3, 4, 50, 171, x[[39]]}}]];

(* constant amplitude: f = 1/2 delta_l0 c_I T0, c = (5, 2, 0) for I = (0, 2, 1) *)
cI = <|0 -> 5, 2 -> 2, 1 -> 0|>;
check["T0 rows, physical", And @@ Flatten@Table[
   With[{r = GTBRow[g, i, l, {"node", 7}]},
     close[r[[1, 1]], If[l == 0, cI[i]/2 g["omega"][[7]], 0], 10^-80] && r[[2, 1]] == 0],
   {i, {0, 2, 1}}, {l, If[i == 1, {1, 3}, {0, 2}]}]];
check["T0 rows, unphysical", And @@ Flatten@Table[
   close[GTBRow[g, i, l, {"point", 3/2}][[1, 1]], If[l == 0, cI[i]/2, 0], 10^-80],
   {i, {0, 2, 1}}, {l, If[i == 1, {1, 3}, {0, 2}]}]];

(* single-entry closed forms, written from the derivation with the built-in LegendreQ *)
d = g["Delta"]; xp = g["xp"]; nu0 = -20;
aHand[l_, s_, m_] := d xp[[m]] Re[LegendreQ[l, 0, 3, 1 + 2 x[[m]]/(s - 4)]]/(Pi (s - 4));
wHand[s_, a_, b_] := d/Pi (xp[[b]]/(x[[a]] + x[[b]] + s - 4) - xp[[b]]/(x[[b]] - nu0));
(* I=0, l=2, node j=7, rho2 variable (a,b)=(3,9): k = 3 fTU + fSU + fST^T, coefficient k(3,9)+k(9,3), times omega_j *)
With[{j = 7, l = 2, a = 3, b = 9},
  Module[{s = x[[j]], ftu, fsu, fstT, k, want, got},
    ftu[m_, n_] := aHand[l, s, m] wHand[s, m, n] + aHand[l, s, n] wHand[s, n, m];
    fsu[m_, n_] := kk[[j, m]] aHand[l, s, n];
    fstT[m_, n_] := kk[[j, n]] aHand[l, s, m];
    k[m_, n_] := 3 ftu[m, n] + fsu[m, n] + fstT[m, n];
    want = g["omega"][[j]] (k[a, b] + k[b, a]);
    got = GTBRow[g, 0, l, {"node", j}][[1, GTBIndexRho2[g, a, b]]];
    check["Re_rho entry I=0 l=2 j=7 rho2(3,9)", close[got, want, 10^-40]]]];
(* I=1, l=1, s*=1/2, rho1 variable (m,n)=(4,6): k(m,n) = 2 Ku_n A_m + fTU(m,n) - fTU(n,m), fTU(m,n)=A_m W(m,n) - A_n W(n,m) *)
With[{l = 1, s = 1/2, m = 4, n = 6},
  Module[{ku, ftu, want, got},
    ku[q_] := d/Pi (xp[[q]]/(x[[q]] - s) - xp[[q]]/(x[[q]] - nu0));
    ftu[p_, q_] := aHand[l, s, p] wHand[s, p, q] - aHand[l, s, q] wHand[s, q, p];
    want = 2 ku[n] aHand[l, s, m] + ftu[m, n] - ftu[n, m];
    got = GTBRow[g, 1, l, {"point", s}][[1, GTBIndexRho1[g, m, n]]];
    check["fchi_rho entry I=1 l=1 s*=1/2 rho1(4,6)", close[got, want, 10^-40]]]];

(* export round trip at M = 10 *)
all10 = GTBBuildAll[10, -20, 80];
GTBExport[all10, "/out/test_M10", 40];
man = Import["/out/test_M10/manifest.json", "RawJSON"];
back = ToExpression /@ StringReplace[Import["/out/test_M10/Re_rho.txt.gz", {"GZIP", "Lines"}], "e" -> "*^"];
check["export shapes", man["matrices"]["Re_rho"]["shape"] == {300, 100 + 55}];
check["export round trip", Max[Abs[back - Flatten[all10["Re_rho"]]]] <= 10^-38 Max[1, Max[Abs[back]]]];

(* App. B variations: L = 24 (12 waves per isospin) and N_B = 46 *)
allv = GTBBuildAll[10, -20, 80, 24, 46];
check["L=24 rows", Dimensions[allv["Re_rho"]] == {3*12*10, 100 + 55}];
check["N_B=46 columns", Dimensions[allv["BSpline3"]] == {3*allv["grid"]["n0"], 3*46}];
check["N_B=46 partition of unity", Max[Abs[Total /@ allv["BSpline3"][[1 ;; allv["grid"]["n0"], 1 ;; 46]] - 1]] < 10^-50];

Print["failed: ", fails];
Exit[fails];
