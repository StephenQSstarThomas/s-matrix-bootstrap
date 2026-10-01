(* ::Package:: *)
(* GTBKernels: our own kernels for the arXiv:2403.10772 finite model (derivation: docs/gtb2403/DERIVATION_ZH.md).
   Independent of the authors' notebook: Q_l from the hypergeometric representation (never LegendreQ), every
   partial-wave entry from its closed form (no M^2 x M^2 products), B-splines by Cox-de Boor, SVZ targets from exact
   rational inputs.  Units m_pi = 1.  Variable order: [T0, sigma1(1..M), sigma2(1..M), rho1 row-major (first index =
   s-channel node, slow), rho2 upper triangle a<=b row-major with coefficient k_ab + k_ba off the diagonal].
   Row order of physical/unphysical matrices: [I=0: l=0,2,..,L-2 | I=2: same | I=1: l=1,3,..,L-1]. *)
BeginPackage["GTBKernels`"];

GTBGrid::usage = "GTBGrid[M, nu0, prec, L:20] -> association: M, L, nu0, prec, phi, x, xp, Delta, omega, s0, n0, K.";
GTBConjugateKernel::usage = "GTBConjugateKernel[g]: M x M discrete conjugate function, Re g_i = Sum_j K_ij Im g_j.";
GTBK0::usage = "GTBK0[g]: K0_m = (Delta/Pi) x'_m nu0/(x_m (x_m - nu0)) (F(0) = 1 with the disk centred at nu0).";
GTBQ::usage = "GTBQ[l, z]: Legendre Q_l(z) for real |z| > 1 from the hypergeometric representation.";
GTBLambda::usage = "GTBLambda[l, s] = ((Sqrt[s]-2)/(Sqrt[s]+2))^(l/2).";
GTBRow::usage = "GTBRow[g, I, l, {\"node\", j} | {\"point\", s}] -> {re, im} coefficient lists over the amplitude \
variables: h = omega f at a physical node, f at an unphysical point (im = 0).";
GTBIndexRho1::usage = "GTBIndexRho1[g, m, n]: position of rho1(x_m, y_n) in the amplitude vector.";
GTBIndexRho2::usage = "GTBIndexRho2[g, a, b]: position of the symmetric rho2(x_a, y_b) in the amplitude vector.";
GTBFormFactorKernel::usage = "GTBFormFactorKernel[g]: KFF3 = I3 (x) (K + 1 K0^T).";
GTBBSplineBasis::usage = "GTBBSplineBasis[g, nb:40]: function s -> {B_1(s), .., B_nb(s)} (cubic, clamped on [0, x_n0], nb-3 equal intervals).";
GTBMoments::usage = "GTBMoments[g] -> {intS0, intP1, intD0}, each 3 x n0, normalised by x_n0^(n + n_l).";
GTBSVZ::usage = "GTBSVZ[] -> the nine O(alpha_s) identity-operator FESR targets {S0(0,1,2), P1(-1,0,1), D0(-2,-1,0)}.";
GTBBuildAll::usage = "GTBBuildAll[M, nu0, prec, L:20, nb:40] -> association of every exported matrix.";
GTBExport::usage = "GTBExport[all, dir, digits]: manifest.json + one gzip text file per matrix.";

Begin["`Private`"];

mpi = 13957/100;                     (* MeV, exact *)
alphaS = 314588744/10^9;             (* alpha_s(2 GeV) of the authors' notebook, taken as an exact decimal *)
mu = (36/10)/mpi; md = (65/10)/mpi;  (* quark masses at 2 GeV in units of m_pi *)

GTBGrid[M_Integer, nu0_, prec_Integer, L_Integer : 20] := Module[{ph, x, s0, g},
  ph = Table[(j - 1/2) Pi/M, {j, M}];
  x = N[nu0 + (8 - 2 nu0)/(1 + Cos[ph]), prec];
  s0 = (2000/mpi)^2;
  g = <|"M" -> M, "L" -> L, "nu0" -> nu0, "prec" -> prec, "phi" -> N[ph, prec], "x" -> x,
        "xp" -> N[(8 - 2 nu0) Sin[ph]/(1 + Cos[ph])^2, prec], "Delta" -> N[Pi/M, prec],
        "omega" -> N[Pi Sqrt[(x - 4)/x], prec], "s0" -> s0, "n0" -> Count[x, v_ /; v <= s0]|>;
  Append[g, "K" -> GTBConjugateKernel[g]]];

GTBConjugateKernel[g_] := Module[{M = g["M"], p = g["prec"], c},
  c[n_] := With[{r = Mod[n, 2 M]}, If[EvenQ[r], 0, N[Cot[r Pi/(2 M)]/M, p]]];
  Table[-c[i - j] + c[i + j - 1], {i, M}, {j, M}]];

GTBK0[g_] := g["Delta"]/Pi g["xp"] g["nu0"]/(g["x"] (g["x"] - g["nu0"]));

qPos[l_, z_] := Sqrt[Pi] Gamma[l + 1]/(Gamma[l + 3/2] (2 z)^(l + 1)) *
   Hypergeometric2F1[(l + 1)/2, (l + 2)/2, l + 3/2, 1/z^2];
GTBQ[l_Integer, z0_] := With[{z = If[Precision[z0] === Infinity, N[z0, 120], z0]},
   If[z > 1, qPos[l, z], (-1)^(l + 1) qPos[l, -z]]];

GTBLambda[l_, s_] := ((Sqrt[s] - 2)/(Sqrt[s] + 2))^(l/2);

GTBIndexRho1[g_, m_, n_] := 1 + 2 g["M"] + (m - 1) g["M"] + n;
GTBIndexRho2[g_, a0_, b0_] := With[{a = Min[a0, b0], b = Max[a0, b0], M = g["M"]},
   1 + 2 M + M^2 + (a - 1) M - (a - 1) (a - 2)/2 + (b - a) + 1];

(* The projections of the seven pieces of A(s,t,u) at one s (physical node or unphysical point).
   Index [m, n] of a double kernel multiplies rho(x_m, y_n). *)
slots[g_, l_, where_] := Module[{M = g["M"], x = g["x"], xp = g["xp"], d = g["Delta"], nu0 = g["nu0"],
    phys, s, dl, sg, a, w0, wm, ksRe, ksIm, fSre, fSim, fT, fU, fTU},
  phys = where[[1]] === "node";
  s = If[phys, x[[where[[2]]]], N[where[[2]], g["prec"]]];
  dl = If[l == 0, 1, 0]; sg = (-1)^l;
  a = d xp (GTBQ[l, #] & /@ (1 + 2 x/(s - 4)))/(Pi (s - 4));          (* A_l(s, x_m) with Delta x'_m *)
  w0 = d xp/(Pi (x - nu0));                                            (* W0_m *)
  wm = Table[d/Pi (xp[[b]]/(x[[aa]] + x[[b]] + s - 4) - xp[[b]]/(x[[b]] - nu0)), {aa, M}, {b, M}];
  If[phys,
    ksRe = g["K"][[where[[2]]]]; ksIm = UnitVector[M, where[[2]]],
    ksRe = d/Pi (xp/(x - s) - xp/(x - nu0)); ksIm = ConstantArray[0, M]];
  fSre = dl/2 ksRe; fSim = dl/2 ksIm;
  fT = a - dl/2 w0; fU = sg a - dl/2 w0;
  fTU = a wm + sg Transpose[a wm] + dl/2 Outer[Times, w0, w0];
  <|"dl" -> dl, "fSre" -> fSre, "fSim" -> fSim, "fT" -> fT, "fU" -> fU,
    "fSTre" -> Outer[Times, ksRe, fT], "fSTim" -> Outer[Times, ksIm, fT],
    "fSUre" -> Outer[Times, ksRe, fU], "fSUim" -> Outer[Times, ksIm, fU], "fTU" -> fTU, "phys" -> phys, "s" -> s|>];

comb[i_, s_, t_, u_] := Switch[i, 0, 3 s + t + u, 2, t + u, 1, t - u];
fold[k_, M_] := Flatten[Table[If[a == b, k[[a, a]], k[[a, b]] + k[[b, a]]], {a, M}, {b, a, M}]];

rowFromSlots[g_, i_, sl_, scale_] := Module[{M = g["M"], z, t0, s1re, s1im, s2re, s2im, r1re, r1im, r2re, r2im},
  z = ConstantArray[0, M];
  t0 = comb[i, sl["dl"]/2, sl["dl"]/2, sl["dl"]/2];
  s1re = comb[i, sl["fSre"], sl["fT"], sl["fU"]]; s1im = comb[i, sl["fSim"], z, z];
  s2re = comb[i, sl["fT"] + sl["fU"], sl["fSre"] + sl["fU"], sl["fSre"] + sl["fT"]];
  s2im = comb[i, z, sl["fSim"], sl["fSim"]];
  r1re = comb[i, sl["fSTre"] + sl["fSUre"], Transpose[sl["fSTre"]] + sl["fTU"], Transpose[sl["fTU"] + sl["fSUre"]]];
  r1im = comb[i, sl["fSTim"] + sl["fSUim"], Transpose[sl["fSTim"]], Transpose[sl["fSUim"]]];
  r2re = comb[i, sl["fTU"], sl["fSUre"], Transpose[sl["fSTre"]]];
  r2im = comb[i, 0 sl["fTU"], sl["fSUim"], Transpose[sl["fSTim"]]];
  {scale Join[{t0}, s1re, s2re, Flatten[r1re], fold[r2re, M]],
   scale Join[{0}, s1im, s2im, Flatten[r1im], fold[r2im, M]]}];

GTBRow[g_, i_, l_, where_] := With[{sl = slots[g, l, where]},
   rowFromSlots[g, i, sl, If[sl["phys"], g["omega"][[where[[2]]]], 1]]];

GTBFormFactorKernel[g_] := KroneckerProduct[IdentityMatrix[3], g["K"] + ConstantArray[GTBK0[g], g["M"]]];

bsp[t_, i_, 0, s_, end_] := If[(t[[i]] <= s < t[[i + 1]]) || (s == end && t[[i]] < t[[i + 1]] == end), 1, 0];
bsp[t_, i_, p_, s_, end_] :=
  If[t[[i + p]] == t[[i]], 0, (s - t[[i]])/(t[[i + p]] - t[[i]]) bsp[t, i, p - 1, s, end]] +
  If[t[[i + p + 1]] == t[[i + 1]], 0, (t[[i + p + 1]] - s)/(t[[i + p + 1]] - t[[i + 1]]) bsp[t, i + 1, p - 1, s, end]];
GTBBSplineBasis[g_, nb_Integer : 40] := With[{e = g["x"][[g["n0"]]]},
   With[{t = Join[{0, 0, 0, 0}, Table[k e/(nb - 3), {k, nb - 4}], {e, e, e, e}]},
     Function[s, Table[bsp[t, i, 3, s, e], {i, 1, nb}]]]];

kap[0, x_] := 3 Sqrt[x - 4]/(256 Pi^5 Sqrt[x]);
kap[1, x_] := (x - 4)^(3/2)/(384 Pi^5 Sqrt[x]);
kap[2, x_] := (x - 4)^(5/2)/(2560 Pi^5 Sqrt[x]);
GTBMoments[g_] := Module[{n0 = g["n0"], x, xp, d = g["Delta"], e, mom},
  x = Take[g["x"], n0]; xp = Take[g["xp"], n0]; e = x[[n0]];
  mom[l_, n_, nl_] := kap[l, x] d xp x^n/e^(n + nl);
  {Table[mom[0, n, 2], {n, 0, 2}], Table[mom[1, n, 2], {n, -1, 1}], Table[mom[2, n, 3], {n, -2, 0}]}];

GTBSVZ[] := N[Join[
   Table[(mu^2 + md^2)/(2 Pi)^4 3/(4 Pi (n + 2)) (1 + 13/3 alphaS/Pi), {n, 0, 2}],
   Table[1/(2 Pi)^4 1/(4 Pi (n + 2)) (1 + alphaS/Pi), {n, -1, 1}],
   Table[1/(2 Pi)^4 1/(4 Pi (n + 3)) (11/10 - 17/18 alphaS/Pi), {n, -2, 0}]], 60];

waves[g_] := Join[Table[{0, l}, {l, 0, g["L"] - 2, 2}], Table[{2, l}, {l, 0, g["L"] - 2, 2}],
   Table[{1, l}, {l, 1, g["L"] - 1, 2}]];

GTBBuildAll[M_Integer, nu0_, prec_Integer, L_Integer : 20, nb_Integer : 40] :=
 Module[{g, ws, keep, phys, chi, ff, bs, mom},
  g = Append[GTBGrid[M, nu0, prec, L], "nb" -> nb]; ws = waves[g]; keep = N[#, Min[prec, 60]] &;
  (* slots depend on (l, where) only; I=0 and I=2 share the even-l slots *)
  phys = Table[With[{sl = slots[g, w[[2]], {"node", j}]},
       keep /@ rowFromSlots[g, w[[1]], sl, g["omega"][[j]]]], {w, ws}, {j, M}];
  chi = Table[With[{sl = slots[g, w[[2]], {"point", s}]}, keep[rowFromSlots[g, w[[1]], sl, 1][[1]]]],
     {w, ws}, {s, {1/2, 1, 3/2, 2}}];
  ff = Table[With[{sl = slots[g, w[[2]], {"point", 3}]}, keep[rowFromSlots[g, w[[1]], sl, 1][[1]]]], {w, ws}];
  phys = Flatten[phys, 1]; bs = GTBBSplineBasis[g, nb]; mom = GTBMoments[g];
  <|"grid" -> g,
    "x" -> keep[g["x"]], "xp" -> keep[g["xp"]], "K" -> keep[g["K"]], "K0" -> keep[GTBK0[g]],
    "Lam" -> keep[Flatten[Table[GTBLambda[w[[2]], g["x"][[j]]], {w, ws}, {j, M}]]],
    "Re_sig" -> phys[[All, 1, 1 ;; 2 M + 1]], "Im_sig" -> phys[[All, 2, 1 ;; 2 M + 1]],
    "Re_rho" -> phys[[All, 1, 2 M + 2 ;;]], "Im_rho" -> phys[[All, 2, 2 M + 2 ;;]],
    "fchi_sig" -> Flatten[chi, 1][[All, 1 ;; 2 M + 1]], "fchi_rho" -> Flatten[chi, 1][[All, 2 M + 2 ;;]],
    "fF_sig" -> ff[[All, 1 ;; 2 M + 1]], "fF_rho" -> ff[[All, 2 M + 2 ;;]],
    "Lambda3" -> keep[Flatten[Table[GTBLambda[l, g["x"][[j]]], {l, 0, 2}, {j, g["n0"]}]]],
    "KFF3" -> keep[GTBFormFactorKernel[g]],
    "BSpline3" -> keep[KroneckerProduct[IdentityMatrix[3], Table[bs[g["x"][[j]]], {j, g["n0"]}]]],
    "intS0" -> keep[mom[[1]]], "intP1" -> keep[mom[[2]]], "intD0" -> keep[mom[[3]]], "SVZSR" -> GTBSVZ[]|>];

fmt[v_, digits_] := If[v == 0, "0", Module[{m, e, d},
   {m, e} = MantissaExponent[v]; d = Round[m 10^digits];
   If[Abs[d] >= 10^digits, d = Round[d/10]; e = e + 1];
   ToString[d] <> "e" <> ToString[e - digits]]];

GTBExport[all_, dir_String, digits_Integer] := Module[{g = all["grid"], names, man = <||>, file},
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  names = DeleteCases[Keys[all], "grid"];
  Do[file = FileNameJoin[{dir, nm <> ".txt.gz"}];
    Export[file, fmt[#, digits] & /@ Flatten[all[nm]], {"GZIP", "Lines"}];
    man[nm] = <|"shape" -> Dimensions[all[nm]], "file" -> nm <> ".txt.gz",
       "sha256" -> FileHash[file, "SHA256", All, "HexString"]|>, {nm, names}];
  Export[FileNameJoin[{dir, "manifest.json"}],
    <|"M" -> g["M"], "L" -> g["L"], "nb" -> Lookup[g, "nb", 40], "nu0" -> g["nu0"], "n0" -> g["n0"], "digits" -> digits, "precision" -> g["prec"],
      "generator" -> "ours", "package" -> "mma/gtb2403/GTBKernels.wl", "matrices" -> man|>, "RawJSON"]];

End[];
EndPackage[];
