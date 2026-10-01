(* Build and export our 2403 kernels.  Usage (container):
     scripts/mma/wolfram.sh mma/gtb2403/build_kernels.wls rw:/playpen1/shiqiu/gtb2403_runs/kernels:/out -- M OUTNAME [PREC [L [NB [NSR]]]]
   writes /out/OUTNAME/{manifest.json, *.txt.gz} at 40 significant digits. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "GTBKernels.wl"}]];
args = Rest[$ScriptCommandLine];
m = ToExpression[args[[1]]]; name = args[[2]]; prec = If[Length[args] >= 3, ToExpression[args[[3]]], 120];
lw = If[Length[args] >= 4, ToExpression[args[[4]]], 20]; nb = If[Length[args] >= 5, ToExpression[args[[5]]], 40];
nsr = If[Length[args] >= 6, ToExpression[args[[6]]], 3];
{t1, all} = AbsoluteTiming[GTBBuildAll[m, -20, prec, lw, nb, nsr]];
Print["build M=", m, " prec=", prec, ": ", t1, " s"];
t2 = First@AbsoluteTiming[GTBExport[all, FileNameJoin[{"/out", name}], 40]];
Print["export: ", t2, " s"];
Exit[0];
