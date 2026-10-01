(* Build and export our 2403 kernels.  Usage (container):
     scripts/mma/wolfram.sh mma/gtb2403/build_kernels.wls rw:/playpen1/shiqiu/gtb2403_runs/kernels:/out -- M OUTNAME [PREC]
   writes /out/OUTNAME/{manifest.json, *.txt.gz} at 40 significant digits. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "GTBKernels.wl"}]];
args = Rest[$ScriptCommandLine];
m = ToExpression[args[[1]]]; name = args[[2]]; prec = If[Length[args] >= 3, ToExpression[args[[3]]], 120];
{t1, all} = AbsoluteTiming[GTBBuildAll[m, -20, prec]];
Print["build M=", m, " prec=", prec, ": ", t1, " s"];
{t2, _} = AbsoluteTiming[GTBExport[all, FileNameJoin[{"/out", name}], 40]];
Print["export: ", t2, " s"];
Exit[0];
