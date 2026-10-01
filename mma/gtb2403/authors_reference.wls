#!/usr/bin/env wolframscript
(* authors_reference.wls -- Task A4 (Phase A plan, 2026-10-01).

   Evaluates the authors' arXiv:2403.10772 notebook GTB_numerics.nb UNCHANGED and
   headless, then writes its matrices in the common kernel-export format, as the
   reference leg of the three-way kernel comparison. The notebook stays in the
   pinned submodule (third_party/gauge-theory-bootstrap @ 801684d); nothing of it
   is copied here.

   usage:
     scripts/mma/wolfram.sh mma/gtb2403/authors_reference.wls \
       ro:third_party/gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772:/nb \
       rw:/playpen1/shiqiu/gtb2403_runs/kernels:/out \
       -- nb=/nb/GTB_numerics.nb out=/out/authors_M50 [M=10] [L=20] [force]
     (key=value arguments: wolframscript itself consumes options such as --nb)

   Evaluation contract
   - Every Input/Code cell is evaluated in document order with
     ToExpression[boxes, StandardForm] in context Global`.
   - Skipped: the SetDirectory[NotebookDirectory[]] cell and the final Export[...]
     cell (the .mat output). No formula is edited; precision stays p = 100.
   - The ONLY modification, used for cheap shape/timing runs: with M=<value> (L=<value>) the
     statement `M = <value>` (`L = <value>`) is evaluated right after the notebook
     cell that sets M = 50; L = 20. An M = 50 run passes no override. The manifest
     records which case applies.
   - $HistoryLength = 0 (a kernel setting, not part of the notebook).

   Export format (Task A3 "Interfaces"): manifest.json + one <name>.txt.gz per
   matrix, row-major, one number per line, 40 significant digits written as
   <signed integer mantissa>e<exponent>, "0" for any x with x == 0. Every number
   is rounded from its exact binary value (SetPrecision[x, Infinity]); machine
   numbers (e.g. SVZSR) therefore also get 40 digits, and the manifest lists them.
   LegendreQ[l, 0, 3, z] returns Complex numbers whose imaginary part is an
   approximate zero; such entries are written as their real part, and the
   manifest counts them (and any imaginary part that is NOT == 0).
   Exported names <- authors' symbols: x <- vx, xp <- vdx\[Phi], K, K0,
   Lam <- \[CapitalLambda]I, Re_sig <- a, Im_sig <- b, Re_rho <- A, Im_rho <- B
   (unrescaled, i.e. before /\[CapitalLambda]I), fchi_sig <- fsigmachi,
   fchi_rho <- frhochi, fF_sig <- fsigmaF, fF_rho <- frhoF,
   Lambda3 <- \[CapitalLambda]3, KFF3, BSpline3 <- BSplineM3, intS0 <- vintS0,
   intP1 <- vintP1, intD0 <- vintD0, SVZSR. *)

Begin["AuthorsRef`"];
$HistoryLength = 0;

arArgs = Rest[$ScriptCommandLine];
arOpt[name_String, default_] := With[{hits = Select[arArgs, StringStartsQ[#, name <> "="] &]},
   If[hits === {}, default, StringDrop[Last[hits], StringLength[name] + 1]]];
arIntOpt[name_String] := With[{s = arOpt[name, None]},
   Which[s === None, None, StringMatchQ[s, DigitCharacter ..], FromDigits[s],
    True, Print["bad ", name, "=", s]; Exit[2]]];

arLogFile = None;
arLog[s_String] := Module[{line, st},
   line = DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}] <> "  " <> s;
   Print[line];
   If[arLogFile =!= None,
    st = OpenAppend[arLogFile, CharacterEncoding -> "UTF-8"]; WriteString[st, line, "\n"]; Close[st]]];
arFail[s_String] := (arLog["FATAL: " <> s]; Exit[1]);

(* value of a notebook symbol, resolved exactly as the notebook resolved it
   (e.g. K is System`K, a is Global`a) *)
arNb[name_String] := Block[{$Context = "Global`"}, ToExpression[name]];

(* ---------- cell collection (numbering identical to aux/nb2wl.wls) ---------- *)
arTextStyles = {"Title", "Subtitle", "Chapter", "Section", "Subsection", "Subsubsection",
   "Subsubsubsection", "Text", "Item", "ItemNumbered", "Subitem", "DisplayFormula",
   "DisplayFormulaNumbered", "Program"};
arCodeStyles = {"Input", "Code"};
arCollect[nb_Notebook] := Module[{k = 0, j = 0, out = {}, walk},
   walk[Cell[CellGroupData[l_List, ___], ___]] := Scan[walk, l];
   walk[Cell[content_, style_String, opts___]] := Which[
     MemberQ[arTextStyles, style], k++,
     MemberQ[arCodeStyles, style], k++; j++;
      AppendTo[out, <|"k" -> k, "j" -> j,
        "label" -> Replace[Quiet@Lookup[{opts}, CellLabel, None], Except[_String] -> ""],
        "content" -> content|>],
     True, Null];
   walk[_] := Null;
   Scan[walk, First[nb]];
   out];
arCellId[c_] := "#" <> ToString[c["k"]] <> " in#" <> ToString[c["j"]] <>
   If[c["label"] =!= "", " " <> StringTrim[c["label"]], ""];

(* ---------- number formatting ---------- *)
arFmtNZ[q_] := Module[{e, m},
   e = Floor[Log10[Abs[N[q, 30]]]];
   m = Round[q*10^(39 - e)];
   If[Abs[m] >= 10^40, e = e + 1; m = Round[q*10^(39 - e)]];
   If[Abs[m] < 10^39, e = e - 1; m = Round[q*10^(39 - e)]];
   ToString[m] <> "e" <> ToString[e - 39]];
arFmtReal[r_] := Which[
   r == 0, If[! ExactNumberQ[r], arNApproxZero++]; "0",
   ExactNumberQ[r], arFmtNZ[r],
   MachineNumberQ[r], arNMachine++; arFmtNZ[SetPrecision[r, Infinity]],
   True, arMinPrec = Min[arMinPrec, Precision[r]]; arFmtNZ[SetPrecision[r, Infinity]]];
arFmt[v_] := Which[
   v === 0, "0",
   Head[v] === Complex,
    arNComplex++;
    If[! (Im[v] == 0), arNComplexBad++;
     arMaxImRel = Max[arMaxImRel, N[Abs[Im[v]]/Max[Abs[Re[v]], 10^-300]]]];
    arFmtReal[Re[v]],
   True, arFmtReal[v]];

arWrite[dir_String, name_String, mat_] := Module[{txt, gz, st, rp, t0 = AbsoluteTime[]},
   arNApproxZero = 0; arNMachine = 0; arNComplex = 0; arNComplexBad = 0;
   arMaxImRel = 0; arMinPrec = Infinity;
   txt = FileNameJoin[{dir, name <> ".txt"}]; gz = txt <> ".gz";
   st = OpenWrite[txt];
   If[ArrayDepth[mat] == 1,
    WriteString[st, StringRiffle[arFmt /@ mat, "\n"], "\n"],
    Scan[WriteString[st, StringRiffle[arFmt /@ #, "\n"], "\n"] &, mat]];
   Close[st];
   If[FileExistsQ[gz], DeleteFile[gz]];
   rp = RunProcess[{"gzip", "-n", "-6", txt}];
   If[rp["ExitCode"] =!= 0, arFail["gzip failed for " <> name <> ": " <> rp["StandardError"]]];
   <|"shape" -> Dimensions[mat], "file" -> name <> ".txt.gz",
     "sha256" -> FileHash[gz, "SHA256", All, "HexString"],
     "stats" -> <|"approx_zero_written_as_0" -> arNApproxZero, "machine_entries" -> arNMachine,
       "complex_entries" -> arNComplex, "complex_nonzero_imag" -> arNComplexBad,
       "max_rel_imag" -> N[arMaxImRel],
       "min_precision_digits" -> If[arMinPrec === Infinity, Null, N[arMinPrec]],
       "seconds" -> Round[AbsoluteTime[] - t0, 0.01]|>|>];

(* ---------- main ---------- *)
arMap = {"x" -> "vx", "xp" -> "vdx\[Phi]", "K" -> "K", "K0" -> "K0", "Lam" -> "\[CapitalLambda]I",
   "Re_sig" -> "a", "Im_sig" -> "b", "Re_rho" -> "A", "Im_rho" -> "B",
   "fchi_sig" -> "fsigmachi", "fchi_rho" -> "frhochi", "fF_sig" -> "fsigmaF", "fF_rho" -> "frhoF",
   "Lambda3" -> "\[CapitalLambda]3", "KFF3" -> "KFF3", "BSpline3" -> "BSplineM3",
   "intS0" -> "vintS0", "intP1" -> "vintP1", "intD0" -> "vintD0", "SVZSR" -> "SVZSR"};

arMain[] := Module[{nbFile, outDir, mOv, lOv, commit, srcRel, nbHash, nb, cells, held, ids,
    skip, mCell, recs = {}, c, t, msgs, tEval, tEval0, tExp0, tExp, mm, ll, n0, nbb, mrho,
    expect, mats = <||>, shapeBad = {}, man, mod},
   nbFile = arOpt["nb", "/nb/GTB_numerics.nb"];
   outDir = arOpt["out", None];
   If[outDir === None, Print["out=DIR required"]; Exit[2]];
   mOv = arIntOpt["M"]; lOv = arIntOpt["L"];
   commit = arOpt["commit", "801684d9ece3de2918a20145178f569459081098"];
   srcRel = arOpt["source", "third_party/gauge-theory-bootstrap/theories/qcd/papers/arxiv-2403.10772/GTB_numerics.nb"];
   If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
   If[FileExistsQ[FileNameJoin[{outDir, "manifest.json"}]] && ! MemberQ[arArgs, "force"],
    Print["manifest.json exists in ", outDir, " (pass force)"]; Exit[2]];
   arLogFile = FileNameJoin[{outDir, "run_log.txt"}];
   If[FileExistsQ[arLogFile], DeleteFile[arLogFile]];
   arLog["authors_reference.wls  " <> $Version <> "  args " <> ToString[arArgs, InputForm]];
   nbHash = FileHash[nbFile, "SHA256", All, "HexString"];
   arLog["notebook " <> nbFile <> " sha256 " <> nbHash];
   nb = Get[nbFile];
   If[Head[nb] =!= Notebook, arFail["not a notebook"]];
   cells = arCollect[nb];
   held = Block[{$Context = "Global`"}, Quiet@Check[ToExpression[#["content"], StandardForm, HoldComplete], $Failed]] & /@ cells;
   If[MemberQ[held, $Failed], arFail["unparsable cell(s): " <> ToString[arCellId /@ Pick[cells, (# === $Failed) & /@ held]]]];
   ids = arCellId /@ cells;
   skip = Join[Flatten@Position[held, h_ /; ! FreeQ[h, SetDirectory], {1}, Heads -> False],
     Flatten@Position[held, h_ /; ! FreeQ[h, Export], {1}, Heads -> False]];
   If[Length[skip] =!= 2 || Last[skip] =!= Length[cells],
    arFail["expected one SetDirectory cell and one final Export cell; got " <> ToString[ids[[skip]]]]];
   mCell = Flatten@Position[held, h_ /; StringContainsQ[ToString[h, InputForm], "M = 50"], {1}, Heads -> False];
   If[Length[mCell] =!= 1, arFail["expected one cell setting M = 50, got " <> ToString[ids[[mCell]]]]];
   mCell = First[mCell];
   arLog[ToString[Length[cells]] <> " code cells; skipping " <> StringRiffle[ids[[skip]], ", "] <>
     "; M-cell " <> ids[[mCell]]];

   (* ---- evaluation ---- *)
   tEval0 = AbsoluteTime[];
   Do[
    c = cells[[i]];
    If[MemberQ[skip, i],
     arLog["skip " <> ids[[i]]];
     AppendTo[recs, <|"cell" -> ids[[i]], "skipped" -> True|>]; Continue[]];
    {t, msgs} = AbsoluteTiming[Block[{$Context = "Global`", $MessageList = {}},
       ToExpression[c["content"], StandardForm]; $MessageList]];
    msgs = ToString[#, InputForm] & /@ msgs;
    arLog["cell " <> ids[[i]] <> "  " <> ToString[Round[t, 0.001]] <> " s  mem " <>
      ToString[Round[MemoryInUse[]/2.^30, 0.01]] <> " GiB" <>
      If[msgs === {}, "", "  MESSAGES " <> StringRiffle[msgs, " "]]];
    AppendTo[recs, <|"cell" -> ids[[i]], "seconds" -> Round[t, 0.001],
      "mem_GiB" -> Round[MemoryInUse[]/2.^30, 0.01], "messages" -> msgs|>];
    If[i === mCell && (mOv =!= None || lOv =!= None),
     If[mOv =!= None, Block[{$Context = "Global`"}, ToExpression["M = " <> ToString[mOv]]]];
     If[lOv =!= None, Block[{$Context = "Global`"}, ToExpression["L = " <> ToString[lOv]]]];
     arLog["OVERRIDE after " <> ids[[i]] <> ": M = " <> ToString[arNb["M"]] <> ", L = " <> ToString[arNb["L"]]];
     AppendTo[recs, <|"cell" -> "override", "M" -> Replace[mOv, None -> Null], "L" -> Replace[lOv, None -> Null]|>]],
    {i, Length[cells]}];
   tEval = AbsoluteTime[] - tEval0;
   arLog["evaluation wall time " <> ToString[Round[tEval, 0.1]] <> " s, max memory " <>
     ToString[Round[MaxMemoryUsed[]/2.^30, 0.01]] <> " GiB"];

   (* ---- export ---- *)
   tExp0 = AbsoluteTime[];
   {mm, ll, n0, nbb} = arNb /@ {"M", "L", "n0", "Nb"};
   mrho = mm^2 + mm (mm + 1)/2;
   expect = <|"x" -> {mm}, "xp" -> {mm}, "K" -> {mm, mm}, "K0" -> {mm}, "Lam" -> {3 mm ll/2},
     "Re_sig" -> {3 mm ll/2, 2 mm + 1}, "Im_sig" -> {3 mm ll/2, 2 mm + 1},
     "Re_rho" -> {3 mm ll/2, mrho}, "Im_rho" -> {3 mm ll/2, mrho},
     "fchi_sig" -> {3 (ll/2) 4, 2 mm + 1}, "fchi_rho" -> {3 (ll/2) 4, mrho},
     "fF_sig" -> {3 ll/2, 2 mm + 1}, "fF_rho" -> {3 ll/2, mrho}, "Lambda3" -> {3 n0},
     "KFF3" -> {3 mm, 3 mm}, "BSpline3" -> {3 n0, 3 nbb}, "intS0" -> {3, n0},
     "intP1" -> {3, n0}, "intD0" -> {3, n0}, "SVZSR" -> {9}|>;
   Do[
    With[{name = r[[1]], val = arNb[r[[2]]]},
     If[! ArrayQ[val, _, NumberQ], arFail[name <> " <- " <> r[[2]] <> " is not a numeric array"]];
     If[Dimensions[val] =!= expect[name], AppendTo[shapeBad, name -> {Dimensions[val], expect[name]}]];
     mats[name] = arWrite[outDir, name, val];
     arLog["wrote " <> name <> " " <> ToString[Dimensions[val]] <> "  " <> ToString[mats[name]["stats"], InputForm]]],
    {r, arMap}];
   tExp = AbsoluteTime[] - tExp0;
   arLog["export wall time " <> ToString[Round[tExp, 0.1]] <> " s"];
   If[shapeBad =!= {}, arLog["SHAPE MISMATCH " <> ToString[shapeBad, InputForm]]];

   mod = If[mOv === None && lOv === None, "none",
     StringRiffle[Join[If[mOv =!= None, {"M = " <> ToString[mOv]}, {}], If[lOv =!= None, {"L = " <> ToString[lOv]}, {}]], "; "] <>
      " evaluated right after notebook cell " <> ids[[mCell]] <> " (which sets M = 50; L = 20); no other change"];
   man = <|"M" -> mm, "L" -> ll, "nu0" -> arNb["\[Nu]0"], "digits" -> 40, "precision" -> arNb["p"],
     "generator" -> "authors",
     "source" -> <|"file" -> srcRel, "sha256" -> nbHash, "commit" -> commit|>,
     "matrices" -> (KeyDrop[#, "stats"] & /@ mats),
     "n0" -> n0,
     "modification" -> mod,
     "skipped_cells" -> ids[[skip]],
     "machine_precision_matrices" -> Keys@Select[mats, #["stats"]["machine_entries"] > 0 &],
     "complex_entries_note" -> "LegendreQ[l,0,3,z] gives Complex values with an approximate-zero imaginary part; such entries are written as their real part (counts per matrix in entry_stats; complex_nonzero_imag counts imaginary parts that are not == 0)",
     "entry_stats" -> (#["stats"] & /@ mats),
     "shape_check" -> If[shapeBad === {}, "ok", ToString[shapeBad, InputForm]],
     "timing_s" -> <|"evaluate" -> Round[tEval, 0.1], "export" -> Round[tExp, 0.1]|>,
     "max_memory_GiB" -> Round[MaxMemoryUsed[]/2.^30, 0.01],
     "wolfram_version" -> $Version|>;
   Export[FileNameJoin[{outDir, "timings.json"}], <|"evaluate_s" -> Round[tEval, 0.1],
      "export_s" -> Round[tExp, 0.1], "cells" -> recs|>, "RawJSON"];
   Export[FileNameJoin[{outDir, "manifest.json"}], man, "RawJSON"];
   arLog["done: " <> FileNameJoin[{outDir, "manifest.json"}]];
   ];

arMain[];
End[];
Exit[0];
