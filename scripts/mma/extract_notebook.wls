(* Convert a notebook's Input cells to *evaluable* WL source.
   A cell holding "a; b;" comes back as Hold[a, b, Null]; emit each element as
   its own statement so the result can be Get[]. *)
args = Rest[$ScriptCommandLine];
If[Length[args] < 2, Print["usage: extract_notebook.wls IN.nb OUT.wl"]; Exit[2]];
{in, out} = args[[1 ;; 2]];
cells = Cases[Get[in], Cell[c_, "Input", ___] :> c, Infinity];
render[i_] := Module[{held, parts},
  held = Quiet@Check[ToExpression[cells[[i]], StandardForm, Hold], $Failed];
  If[held === $Failed, Return["(* cell " <> ToString[i] <> ": unconvertible *)"]];
  parts = List @@ Map[Hold, held];           (* one Hold per top-level expr *)
  StringRiffle[
   Table[StringReplace[ToString[p, InputForm],
      {StartOfString ~~ "Hold[" -> "", "]" ~~ EndOfString -> ""}], {p, parts}],
   ";\n"] <> ";"];
Export[out, StringRiffle[
   Table["(* ---- cell " <> ToString[i] <> " ---- *)\n" <> render[i],
         {i, Length[cells]}], "\n\n"], "Text"];
Print["cells: ", Length[cells], " -> ", out];
