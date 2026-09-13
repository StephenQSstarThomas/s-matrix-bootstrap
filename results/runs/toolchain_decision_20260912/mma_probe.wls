Export["/out/mma_probe.json", <|
  "MachineID" -> $MachineID,
  "Version" -> $Version,
  "LicenseExpirationDate" -> ToString[$LicenseExpirationDate],
  "LicenseType" -> ToString[$LicenseType],
  "ProcessorCount" -> $ProcessorCount,
  "symbolic_integrate_1_over_x_minus_3_4_to_10" -> ToString[N[Integrate[1/(x - 3), {x, 4, 10}], 20]],
  "LegendreQ_2_of_3p5" -> ToString[N[LegendreQ[2, 7/2], 20]],
  "pi_40_digits" -> ToString[N[Pi, 40]],
  "author_gen_matrices_input_cells" -> Length[Cases[
      Get["/nb/src/mathematica/gtb_qcd_01_generate_matrices.nb"],
      Cell[c_, "Input", ___] :> c, Infinity]],
  "author_2403_input_cells" -> Length[Cases[
      Get["/nb/papers/arxiv-2403.10772/GTB_numerics.nb"],
      Cell[c_, "Input", ___] :> c, Infinity]]
|>, "RawJSON"];
Print["probe ok"];
