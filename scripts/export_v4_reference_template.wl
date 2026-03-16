(*
Template for exporting V4 benchmark fixtures used by scripts/run_parity.py.
Run inside Mathematica/Wolfram Language after setting correct paths.
*)

baseV4 = "/path/to/GCascade";
outputBase = "/path/to/GCascadeV5/benchmarks/v4_reference";

SetDirectory[baseV4];
Get["GCascadeV4.wl"];

exportCase[caseName_, metaAssoc_, inj_, expected_, zDistrib_:None, inj2d_:None] := Module[
  {dir = FileNameJoin[{outputBase, caseName}]},
  CreateDirectory[dir, CreateIntermediateDirectories -> True];

  Export[FileNameJoin[{dir, "meta.json"}], ExportString[metaAssoc, "JSON"], "Text"];
  Export[FileNameJoin[{dir, "inj.csv"}], inj, "CSV"];
  Export[FileNameJoin[{dir, "expected.csv"}], expected, "CSV"];

  If[zDistrib =!= None, Export[FileNameJoin[{dir, "z_distrib.csv"}], zDistrib, "CSV"]];
  If[inj2d =!= None, Export[FileNameJoin[{dir, "inj2d.csv"}], inj2d, "CSV"]];
];

(* Example: RedshiftPoint case *)
inj = Table[cutoffPowerLaw[e, 2.2, 10^7, 10^40], {e, energies}];
z0 = 0.3;
expected0 = RedshiftPoint[inj, z0];

exportCase[
  "redshift_point_case_01",
  <|"function" -> "RedshiftPoint", "z_start" -> z0, "ebl_index" -> EBLindex|>,
  inj,
  expected0
];

(* Add additional canonical cases for each milestone before running Python parity checks. *)
