(*
Export expected parity outputs from GCascadeV4 for pre-generated fixture inputs.

Usage:
  wolframscript -file scripts/export_v4_parity_fixtures.wl \
    /path/to/GCascadeV5/benchmarks/v4_reference \
    /path/to/GCascade \
    [optional_case_filter_substring] \
    [optional_case_limit]

Arguments:
  1) fixtures root (contains case dirs with meta.json and input CSV files)
  2) GCascadeV4 package root (contains GCascadeV4.wl and LibrariesV4)
  3) optional case name substring filter
  4) optional positive integer limit
*)

args = Rest[$ScriptCommandLine];
fixturesRoot = If[Length[args] >= 1, args[[1]], FileNameJoin[{Directory[], "benchmarks", "v4_reference"}]];
v4Base = If[Length[args] >= 2, args[[2]], "/path/to/GCascade"];
caseFilter = If[Length[args] >= 3, args[[3]], ""];
caseLimit = If[Length[args] >= 4, ToExpression[args[[4]]], -1];

metaPath[dir_] := FileNameJoin[{dir, "meta.json"}];
expectedPath[dir_] := FileNameJoin[{dir, "expected.csv"}];
expectedSparsePath[dir_] := FileNameJoin[{dir, "expected_cycle_sparse.csv"}];

loadMeta[dir_] := Import[metaPath[dir], "RawJSON"];
toVector[x_] := Flatten[x];

(* Explicit GCascadeV4 context access to avoid Global` shadowing issues. *)
v4[name_String] := Symbol["GCascadeV4`" <> name];

groupPriority[group_] := Which[
  group === "propagation", 0,
  group === "change_ebl", 1,
  group === "change_bfield", 2,
  True, 3
];

safeChangeEBL[target_Integer] := Module[{current},
  current = v4["EBLindex"];
  If[target =!= current,
    v4["changeEBLModel"][target];
  ];
];

runFunction[fnName_, inj_, zStart_] := Switch[fnName,
  "RedshiftPoint", v4["RedshiftPoint"][inj, zStart],
  "AttenuatePoint", v4["AttenuatePoint"][inj, zStart],
  "CascadePoint", v4["CascadePoint"][inj, zStart],
  _, (Print["Unsupported point function: ", fnName]; Abort[])
];

runDiffuseFunction[fnName_, inj_, zStart_, zDistrib_] := Switch[fnName,
  "RedshiftDiffuse", v4["RedshiftDiffuse"][inj, zStart, zDistrib],
  "AttenuateDiffuse", v4["AttenuateDiffuse"][inj, zStart, zDistrib],
  "CascadeDiffuse", v4["CascadeDiffuse"][inj, zStart, zDistrib],
  _, (Print["Unsupported diffuse function: ", fnName]; Abort[])
];

runEvolvingFunction[fnName_, inj2d_, zStart_, zDistrib_] := Switch[fnName,
  "RedshiftEvolving", v4["RedshiftEvolving"][inj2d, zStart, zDistrib],
  "AttenuateEvolving", v4["AttenuateEvolving"][inj2d, zStart, zDistrib],
  "CascadeEvolving", v4["CascadeEvolving"][inj2d, zStart, zDistrib],
  _, (Print["Unsupported evolving function: ", fnName]; Abort[])
];

applyAction[actionAssoc_] := Module[{actionName, actionArgs},
  actionName = actionAssoc["action"];
  actionArgs = actionAssoc["args"];

  Switch[actionName,
    "changeEBLModel",
      safeChangeEBL[Round[actionArgs[[1]]]],
    "changeMagneticField",
      v4["changeMagneticField"][actionArgs[[1]], actionArgs[[2]], Round[actionArgs[[3]]]],
    _,
      (Print["Unknown pre_action: ", actionName]; Abort[])
  ];
];

exportCycleSparse[dir_, sparseIndices_] := Module[{zIdx, eInIdx, eOutIdx, cycle, rows},
  zIdx = sparseIndices["z_idx"];
  eInIdx = sparseIndices["e_in_idx"];
  eOutIdx = sparseIndices["e_out_idx"];
  cycle = v4["cycleSpec"];

  If[Depth[cycle] < 4,
    Print["cycleSpec is not rank-3 for case ", FileNameTake[dir], ". Depth=", Depth[cycle]];
    Abort[];
  ];

  rows = Flatten[
    Table[
      {
        z, eIn, eOut,
        cycle[[z + 1, eIn + 1, eOut + 1]]
      },
      {z, zIdx}, {eIn, eInIdx}, {eOut, eOutIdx}
    ],
    2
  ];

  Export[expectedSparsePath[dir], rows, "CSV"];
];

Print["Using fixtures root: ", fixturesRoot];
Print["Using V4 base: ", v4Base];
If[StringLength[caseFilter] > 0, Print["Applying case filter: ", caseFilter]];
If[caseLimit > 0, Print["Applying case limit: ", caseLimit]];

caseDirs = Select[
  FileNames["*", fixturesRoot],
  DirectoryQ[#] && FileExistsQ[metaPath[#]] &
];

If[StringLength[caseFilter] > 0,
  caseDirs = Select[caseDirs, StringContainsQ[FileNameTake[#], caseFilter] &];
];

casePairs = Table[{dir, loadMeta[dir]}, {dir, Sort[caseDirs]}];
casePairs = SortBy[
  casePairs,
  {
    groupPriority[If[KeyExistsQ[#[[2]], "case_group"], #[[2]]["case_group"], ""]] &,
    FileNameTake[#[[1]]] &
  }
];

If[caseLimit > 0,
  casePairs = Take[casePairs, UpTo[caseLimit]];
];

If[Length[casePairs] == 0,
  Print["No fixture case directories found under ", fixturesRoot];
  Exit[0];
];

SetDirectory[v4Base];
Get[FileNameJoin[{v4Base, "GCascadeV4.wl"}]];

Print["Loaded GCascadeV4 once. Exporting ", Length[casePairs], " cases..."];

caseCounter = 0;
For[i = 1, i <= Length[casePairs], i++,
  dir = casePairs[[i, 1]];
  meta = casePairs[[i, 2]];

  fnName = meta["function"];
  zStart = N[meta["z_start"]];
  inputsKind = meta["inputs_kind"];
  preActions = If[KeyExistsQ[meta, "pre_actions"], meta["pre_actions"], {}];
  parityTargets = If[KeyExistsQ[meta, "parity_targets"], meta["parity_targets"], {"output"}];
  targetEBL = If[KeyExistsQ[meta, "ebl_index"], Round[meta["ebl_index"]], 1];

  If[Length[preActions] > 0,
    Do[applyAction[preActions[[j]]], {j, 1, Length[preActions]}],
    safeChangeEBL[targetEBL]
  ];

  expected = Switch[inputsKind,
    "point",
      inj = toVector@Import[FileNameJoin[{dir, "inj.csv"}], "CSV"];
      runFunction[fnName, inj, zStart],
    "diffuse",
      inj = toVector@Import[FileNameJoin[{dir, "inj.csv"}], "CSV"];
      zDistrib = toVector@Import[FileNameJoin[{dir, "z_distrib.csv"}], "CSV"];
      runDiffuseFunction[fnName, inj, zStart, zDistrib],
    "evolving",
      inj2d = Import[FileNameJoin[{dir, "inj2d.csv"}], "CSV"];
      zDistrib = toVector@Import[FileNameJoin[{dir, "z_distrib.csv"}], "CSV"];
      runEvolvingFunction[fnName, inj2d, zStart, zDistrib],
    _,
      (Print["Unknown inputs_kind in ", dir, ": ", inputsKind]; Abort[])
  ];

  Export[expectedPath[dir], expected, "CSV"];

  If[MemberQ[parityTargets, "cycle_table_sparse"],
    If[KeyExistsQ[meta, "cycle_sparse_indices"],
      exportCycleSparse[dir, meta["cycle_sparse_indices"]],
      (Print["Missing cycle_sparse_indices for case ", dir]; Abort[])
    ]
  ];

  caseCounter++;
  If[Mod[caseCounter, 10] == 0 || caseCounter == Length[casePairs],
    Print["Processed ", caseCounter, "/", Length[casePairs], " cases..."];
  ];
];

Print["Completed export for ", caseCounter, " fixture cases."];
