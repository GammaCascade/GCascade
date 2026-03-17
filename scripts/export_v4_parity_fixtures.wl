(*
Export expected parity outputs from GCascadeV4 for pre-generated fixture inputs.

Usage:
  wolframscript -file scripts/export_v4_parity_fixtures.wl \
    /path/to/GCascadeV5/benchmarks/v4_reference \
    /path/to/GCascade

Arguments:
  1) fixtures root (contains case dirs with meta.json and input CSV files)
  2) GCascadeV4 package root (contains GCascadeV4.wl and LibrariesV4)
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

runFunction[fnName_, inj_, zStart_] := Switch[fnName,
  "RedshiftPoint", RedshiftPoint[inj, zStart],
  "AttenuatePoint", AttenuatePoint[inj, zStart],
  "CascadePoint", CascadePoint[inj, zStart],
  _, (Print["Unsupported point function: ", fnName]; Abort[])
];

runDiffuseFunction[fnName_, inj_, zStart_, zDistrib_] := Switch[fnName,
  "RedshiftDiffuse", RedshiftDiffuse[inj, zStart, zDistrib],
  "AttenuateDiffuse", AttenuateDiffuse[inj, zStart, zDistrib],
  "CascadeDiffuse", CascadeDiffuse[inj, zStart, zDistrib],
  _, (Print["Unsupported diffuse function: ", fnName]; Abort[])
];

runEvolvingFunction[fnName_, inj2d_, zStart_, zDistrib_] := Switch[fnName,
  "RedshiftEvolving", RedshiftEvolving[inj2d, zStart, zDistrib],
  "AttenuateEvolving", AttenuateEvolving[inj2d, zStart, zDistrib],
  "CascadeEvolving", CascadeEvolving[inj2d, zStart, zDistrib],
  _, (Print["Unsupported evolving function: ", fnName]; Abort[])
];

applyAction[actionAssoc_] := Module[{actionName, actionArgs},
  actionName = actionAssoc["action"];
  actionArgs = actionAssoc["args"];
  Switch[actionName,
    "changeEBLModel", changeEBLModel[Round[actionArgs[[1]]]],
    "changeMagneticField", changeMagneticField[actionArgs[[1]], actionArgs[[2]], Round[actionArgs[[3]]]],
    _, (Print["Unknown pre_action: ", actionName]; Abort[])
  ];
];

exportCycleSparse[dir_, sparseIndices_] := Module[{zIdx, eInIdx, eOutIdx, rows},
  zIdx = sparseIndices["z_idx"];
  eInIdx = sparseIndices["e_in_idx"];
  eOutIdx = sparseIndices["e_out_idx"];
  rows = Flatten[
    Table[
      {
        z, eIn, eOut,
        cycleSpec[[z + 1, eIn + 1, eOut + 1]]
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
caseDirs = Sort[caseDirs];
If[StringLength[caseFilter] > 0, caseDirs = Select[caseDirs, StringContainsQ[FileNameTake[#], caseFilter] &]];
If[caseLimit > 0, caseDirs = Take[caseDirs, UpTo[caseLimit]]];

If[Length[caseDirs] == 0,
  Print["No fixture case directories found under ", fixturesRoot];
  Exit[0];
];

SetDirectory[v4Base];

caseCounter = 0;
For[i = 1, i <= Length[caseDirs], i++,
  dir = caseDirs[[i]];
  meta = loadMeta[dir];
  fnName = meta["function"];
  zStart = N[meta["z_start"]];
  inputsKind = meta["inputs_kind"];
  preActions = If[KeyExistsQ[meta, "pre_actions"], meta["pre_actions"], {}];
  parityTargets = If[KeyExistsQ[meta, "parity_targets"], meta["parity_targets"], {"output"}];

  Get[FileNameJoin[{v4Base, "GCascadeV4.wl"}]];

  If[Length[preActions] > 0,
    Do[applyAction[preActions[[j]]], {j, 1, Length[preActions]}];
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
  If[Mod[caseCounter, 10] == 0 || caseCounter == Length[caseDirs],
    Print["Processed ", caseCounter, "/", Length[caseDirs], " cases..."];
  ];
];

Print["Completed export for ", caseCounter, " fixture cases."];
