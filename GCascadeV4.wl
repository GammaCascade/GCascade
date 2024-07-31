(* ::Package:: *)

BeginPackage["GCascadeV4`"]


energies::usage = "An array of 300 logarithmically-spaced energies from \!\(\*SuperscriptBox[\(10\), \(-1\)]\) to \!\(\*SuperscriptBox[\(10\), \(12\)]\) GeV, which serve as the energies of gamma rays.";
specPlot::usage = "specPlot[]: specPlot is a function which plots \!\(\*SuperscriptBox[\(E\), \(2\)]\)f(E) given a grid with some function f(E) evaluated at the values of the energies array.";
cutoffPowerLaw::usage = "cutoffPowerLaw[energy E, spectral index \[CapitalGamma], cutoff energy \!\(\*SubscriptBox[\(E\), \(cut\)]\), normalization N]: Function that generates spectrum values of a broken powerlaw \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) = N \[Times] \!\(\*SuperscriptBox[\(E\), \(\[CapitalGamma]\)]\) \[Times] Exp[-E/\!\(\*SubscriptBox[\(E\), \(cut\)]\)] for a given energy, spectral index, cutoff energy, and normalization.";
(*zReg::usage = "Array which holds the values of z which will be treated as the regions of constant CMB/EBL. These range from 0 to 10 in steps of 0.01";*)
diffuseDistances::usage = "Array of redshift values for defining source distribution grids and for cosmological redshifting of spectra during cascade evaluation.";
(*diffuseSteps::usage = "Array which holds values of \[CapitalDelta]z which will be treated as the limits of cascade evaluation";*)
EBLindex::usage = "Index referring to the current EBL model active in \[Gamma]-Cascade.";


RedshiftPoint::usage = "RedshiftPoint[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], source redshift z]:
RedshiftPoint is a function which takes in a gamma-ray injected spectrum from a point source located at a given redshift, 
and produces the corresponding redshifted gamma-ray flux at the Earth, \[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\)].
It only takes into account cosmological redshifing effects during propagation.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
Source redshifts beyond z=10 are not supported by \[Gamma]-Cascade.";

AttenuatePoint::usage = "AttenuatePoint[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], source redshift z]:
AttenuatePoint is a function which takes in a gamma-ray injected spectrum from a point source located at a given redshift, 
and produces the corresponding attenuated gamma-ray flux at the Earth, \[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\)].
It takes into account cosmological redshifing effects and atteunuation of gamma rays due to pair production with the CMB/EBL.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
Source redshifts beyond z=10 are not supported by \[Gamma]-Cascade.";

CascadePoint::usage = "CascadePoint[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], source redshift z]:
CascadePoint is a function which takes in a gamma-ray injected spectrum from a point source located at a given redshift, 
and produces the corresponding total (attenuated + cascaded) gamma-ray flux at the Earth, \[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\)], after an electromagnetic cascade.
It takes into account cosmological redshifing effects, atteunuation due to pair production, and regeneration of gamma rays
by inverse Compton scattering with the CMB/EBL. It also accounts for synchrotron energy loss due to intergalactic magnetic fields.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
Source redshifts beyond z=10 are not supported by \[Gamma]-Cascade.";

RedshiftDiffuse::usage = "RedshiftDiffuse[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], maximum source redshift \!\(\*SubscriptBox[\(z\), \(max\)]\), redshift distribution (density per comobing volume) of sources \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)]]:
RedshiftDiffuse is a function which takes in a gamma-ray injected spectrum from a population of identical sources, 
with intrinsic spectra \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], following a given comoving density distribution in redshift, \[Rho](z) = \!\(\*FractionBox[SubscriptBox[\(dN\), \(sources\)], SubscriptBox[\(dV\), \(com\)]]\)(z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)], 
up to a redshift of \!\(\*SubscriptBox[\(z\), \(max\)]\), and produces the corresponding redshifted gamma-ray flux at the Earth, \[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\) \!\(\*SuperscriptBox[\(sr\), \(-1\)]\)].
It only takes into account cosmological redshifing effects during propagation.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
The density distribution should be formatted as an array of values \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)] evaluated at the redshifts given by the diffuseDistances array.
Distributions with \!\(\*SubscriptBox[\(z\), \(max\)]\) < 10 are not supported by \[Gamma]-Cascade.";

AttenuateDiffuse::usage = "AttenuateDiffuse[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], maximum source redshift \!\(\*SubscriptBox[\(z\), \(max\)]\), redshift distribution (density per comobing volume) of sources \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)]]:
AttenuateDiffuse is a function which takes in a gamma-ray injected spectrum from a population of identical sources, 
with intrinsic spectra \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], following a given comoving density distribution in redshift, \[Rho](z) = \!\(\*FractionBox[SubscriptBox[\(dN\), \(sources\)], SubscriptBox[\(dV\), \(com\)]]\)(z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)], 
up to a redshift of \!\(\*SubscriptBox[\(z\), \(max\)]\), and produces the corresponding attenuated gamma-ray flux at the Earth, \[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\) \!\(\*SuperscriptBox[\(sr\), \(-1\)]\)].
It takes into account cosmological redshifing effects and atteunuation of gamma rays due to pair production with the CMB/EBL.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
The density distribution should be formatted as an array of values \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)] evaluated at the redshifts given by the diffuseDistances array.
Distributions with \!\(\*SubscriptBox[\(z\), \(max\)]\) < 10 are not supported by \[Gamma]-Cascade.";

CascadeDiffuse::usage = "CascadeDiffuse[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], maximum source redshift \!\(\*SubscriptBox[\(z\), \(max\)]\), redshift distribution (density per comobing volume) of sources \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)]]:
CascadeDiffuse is a function which takes in a gamma-ray injected spectrum from a population of identical sources, 
with intrinsic spectra \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], following a given comoving density distribution in redshift, \[Rho](z) = \!\(\*FractionBox[SubscriptBox[\(dN\), \(sources\)], SubscriptBox[\(dV\), \(com\)]]\)(z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)], 
up to a redshift of \!\(\*SubscriptBox[\(z\), \(max\)]\), and produces the corresponding total (attenuated + cascaded) isotropic gamma-ray flux at the Earth, 
\[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\) \!\(\*SuperscriptBox[\(sr\), \(-1\)]\)], after electromagnetic cascades have taken place.
It takes into account cosmological redshifing effects, atteunuation due to pair production, and regeneration of gamma rays
by inverse Compton scattering with the CMB/EBL. It also accounts for synchrotron energy loss due to intergalactic magnetic fields.
The injected spectrum should be formatted as an array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the energies given by energies array.
The density distribution should be formatted as an array of values \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)] evaluated at the redshifts given by the diffuseDistances array.
Distributions with \!\(\*SubscriptBox[\(z\), \(max\)]\) < 10 are not supported by \[Gamma]-Cascade.";

CascadeEvolving::usage = "CascadeEvolving[injected spectrum \!\(\*FractionBox[\(dN\), \(dE\)]\)(E,z) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], maximum source redshift \!\(\*SubscriptBox[\(z\), \(max\)]\), redshift distribution (density per comobing volume) of sources \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)]]:
CascadeEvolving is a function which takes in gamma-ray injected spectra from a population of evolving sources, 
that is, sources whose intrinsic spectra are not identical, but instead evolve with redshift as \!\(\*FractionBox[\(dN\), \(dE\)]\)(E,z) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)], 
which follow a given comoving density distribution in redshift, \[Rho](z) = \!\(\*FractionBox[SubscriptBox[\(dN\), \(sources\)], SubscriptBox[\(dV\), \(com\)]]\)(z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)], up to a redshift of \!\(\*SubscriptBox[\(z\), \(max\)]\), 
and produces the corresponding total (attenuated + cascaded) isotropic gamma-ray flux at the Earth, 
\[CapitalPhi](E) [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\) \!\(\*SuperscriptBox[\(cm\), \(-2\)]\) \!\(\*SuperscriptBox[\(sr\), \(-1\)]\)], after electromagnetic cascades have taken place.
It takes into account cosmological redshifing effects, atteunuation due to pair production, and regeneration of gamma rays
by inverse Compton scattering with the CMB/EBL. It also accounts for synchrotron energy loss due to intergalactic magnetic fields.
The injected spectrum should be formatted as a 2D array of values \!\(\*FractionBox[\(dN\), \(dE\)]\)(E,z) [GeV \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] evaluated at the redshifts given by the diffuseDistances array and at the energies given by energies array.
The density distribution should be formatted as an array of values \[Rho](z) [\!\(\*SuperscriptBox[\(cm\), \(-3\)]\)] evaluated at the redshifts given by the diffuseDistances array.
Distributions with \!\(\*SubscriptBox[\(z\), \(max\)]\) < 10 are not supported by \[Gamma]-Cascade.";

changeMagneticField::usage = "changeMagneticField[BField, name] is a function which changes the intergalactic magnetic field of \[Gamma]-Cascade by exporting new libraries. See the tutorial notebook for complete instructions.";

changeEBLModel::usage = "changeEBLModel[new EBL index] is a function which changes the EBL model used in \[Gamma]-Cascade. The default model is Saldana-Lopez et al. 2021 (EBL index = 1).";


Begin["`Private`"]
Print["You have imported \[Gamma]-Cascade V4. Setup will take about a minute."]; 


(*Initialization: physical parameters*)

Mpc = 3.08568*10^24(*[cm/Mpc]*); c = 299792.458(*[km/s]*); t = 2.725(*CMB temp @ z=0 [K]*);
mu0 = 1.256637062*10^-6(*[N/A^2]*); elecmass = 0.51099895*10^6(*[eV]*); echarge = 1.602176634*10^-19(*[C]*);
sigmaTe = 6.65245873*10^-29(*[m^2]*);

(*Cosmological parameters from Planck [arXiv:1807.06209, A&A 641, A6 (2020)]*)
H0=67.4 (*[(km/s)/Mpc]*);
\[CapitalOmega]\[CapitalLambda]=0.685;
\[CapitalOmega]m=0.315;
hubble[z_]:=H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m*(1+z)^3](*[(km/s)/Mpc]*);



(*Initialization: useful functions and tables*)

(*Setup diffuseDistances array, holding redshift values from 10^-6 to 10 in irregular spacings \[CapitalDelta]z (10^-6 <= \[CapitalDelta]z <= 10^-2)*)
firstStep=Table[If[i==0,Table[10^-6*10^i,{j,1,10}],Table[10^-6*10^i,{j,1,9}]]//N,{i,0,4}]//Flatten ;
diffuseSteps={firstStep,ConstantArray[0.01,990]}//Flatten;
diffuseDistances=Table[Total[diffuseSteps[[;;i]]],{i,1,Length[diffuseSteps]}];

(*Finding the position in diffuseDistances closest to a given z, used to find position of Subscript[z, max] in the code*)
diffuseDistancesIndex[x_]:=Position[diffuseDistances,Nearest[diffuseDistances,x][[1]]][[1,1]];

(*Array which holds the values of z which will be treated as the regions of constant CMB/EBL*)
zReg=Range[0.,10.,0.01];

(*Generates an array of n logarithmically-spaced values from 10^a to 10^b*)
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)];

(*Generates an array of 300 logarithmically-spaced energies from 0.1 to 10^12 GeV. These serve as the energies of gamma rays*)
energies = logspace[-1,12,300](*[GeV]*);

(*Differences between the entries in energies (divided by 2, for trapezoidal integration), has length Length[energies]-1*)
dEnergiesGamma = Differences[energies]/2;

(*Function to generate bare plots of E^2dN/dE*)
specPlot[spec_]:=Module[{logSpec=Log10@(energies^2 spec),specInt},
specInt=Interpolation[Thread[{Log10@energies,logSpec}]];
Return[LogLogPlot[10.0^specInt[Log10[x]],{x,0.1,10^12}]];
];

(*Function that evaluates a broken power law for a given energy, spectral index, cutoff, and normalization*)
cutoffPowerLaw[enerGamma_,gamma_,cutoff_,amp_] := (amp*enerGamma^-gamma)*Exp[-enerGamma/cutoff]; 



(*Initialization: directories & arrays*)

packageName="GCascadeV4";
libraryLocation = FileNameJoin[{NotebookDirectory[],"LibrariesV4"}];

(*Default EBL model is Saldana-Lopez 2021*)
EBLindex=1;

(*Importing inverse mean free paths [cm^-1] for each photon background*)
IMFPcmb=Import[FileNameJoin[{libraryLocation,"pp-IMFPs","IMFPcmb.csv"}],"Table"];
IMFPebl=Import[FileNameJoin[{libraryLocation,"pp-IMFPs","IMFPeblSL.csv"}],"Table"];
IMFP=IMFPcmb+IMFPebl;

(*exp(-1/\[Lambda] [Mpc]) for PP attenuation, will later be taken to the power of lightDistance [Mpc] in each \[Delta]z (to become optical depths)*)
extinctionCoeffs=Table[Exp[-Mpc*IMFP[[j,i]]],{j,1,Length[zReg]},{i,1,Length[energies]}]//Quiet;

(*Photon spectrum after one PP + on-the-spot ICS cycle [GeV^-1]*)
cycleSpec=Import[FileNameJoin[{libraryLocation,"cycle-spec","cyclespecSL.mat"}]][[1]]*10^9;



(*Table of light-travel distances [Mpc], dimensions 1036 \[Times] 10000*) 
(*Entry i has distances from diffuseDistances[i] to diffuseDistances[i-1] redshift in steps of 10^-6, completed with zeros to become rectangular*)
stepSizeArray=Import[FileNameJoin[{libraryLocation,"stepSizeArray.mat"}]][[1]];
stepSizeArray=DeleteCases[stepSizeArray,0.,Infinity];

(*How to calculate stepSizeArray*)
(*lightDistanceTravelled[z1_,z2_]:=c*NIntegrate[1/(hubble[zVal]*(1+zVal)),{zVal,z1,z2},Method->{"Trapezoidal"},PrecisionGoal->2](*[Mpc]*);
lightDistancetable=ParallelTable[lightDistanceTravelled[z,z+10^-6],{z,0.,10.,10^-6}];
stepsizearray=Table[0,Length[diffuseDistances],Round[(diffuseDistances[[-1]]-diffuseDistances[[-2]])/10^-6]];
counter=0; For[i=1,i\[LessEqual]Length[stepsizearray],i++,jmax=Round[(diffuseDistances[[i]]-If[i==1,0,diffuseDistances[[i-1]]])/10^-6];counter=counter+jmax;For[j=1,j\[LessEqual]jmax,j++,stepsizearray[[i,j]]=lightDistancetable[[counter-(j-1)]]]];*)

(*Conversion of entries between zReg and diffuseDistances, dimensions 1036 \[Times] 10000*)
zRegIndexArray=Import[FileNameJoin[{libraryLocation,"zRegIndexArray.mat"}]][[1]];
zRegIndexArray=DeleteCases[zRegIndexArray,0.,Infinity];

(*How to calculate zRegIndexArray*)
(*zregindexarray=Table[0,Length[diffuseDistances],Round[(diffuseDistances[[-1]]-diffuseDistances[[-2]])/10^-6]];
For[i=1,i\[LessEqual]Length[zregindexarray],i++,For[j=1,j\[LessEqual]Length[zregindexarray[[1]]],j++,If[diffuseDistances[[i]]-If[i==1,0,diffuseDistances[[i-1]]]\[GreaterEqual]j*10^-6,zregindexarray[[i,j]]=Position[zReg,Nearest[zReg,diffuseDistances[[i]]-j*10^-6][[1]]][[1,1]]]]];*)


changeMagneticField[BField_(*[Gauss]*),gamma_,EBL_]:=Module[{
newB=BField/10000.(*[Tesla]*),
dEdtICS,
dEdtsync,
fICS,
EBLname,
ppspec,
OTSspec,
cyclespec=Table[0.,Length[zReg],Length[energies],Length[energies]],
changelog
},

If[MemberQ[{0,1,2,3,4,5,6},EBL]==False,Print["Invalid EBL index."];Abort[];];
If[EBL==0,EBLname="CMB",If[EBL==1,EBLname="SL",If[EBL==2,EBLname="SLhigh",If[EBL==3,EBLname="SLlow",If[EBL==4,EBLname="Finke",If[EBL==5,EBLname="Franc",EBLname="Dom"]]]]]];

Print["Changing default magnetic field in \[Gamma]-Cascade for chosen EBL model. This should take around 6 minutes."];

(*ICS and Synchrotron energy loss ratea [eV/s]*)
dEdtsync=Table[(1/echarge)*(4/3)*sigmaTe*((newB*(1+i)^gamma)^2/2/mu0)*c*1000*(Sqrt[(j*10^9)^2 - elecmass^2]/elecmass)^2,{i,zReg},{j,energies}];
dEdtICS=Import[FileNameJoin[{libraryLocation,"ics-E-loss-rates","dEdt"<>EBLname<>".mat"}]][[1]];
(*Calculating Subscript[f, ICS]*)
fICS=Table[1.,Length[zReg],Length[energies]];
fICS[[All,50;;]]=dEdtICS[[All,50;;]]/(dEdtsync[[All,50;;]]+dEdtICS[[All,50;;]]);

(*Renaming old spectral table*)
RenameFile[FileNameJoin[{libraryLocation,"cycle-spec","cyclespec"<>EBLname<>".mat"}],FileNameJoin[{libraryLocation,"cycle-spec","bak_cyclespec"<>EBLname<>".mat"}]];
Print["The name of the previous table containing PP-ICS cycle spectra for the chosen EBL model was prepended with ``bak_''. Rename it or keep it in a separate directory to avoid overwriting it the next time the magnetic field is changed."];

(*Importing PP and on-the-spot ICS spectral tables (multiplied by Subscript[f, ICS])*)
ppspec=Import[FileNameJoin[{libraryLocation,"pp-spec","normalizedPPspec"<>EBLname<>".mat"}]][[1]];
OTSspec=Import[FileNameJoin[{libraryLocation,"on-the-spot-ics-spec","onthespotICSspec"<>EBLname<>".mat"}]][[1]]*fICS;

(*Calculating and exporting new PP-ICS cycle spectral table*)
For[i=1,i<=Length[zReg],i++,For[j=1,j<=Length[energies],j++,cyclespec[[i,j]]=10^9*Total[dEnergiesGamma*(ppspec[[i,j,;;-2]]*OTSspec[[i,;;-2]]+ppspec[[i,j,2;;]]*OTSspec[[i,2;;]])]]];
Export[FileNameJoin[{libraryLocation,"cycle-spec","cyclespec"<>EBLname<>".mat"}],cyclespec];

(*Changing magnetic field (and EBL model, if not the current one) in the active kernel*)
If[EBL==EBLindex,cycleSpec=Import[FileNameJoin[{libraryLocation,"cycle-spec","cyclespec"<>EBLname<>".mat"}]][[1]]*10^9,changeEBLModel[EBL]];

(*Updating magnetic field change log*)
changelog=Import[FileNameJoin[{libraryLocation,"BfieldChangelog.txt"}],"List"];
AppendTo[changelog,"On "<>DateString[]<>", the intergalactic magnetic field for EBL model "<>ToString[EBL]<>" was changed.
The new field strength is B(z) = "<>ToString[BField//N,InputForm]<>"*(1+z)^"<>ToString[gamma]<>" Gauss."];
Export[FileNameJoin[{libraryLocation,"BfieldChangelog.txt"}],changelog];

Print["New magnetic field successfully implemented for the chosen EBL model. BfieldChangelog.txt has been updated."];

Return[Null];
];


changeEBLModel[EBL_]:=Block[{
oldEBL=If[EBLindex==0,"CMB only",If[EBLindex==1,"Saldana-Lopez et al. (2021)",If[EBLindex==2,"Saldana-Lopez et al. (2021, 1\[Sigma] upper uncertainty)",If[EBLindex==3,"Saldana-Lopez et al. (2021, 1\[Sigma] lower uncertainty)",If[EBLindex==4,"Finke et al. (2022)",If[EBLindex==5,"Franceschini & Rodighiero (2018)",If[EBLindex==6,"Dom\[IAcute]nguez et al. (2011)"]]]]]]],
newEBL=If[EBL==0,"CMB only",If[EBL==1,"Saldana-Lopez et al. (2021)",If[EBL==2,"Saldana-Lopez et al. (2021, 1\[Sigma] upper uncertainty)",If[EBL==3,"Saldana-Lopez et al. (2021, 1\[Sigma] lower uncertainty)",If[EBL==4,"Finke et al. (2022)",If[EBL==5,"Franceschini & Rodighiero (2018)",If[EBL==6,"Dom\[IAcute]nguez et al. (2011)",Print["Invalid EBL index."];Abort[];]]]]]]],
EBLname=If[EBL==0,"CMB",If[EBL==1,"SL",If[EBL==2,"SLhigh",If[EBL==3,"SLlow",If[EBL==4,"Finke",If[EBL==5,"Franc","Dom"]]]]]]
},

(*CMB only, EBLindex=0*)
(*SL, EBLindex=1    ->    default*)
(*SL high, EBLindex=2*)
(*SL low, EBLindex=3*)
(*Finke, EBLindex=4*)
(*Franc, EBLindex=5*)
(*Dom, EBLindex=6*)

If[EBL==EBLindex,Print[newEBL<>" is already the current EBL model."];Abort[];];

EBLindex=EBL;
If[EBL==0,IMFP=IMFPcmb,IMFPebl=Import[FileNameJoin[{libraryLocation,"pp-IMFPs","IMFPebl"<>EBLname<>".csv"}],"Table"];IMFP=IMFPcmb+IMFPebl;];
extinctionCoeffs=Table[Exp[-Mpc*IMFP[[j,i]]],{j,1,Length[zReg]},{i,1,Length[energies]}]//Quiet;
cycleSpec=Import[FileNameJoin[{libraryLocation,"cycle-spec","cyclespec"<>EBLname<>".mat"}]][[1]]*10^9;

Print["EBL model changed from "<>oldEBL<>" to "<>newEBL<>"."];

Return[Null];
];


RedshiftingCycle[injSpectra_,zArrayLocal_]:=Block[
{func,
 funcfix,
 logfunc,
 strechedResult,
 stretchedEnergies = energies*((1. + zArrayLocal[[1]])/(1. + zArrayLocal[[-1]]))
 },

(*Interpolation and redshifting*)
logfunc = ReplaceAll[Log10[injSpectra],{Indeterminate -> -200.}];
func = Interpolation[Thread[{energies, logfunc}]];
strechedResult = Append[Table[func[stretchedEnergies[[i]]],{i,1,Length[stretchedEnergies]-1}],func[energies[[-1]]]];(*careful: no redshifting of last energy bin*)
funcfix = Table[If[x>=-199.,10^x,0.],{x,strechedResult}];

Return[funcfix];
];


AttenuationCycle[injSpectra_,zArrayLocal_,stepSizeArrayLocal_,zRegIndexArrayLocal_]:=Block[
{func,
 funcfix,
 logfunc,
 strechedResult,
 finalResultLocal,
 paramsLocal,
 singleCycle,
 stretchedEnergies = energies*((1. + zArrayLocal[[1]])/(1. + zArrayLocal[[-1]]))
 },

(*Single cycle/step of size \[Delta]z = 10^-6*)
 singleCycle[spectrum_,zRegionIndex_,thisStepSize_]:=Block[{attenuatedSpec},         
        attenuatedSpec = (extinctionCoeffs[[zRegionIndex]]^thisStepSize)*spectrum;
		Return[attenuatedSpec];
];

paramsLocal=Table[{zRegIndexArrayLocal[[i]],stepSizeArrayLocal[[i]]},{i,1,Length[stepSizeArrayLocal]}];
finalResultLocal = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, paramsLocal];

(*Interpolation and redshifting*)
logfunc = ReplaceAll[Log10[finalResultLocal],{Indeterminate -> -200.}];
func = Interpolation[Thread[{energies, logfunc}]];
strechedResult = Append[Table[func[stretchedEnergies[[i]]],{i,1,Length[stretchedEnergies]-1}],func[energies[[-1]]]];
funcfix = Table[If[x>=-199.,10^x,0.],{x,strechedResult}];

Return[funcfix];
];


CascadeCycle[injSpectra_,zArrayLocal_,stepSizeArrayLocal_,zRegIndexArrayLocal_]:=Block[
{func,
 funcfix,
 logfunc,
 strechedResult,
 finalResultLocal,
 paramsLocal,
 singleCycle,
 stretchedEnergies = energies*((1. + zArrayLocal[[1]])/(1. + zArrayLocal[[-1]]))
 },

(*Single cycle/step of size \[Delta]z = 10^-6*)
 singleCycle[spectrum_,zRegionIndex_,thisStepSize_]:=Block[{attenuatedSpec,result,kern},
        attenuatedSpec = (extinctionCoeffs[[zRegionIndex]]^thisStepSize)*spectrum;
        kern=((spectrum - attenuatedSpec)*cycleSpec[[zRegionIndex]])[[2;;]]+((spectrum - attenuatedSpec)*cycleSpec[[zRegionIndex]])[[;;-2]]; (*trapezoidal integration*)
        result = Total[kern*dEnergiesGamma] + attenuatedSpec;
		Return[result];
];

paramsLocal=Table[{zRegIndexArrayLocal[[i]],stepSizeArrayLocal[[i]]},{i,1,Length[stepSizeArrayLocal]}];
finalResultLocal = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, paramsLocal];

(*Interpolation and redshifting*)
logfunc = ReplaceAll[Log10[finalResultLocal],{Indeterminate -> -200.}];
func = Interpolation[Thread[{energies, logfunc}]];
strechedResult = Append[Table[func[stretchedEnergies[[i]]],{i,1,Length[stretchedEnergies]-1}],func[energies[[-1]]]];
funcfix = Table[If[x>=-199.,10^x,0.],{x,strechedResult}];

Return[funcfix];
];


RedshiftPoint[injSpectraPre_,zStart_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6,
 dL
 },

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

dL=(1+zStart)*NIntegrate[c/hubble[zVal],{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2](*[Mpc]*);

(*Calculates an array where every member is a list of z-values, in stepsizes of \[Delta]z = stepSize, subdividing each \[CapitalDelta]z window in diffuseDistances*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,zMaxIndex-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

finalResult=((1+zStart)^2*Fold[RedshiftingCycle[#1,#2]&,injSpectra,Reverse[zArray]])/(4*\[Pi]*(dL*Mpc)^2)(*[GeV^-1 s^-1 cm^-2]*);

Return[finalResult];
Clear["zFunc$*",
"dL$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


AttenuatePoint[injSpectraPre_,zStart_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6,
 dL
 },

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

dL=(1+zStart)*NIntegrate[c/hubble[zVal],{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2](*[Mpc]*);

(*Calculates an array where every member is a list of z-values, in stepsizes of \[Delta]z = stepSize, subdividing each \[CapitalDelta]z window in diffuseDistances*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,zMaxIndex-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=((1+zStart)^2*Fold[AttenuationCycle[#1,#2[[1]],#2[[2]],#2[[3]]]&,injSpectra,Reverse[params]])/(4*\[Pi]*(dL*Mpc)^2)(*[GeV^-1 s^-1 cm^-2]*);

Return[finalResult];
Clear["zFunc$*",
"dL$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


CascadePoint[injSpectraPre_,zStart_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6,
 dL
 },

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

dL=(1+zStart)*NIntegrate[c/hubble[zVal],{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2](*[Mpc]*);

(*Calculates an array where every member is a list of z-values, in stepsizes of \[Delta]z = stepSize, subdividing each \[CapitalDelta]z window in diffuseDistances*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,zMaxIndex-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=((1+zStart)^2*Fold[CascadeCycle[#1,#2[[1]],#2[[2]],#2[[3]]]&,injSpectra,Reverse[params]])/(4*\[Pi]*(dL*Mpc)^2)(*[GeV^-1 s^-1 cm^-2]*);

Return[finalResult];
Clear["zFunc$*",
"dL$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


RedshiftDiffuse[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Comoving density distribution properly formatted."];,Print["Comoving density distribution is not formatted correctly. Please use a list of length "<>ToString[Length[diffuseDistances]]<>" containing values of {z,\[Rho](z)} for each redshhift in the ``diffuseDistances'' array."];Abort[];];

(*Interpolates \[Rho](z) [cm^-3] given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates \[Rho](z)Subscript[dV, c](z)/4Subscript[\[Pi]d, c](z) = \[Rho](z)dz/H(z), which is the redshift integrand*)
volumeNorms=Table[((c*Mpc*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/hubble[diffuseDistances[[i]]],{i,1,zMaxIndex}](*[cm^-2]*);

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]],zArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[RedshiftingCycle[#1+(#2[[1]]*injSpectra),#2[[2]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]]/(4*\[Pi])(*[GeV^-1 s^-1 cm^-2 sr^-1]*);

Return[finalResult];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


AttenuateDiffuse[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Comoving density distribution properly formatted."];,Print["Comoving density distribution is not formatted correctly. Please use a list of length "<>ToString[Length[diffuseDistances]]<>" containing values of {z,\[Rho](z)} for each redshhift in the ``diffuseDistances'' array."];Abort[];];

(*Interpolates \[Rho](z) [cm^-3] given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates \[Rho](z)Subscript[dV, c](z)/4Subscript[\[Pi]d, c](z) = \[Rho](z)dz/H(z), which is the redshift integrand*)
volumeNorms=Table[((c*Mpc*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/hubble[diffuseDistances[[i]]],{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]],zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[AttenuationCycle[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]]/(4*\[Pi])(*[GeV^-1 s^-1 cm^-2 sr^-1]*);

Return[finalResult];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


CascadeDiffuse[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energies],Print["Injected spectrum properly formatted."],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of dN/dE [\!\(\*SuperscriptBox[\(GeV\), \(-1\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\)] for your source for each energy in the ``energies'' array."];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Comoving density distribution properly formatted."];,Print["Comoving density distribution is not formatted correctly. Please use a list of length "<>ToString[Length[diffuseDistances]]<>" containing values of {z,\[Rho](z)} for each redshhift in the ``diffuseDistances'' array."];Abort[];];

(*Interpolates \[Rho](z) [cm^-3] given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates \[Rho](z)Subscript[dV, c](z)/4Subscript[\[Pi]d, c](z) = \[Rho](z)dz/H(z), which is the redshift integrand*)
volumeNorms=Table[((c*Mpc*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/hubble[diffuseDistances[[i]]],{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]],zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[CascadeCycle[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]]/(4*\[Pi])(*[GeV^-1 s^-1 cm^-2 sr^-1]*);

Return[finalResult];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


CascadeEvolving[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre(*[GeV^-1 s^-1]*),
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=diffuseDistancesIndex[zStart],
 stepSize=1.*10^-6
 },

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]=={Length[zFunction],Length[energies]},Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energies]]<>" containing values of luminosity as a function of energy for your source"];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Successfully incorporated luminosity function. GCascade will produce spectra for diffuse background generated by your injected spectra."];,Print["Luminosity function is not formatted correctly. Please use a list of length "<>ToString[Length[diffuseDistances]]<>" containing values of luminosity density as a function of z from z=0 to z=10"];Abort[];];

(*Interpolates \[Rho](z) given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates \[Rho](z)Subscript[dV, c](z)/4Subscript[\[Pi]d, c](z) = \[Rho](z)dz/H(z), which is the redshift integrand (without the redshift-dependent dN/dE)*)
volumeNorms=Table[((c*Mpc*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/hubble[diffuseDistances[[i]]],{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]]*injSpectra[[i]],zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[CascadeCycle[#1+#2[[1]],#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[energies]}],Reverse[params]]/(4*\[Pi])(*[GeV^-1 s^-1 cm^-2 sr^-1]*);

Return[finalResult];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


End[ ]

EndPackage[ ]
