(* ::Package:: *)

BeginPackage["GCascadeV3`"]



GCascadeAttenuate::usage = "CosmoAttenuate[injSpectra,zStart]: CAttenuate is a function which takes in a gamma-ray injected spectrum, starting z. This function produces an attenuated differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)] \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\) , without including cascade  evolution but taking into account cosmological energy shifts. The injected spectrum should be formatted as an array of values of differential flux \!\(\*FormBox[\((\),
TraditionalForm]\)\!\(\*FormBox[\(\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\)\),
TraditionalForm]\)) evaluated at the energies given by the array, energiesGamma
TraditionalForm]\)).";



energies::usage="a logarithmic array of 300 energies from 10^-1 to 10^12 GeV, these serve as the energies of gamma-rays";
specPlot::usage="Plots a primitive frame of \!\(\*SuperscriptBox[\(E\), \(2\)]\)\!\(\*FractionBox[\(dN\), \(dE\)]\) given some \!\(\*FractionBox[\(dN\), \(dE\)]\) as a single array evaluated at the values of energiesGamma";
cutoffPowerLaw::usage = "cutoffPowerLaw[energy (E),index (\[CapitalGamma]),cutoff (\!\(\*SubscriptBox[\(E\), \(cutoff\)]\)),normalization (N)]: Function that generates spectrum values of broken powerlaws (\!\(\*FractionBox[\(dN\), \(dE\)]\)= N*\!\(\*SuperscriptBox[\(E\), \(-\[CapitalGamma]\)]\)*Exp[-\!\(\*FractionBox[\(E\), SubscriptBox[\(E\), \(cutoss\)]]\)]) for a given energy, spectral index, cutoff, and normalization";
zReg::usage = "Array which holds the values of z which will be treated as the regions of constant EBL. These range from 0 to 10 in steps of 0.01";
diffuseDistances::usage="Array which holds values of Z which will be treated as the points of cascade evaluation";
diffuseSteps::usage="Array which holds values of \[CapitalDelta]Z which will be treated as the limits of cascade evaluation";


GCascadeDiffuse::usage = "GCascadeDiffuse[injected spectrum ,max Z (\!\(\*SubscriptBox[\(Z\), \(max\)]\)),Luminosity Function (L(z))]:GCascadeDiffuse is a function which takes in a gamma-ray injected luminosity spectrum, max z, and a luminosity function. 'GCascadeDiffuse' produces a cascaded diffuse differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)] \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)]\\\ \*SuperscriptBox[\(sr\), \(-1\)])\))\),
TraditionalForm]\) , taking into account cosmological and electromagnetic energy losses
TraditionalForm]\). The injected differential spectrum \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)])\),
TraditionalForm]\) should be formatted as an array of values of power per unit area \!\(\*FormBox[\((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\),
TraditionalForm]\) \!\(\*FormBox[SuperscriptBox[\(s\), \(-1\)],
TraditionalForm]\)) evaluated at the energies given by the array, 'energiesGamma'.  The Luminosity distribution function is assumed to be a separable function of the form: \!\(\*FormBox[\(\*FractionBox[\(dN\), \(\*SubscriptBox[\(dV\), \(c\)] dE\)] \((\*SubscriptBox[\(E\), \(\[Gamma]\)], z)\) = \*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SubscriptBox[\(E\), \(\[Gamma]\)])\) \(\[Rho](z)\)\),
TraditionalForm]\) where the left-hand factor is the intrinsic differential spectrum of the class of sources \!\(\*FormBox[\((\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\) and \[Rho](z) is the spacial distribution of the sources in units of density per comoving volume, so called 'luminosity function'. 'GCascadeDiffuse' uses the luminosity function as an input which must be formatted as list of ordered pairs (z,\!\(\*FormBox[\(\[Rho](z)\),
TraditionalForm]\)), of values of z and comoving density \[Rho](z). The values of z can be any finite list of z from 0 to 10. Additionally, GCascadeDiffuse uses the evolving differential spectrum, injected spectrum = \!\(\*FormBox[\(\*FractionBox[\(dN\), \(\*SubscriptBox[\(dV\), \(c\)] dE\)] \((\*SubscriptBox[\(E\), \(\[Gamma]\)], z)\)\),
TraditionalForm]\), as an input. injected spectrum must be a list of spectra, each formatted as in the examples above. Each member of, injected spectrum, must be the spectrum injected at the z-values given by, diffuseDistances. ";

GCascadeDiffuseConstant::usage = "GCascadeDiffuse[injected spectrum ,max Z (\!\(\*SubscriptBox[\(Z\), \(max\)]\))),Luminosity Function (L(z))]:GCascadeDiffuse is a function which takes in a gamma-ray injected luminosity spectrum, max z, step size in z, and a luminosity function. 'GCascadeDiffuse' produces a cascaded diffuse differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)] \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)]\\\ \*SuperscriptBox[\(sr\), \(-1\)])\))\),
TraditionalForm]\) , taking into account cosmological and electromagnetic energy losses
TraditionalForm]\). The injected differential spectrum \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)])\),
TraditionalForm]\) should be formatted as an array of values of power per unit area \!\(\*FormBox[\((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\),
TraditionalForm]\) \!\(\*FormBox[SuperscriptBox[\(s\), \(-1\)],
TraditionalForm]\)) evaluated at the energies given by the array, 'energiesGamma'.  The Luminosity distribution function is assumed to be a separable function of the form: \!\(\*FormBox[\(\*FractionBox[\(dN\), \(\*SubscriptBox[\(dV\), \(c\)] dE\)] \((\*SubscriptBox[\(E\), \(\[Gamma]\)], z)\) = \*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SubscriptBox[\(E\), \(\[Gamma]\)])\) \(\[Rho](z)\)\),
TraditionalForm]\) where the left-hand factor is the intrinsic differential spectrum of the class of sources \!\(\*FormBox[\((\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\) and \[Rho](z) is the spacial distribution of the sources in units of density per comoving volume, so called 'luminosity function'. 'GCascadeDiffuse' uses the luminosity function as an input which must be formatted as list of ordered pairs (z,\!\(\*FormBox[\(\[Rho](z)\),
TraditionalForm]\)), of values of z and comoving density \[Rho](z). The values of z can be any finite list of z from 0 to 10.";


GCascadePoint::usage = "GCascadePoint[injected spectrum, source distance (Z)]: GCascadePoint is a function which takes in a gamma-ray injected spectrum, source distance, and step size in Z. 'GCascadePoint' produces a cascaded differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\), taking into account cosmological and electromagnetic energy losses
TraditionalForm]\). The injected differential spectrum \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)])\),
TraditionalForm]\) should be formatted as an array of values of power per unit area \!\(\*FormBox[\((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\),
TraditionalForm]\) \!\(\*FormBox[SuperscriptBox[\(s\), \(-1\)],
TraditionalForm]\)) evaluated at the energies given by the array, 'energiesGamma'.";

changeMagneticField::usage = "changeMagneticField[BField_Rational,name_String] is a function which changes the IG magnetic field of GCascade by exporting new libraries. See tutorial for complete instructions.";



Begin["`Private`"]
Print["You Have Imported GCascade. Setup will take about a minute."]; 


libraryLocation = FileNameJoin[{NotebookDirectory[],"LibrariesV3"}];
(*libraryLocation = FileNameJoin[{"/home/carlos/Software/GCascade","LibrariesV2"}];*)
libraryNames="DiffCascadeSpectralTable";

kpc = 3.086*10^19; h = 4.135667/10.0^15; joule2eV = 6.242*10.0^18.0; c = 299792458.0; 
k = 8.6173324/10^5; 
t = 2.725; electronMass = 510998.; 
classicalRadius = 2.81794/10^15; 
thomp = (8./3)*Pi*classicalRadius^2; 
\[Mu] = 1.2566/10^6; interGalacticMagneticField = 1./10^13/10000; 


(*Cosmological parameters from Planck '15*)

packageName="GCascadeV3";
H0=67.7 (*(km/s)/Mpc*);
\[CapitalOmega]\[CapitalLambda]=0.691;
\[CapitalOmega]m=0.308;
 
(*Function which calculates the distance light travels from zCurrent to some redshift zCurrent-dz*)
(*lightDistanceTravelled[zCurrent_,dz_]:=NIntegrate[c/(H0(1+zVal)Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,zCurrent-dz,zCurrent},Method->{"Trapezoidal"},PrecisionGoal->2];*)
lightDistanceTravelled[z1_,z2_]:=NIntegrate[c/(H0(1+zVal)Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,z1,z2},Method->{"Trapezoidal"},PrecisionGoal->2];

(*Black Body differential Number Density*)

dndEtar[Etar_,t_] := ((8.0*Pi)/(h*c)^3)*(Etar^2/(Exp[Etar/(k*t)]-1.))



(*Generates an array of n logarithmically spaced values from 10^a to 10^b*)
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)]

(*Generates a logarithmic array of 300 energies from 10^8 to 10^20 eV, these serve as the energies of gamma-rays *)
energiesGamma = logspace[8,21,300];
energies = 10^-9*energiesGamma;

(*Function to calculate the spacing between energies in the gamma-ray energy array*)
deltaGamma[enerGamma_] := enerGamma*((10.0^((21.-8.)/(300.-1.))) - 1);

(*Generates an array which holds the logarithmic spacing in the gamma-ray energy array, i.e. an array of values Subscript[dE, \[Gamma]]*)
(*dEnergiesGamma= Table[deltaGamma[energiesGamma[[i]]],{i,1,Length[energiesGamma]}];*)
dEnergiesGamma = Table[(energiesGamma[[i+1]]-energiesGamma[[i]])/2.,{i,1,Length[energiesGamma]-1}];

(*Function to calculate the energy integral given an array, "spectra" which must be values of dN/dE evaluated at values of energiesGamma. Returns scalar*)
enerInt[spectra_]:=Total[energiesGamma*spectra*dEnergiesGamma]

(*Function to generate bare plots of E^2dN/dE*)
(*specPlot[spec_]:=ListLogLogPlot[Thread[{energiesGamma*10^-9,energiesGamma^2 spec}],PlotRange->All]*)

specPlot[spec_]:=Module[{logSpec=Log10@(energies^2 spec),specInt},
specInt=Interpolation[Thread[{Log10@energies,logSpec}]];
Return[LogLogPlot[10.0^specInt[Log10[x]],{x,0.1,10^12}]];
];

(*Function that generates broken powerlaws for a given energy, spectral index, cutoff, and normalization*)
cutoffPowerLaw[enerGamma_,gamma_,cutoff_,amp_] := (amp*enerGamma^-gamma)*Exp[-enerGamma/cutoff]; 
(*Generates an array which holds the values of z which will be treated as the regions of constant EBL. These range from 0 to 3.9 in steps of 0.01*)
zReg=Range[0.,10.,0.01];

findIndex[val_]:=FirstPosition[zReg,SelectFirst[zReg,#>val&]][[1]]-1;

(*Import Differential spetral tables & inverse mean free paths*)
(*Declare containers for the diff spectra and inverse mean free path*)

(*Import inverse mean free path in Units of m^-1*)
IMFP=Import[FileNameJoin[{libraryLocation,"IMFP.csv"}]];

diffSpectra=Import[FileNameJoin[{libraryLocation,libraryNames<>".mat"}]];

(*
(*declare container for specctral table*)
diffSpectra=ConstantArray[0,Length[zReg]];

(*Initiate iterator*)
counter=1;

(*Loops of zVal and imports all the relevant tables storing them into the above containers*)
Do[
diffSpectra[[counter]]= Import[FileNameJoin[{libraryLocation,libraryNames<>"_z_"<>ToString[zVal]<>".csv"}]];
counter++;
,{zVal,zReg}]
*)
(*Generates a 3 dimensional array of extinction coefficients as a function of z and gamma energy using values from the above arrays and a distance of 1 kpc, \[Tau][zVal,energiesGamma]*)
(*extinctionCoeffs=Table[Exp[-kpc(**(1+zReg[[j]])*)*IMFP[[j,i]]],{j,1,Length[zReg]},{i,1,Length[energiesGamma]}];*)
extinctionCoeffs=Table[Exp[-kpc*IMFP[[j,i]]],{j,1,Length[zReg]},{i,1,Length[energiesGamma]}];


stepSizeArray=Import[FileNameJoin[{libraryLocation,"testzArray.mat"}]][[1]];
stepSizeArray=DeleteCases[stepSizeArray,0.,Infinity];

zRegIndexArray=Import[FileNameJoin[{libraryLocation,"zRegIndexArray.mat"}]][[1]];
zRegIndexArray=DeleteCases[zRegIndexArray,0.,Infinity];



$HistoryLength=0;
(*Setup diffuse background distance arrays*)
firstStep=Table[If[i==0,Table[10^-6*10^i,{j,1,10}],Table[10^-6*10^i,{j,1,9}]]//N,{i,0,4}]//Flatten ;
diffuseSteps={firstStep,ConstantArray[0.01,990]}//Flatten;
diffuseDistances=Table[Total[diffuseSteps[[;;i]]],{i,1,Length[diffuseSteps]}];
maxXIndex[x_]:=Position[diffuseDistances,Nearest[diffuseDistances,x][[1]]][[1,1]];
zMaxFunc[x_]:=diffuseDistances[[maxXIndex[x]]];



changeMagneticField[BField_Rational,name_String]:=Module[{
newIntergGalacticMagneticField=BField/10000,
powerProp={} ,(*(1.+zReg[[zRegionIndex]])**)(*(1.+zReg[[zRegionIndex]])**)
synchCorrectedSecGammaTable={},
NormalizedPairProductionSpectralTable,
NormalizedCascadedSecondaryGammaSpectralTable
},

Print["Importing pair production & ICS libraries"];
(*Import Pair Production Table*)

NormalizedPairProductionSpectralTable=ParallelTable[Import[libraryLocation<>"NormalizedPairProductionSpectralTable_z_"<>ToString[zVal]<>".csv"],{zVal,zReg},DistributedContexts -> {packageName<>"`Private`"}];

(*Import ICS Table*)
NormalizedCascadedSecondaryGammaSpectralTable=ParallelTable[Import[libraryLocation<>"NormalizedCascadedSecondaryGammaSpectralTable_z_"<>ToString[ToString[zVal]]<>".csv"],{zVal,zReg},DistributedContexts -> {packageName<>"`Private`"}];

(*series of arrays that parametrize this in oder to account for it as a function of redshift during propagation*)

powerProp=Table[powerProportionGen[ener*(1+zVal),t*(1+zVal),newIntergGalacticMagneticField*(1+zVal)^2],{zVal,zReg},{ener,energiesGamma}];

synchCorrectedSecGammaTable=powerProp*NormalizedCascadedSecondaryGammaSpectralTable;

Print["Exporting reweighed cascade libraries"];

Do[
Export[FileNameJoin[{libraryLocation,name<>"DiffCascadeSpectralTable_z_"<>ToString[zReg[[zValIndex]]]<>".csv"}],ParallelTable[Total[NormalizedPairProductionSpectralTable[[zValIndex,j]]*synchCorrectedSecGammaTable[[zValIndex]]*dEnergiesGamma],{j,1,Length[energiesGamma]},DistributedContexts -> {packageName<>"`Private`"}]];
,{zValIndex,1,Length[zReg]}];

Return["New magnetic field successfully implemented"];
];




CosmoCascadeV2[injSpectra_,zArrayLocal_,stepSizeArrayLocal_,zRegIndexArrayLocal_]:=Block[
{func,
 strechedResult,
 finalResultLocal,
 paramsLocal,
 singleCycle,
 stretchedEnergies = energiesGamma*((1. + zArrayLocal[[1]])/(1. + zArrayLocal[[-1]])),
 last=FirstPosition[Reverse@injSpectra,x_/;x!=0,1][[1]]
 },

 singleCycle[spectrum_,zRegionIndex_,thisStepSize_]:=Block[{attenuatedSpec,result,kern},         
        attenuatedSpec = (extinctionCoeffs[[zRegionIndex]]^((1.+zReg[[zRegionIndex]])*thisStepSize))*spectrum;
        kern=((spectrum - attenuatedSpec)*diffSpectra[[zRegionIndex]])[[2;;]]+((spectrum - attenuatedSpec)*diffSpectra[[zRegionIndex]])[[;;-2]];
        result = Total[kern*dEnergiesGamma] + attenuatedSpec;    
        (*result = Total[((spectrum - attenuatedSpec)*diffSpectra[[zRegionIndex]])[[;;-2]]*dEnergiesGamma] + attenuatedSpec;*)
		Return[result];
];

paramsLocal=Table[{zRegIndexArrayLocal[[i]],stepSizeArrayLocal[[i]]},{i,1,Length[stepSizeArrayLocal]}];
finalResultLocal = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, paramsLocal];
func = Interpolation[Thread[{energiesGamma[[;;-last]], finalResultLocal[[;;-last]]}]]; 
strechedResult = Table[If[stretchedEnergies[[i]]<=energiesGamma[[-last]], func[stretchedEnergies[[i]]],0.],{i,1,Length[stretchedEnergies]}];   
Return[strechedResult];
];


CosmoCascadeV3[injSpectra_,zArrayLocal_,stepSizeArrayLocal_,zRegIndexArrayLocal_]:=Block[
{func,
 strechedResult,
 finalResultLocal,
 paramsLocal,
 singleCycle,
 stretchedEnergies = energiesGamma*((1. + zArrayLocal[[1]])/(1. + zArrayLocal[[-1]])),
 last=FirstPosition[Reverse@injSpectra,x_/;x!=0,1][[1]]
 },

 singleCycle[spectrum_,zRegionIndex_,thisStepSize_]:=Block[{attenuatedSpec},         
        attenuatedSpec = (extinctionCoeffs[[zRegionIndex]]^((1.+zReg[[zRegionIndex]])*thisStepSize))*spectrum;
		Return[attenuatedSpec];
];

paramsLocal=Table[{zRegIndexArrayLocal[[i]],stepSizeArrayLocal[[i]]},{i,1,Length[stepSizeArrayLocal]}];
finalResultLocal = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, paramsLocal];
func = Interpolation[Thread[{energiesGamma[[;;-last]], finalResultLocal[[;;-last]]}]]; 
strechedResult = Table[If[stretchedEnergies[[i]]<=energiesGamma[[-last]], func[stretchedEnergies[[i]]],0.],{i,1,Length[stretchedEnergies]}];   
Return[strechedResult];
];


GCascadeAttenuate[injSpectraPre_,zStart_]:=Block[
{injSpectra=injSpectraPre*10^9,
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=maxXIndex[zStart],
 stepSize=1.*10^-6,
 luminosityLength
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energiesGamma],Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of dN/dE as a function of energy for your source"];Abort[];];

luminosityLength=(1+zStart)*NIntegrate[c/(H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2];


(*Calculates dV(z)*)
volumeNorms=Table[1.,{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

(*finalResult=Fold[CosmoCascadeV2[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]];*)
finalResult=((1+zStart)^2*Fold[CosmoCascadeV3[#1,#2[[1]],#2[[2]],#2[[3]]]&,injSpectra,Reverse[params]])/(4*\[Pi]*(luminosityLength*kpc*100)^2);

Return[finalResult*10^-9];
Clear["zFunc$*",
"luminosityLength$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


GCascadePoint[injSpectraPre_,zStart_]:=Block[
{injSpectra=injSpectraPre*10^9,
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=maxXIndex[zStart],
 stepSize=1.*10^-6,
 luminosityLength
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energiesGamma],Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of dN/dE as a function of energy for your source"];Abort[];];

luminosityLength=(1+zStart)*NIntegrate[c/(H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2];


(*Calculates dV(z)*)
volumeNorms=Table[1.,{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

(*finalResult=Fold[CosmoCascadeV2[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]];*)
finalResult=((1+zStart)^2*Fold[CosmoCascadeV2[#1,#2[[1]],#2[[2]],#2[[3]]]&,injSpectra,Reverse[params]])/(4*\[Pi]*(luminosityLength*kpc*100)^2);

Return[finalResult*10^-9];
Clear["zFunc$*",
"luminosityLength$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


GCascadeDiffuseConstant[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre*10^9,
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=maxXIndex[zStart],
 stepSize=1.*10^-6
 },
 
If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energiesGamma],Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of luminosity as a function of energy for your source"];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Successfully incorporated luminosity function. GCascade will produce spectra for diffuse background generated by your injected spectra."];,Print["Luminosity function is not formatted correctly. Please use a list of length "<>ToString[Length[zReg]]<>" containing values of luminosity density as a function of z from z=0 to z=10"];Abort[];];

(*Interpolates \[Rho](z) given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates dV(z)*)
volumeNorms=Table[((c*kpc*100*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/(4.*Pi*H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m*(1+diffuseDistances[[i]])^3]),{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]],zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[CosmoCascadeV2[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]];

Return[finalResult*10^-9];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*",
"injSpectra$*"
];
];


GCascadeDiffuse[injSpectraPre_,zStart_,zFunction_]:=Block[
{injSpectra=injSpectraPre*10^9,
 zArray,
 finalResult,
 zFunc,
 volumeNorms,
 params,
 zMaxIndex=maxXIndex[zStart],
 stepSize=1.*10^-6
 },

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]=={Length[zFunction],Length[energiesGamma]},Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of luminosity as a function of energy for your source"];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Successfully incorporated luminosity function. GCascade will produce spectra for diffuse background generated by your injected spectra."];,Print["Luminosity function is not formatted correctly. Please use a list of length "<>ToString[Length[zReg]]<>" containing values of luminosity density as a function of z from z=0 to z=10"];Abort[];];

(*Interpolates \[Rho](z) given by the user as zFunction*)
zFunc:=Interpolation[zFunction,InterpolationOrder->1];

(*Calculates dV(z)*)
volumeNorms=Table[((c*kpc*100*zFunc[diffuseDistances[[i]]])* diffuseSteps[[i]])/(4.*Pi*H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m*(1+diffuseDistances[[i]])^3]),{i,1,zMaxIndex}];

(*Calculates an array of list where every member is a list of z-values which are the step limits for the next function*)
zArray= Table[Range[diffuseDistances[[i+1]], diffuseDistances[[i]],-stepSize],{i,1,Length[volumeNorms]-1}];
zArray= Prepend[zArray,Range[diffuseDistances[[1]], 0.0,-stepSize]];

(*Creates a list of parameters to pass on to the next function*)
params = Table[{volumeNorms[[i]]*injSpectra[[i]],zArray[[i]],stepSizeArray[[i]],zRegIndexArray[[i]]},{i,1,Length[volumeNorms]}];

finalResult=Fold[CosmoCascadeV2[#1+#2[[1]],#2[[2]],#2[[3]],#2[[4]]]&,Table[0,{i,1,Length[energiesGamma]}],Reverse[params]];
Return[finalResult*10^-9];
Clear["zFunc$*",
"finalResult$*",
"volumeNorms$*",
"params$*",
"zArray$*",
"zMaxIndex$*"
"injSpectra$*"
];
];


End[ ]

EndPackage[ ]
