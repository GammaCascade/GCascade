(* ::Package:: *)

BeginPackage["GCascadeV1`"]



GAttenuate::usage = "CosmoAttenuate[injSpectra,zStart,stepSize]: CAttenuate is a function which takes in a gamma-ray injected spectrum, starting z, and step size in z. This function produces an attenuated differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)] \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\) , without including cascade  evolution but taking into account cosmological energy shifts. The injected spectrum should be formatted as an array of values of differential flux \!\(\*FormBox[\((\),
TraditionalForm]\)\!\(\*FormBox[\(\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\)\),
TraditionalForm]\)) evaluated at the energies given by the array, energiesGamma. It is recommended to use a step size of up to \!\(\*FormBox[\((\*SuperscriptBox[\(10\), \(-5\)]\),
TraditionalForm]\)).";



energiesGamma::usage="a logarithmic array of 300 energies from 10^8 to 10^21 eV, these serve as the energies of gamma-rays";
specPlot::usage="Plots a primitive frame of \!\(\*SuperscriptBox[\(E\), \(2\)]\)\!\(\*FractionBox[\(dN\), \(dE\)]\) given some \!\(\*FractionBox[\(dN\), \(dE\)]\) as a single array evaluated at the values of energiesGamma";
cutoffPowerLaw::usage = "cutoffPowerLaw[energy (E),index (\[CapitalGamma]),cutoff (\!\(\*SubscriptBox[\(E\), \(cutoff\)]\)),normalization (N)]: Function that generates spectrum values of broken powerlaws (\!\(\*FractionBox[\(dN\), \(dE\)]\)= N*\!\(\*SuperscriptBox[\(E\), \(-\[CapitalGamma]\)]\)*Exp[-\!\(\*FractionBox[\(E\), SubscriptBox[\(E\), \(cutoss\)]]\)]) for a given energy, spectral index, cutoff, and normalization";
zReg::usage = "Array which holds the values of z which will be treated as the regions of constant EBL. These range from 0 to 10 in steps of 0.01";
diffuseDistances::usage="Array which holds values of Z which will be treated as the points of cascade evaluation";
diffuseSteps::usage="Array which holds values of \[CapitalDelta]Z which will be treated as the limits of cascade evaluation";


GCascadeDiffuse::usage = "CCascadeDiffuse[injected spectrum ,max Z (\!\(\*SubscriptBox[\(Z\), \(max\)]\)),step size (dZ),Luminosity Function (L(z))]:CCascadeDiffuse is a function which takes in a gamma-ray injected luminosity spectrum, max z, step size in z, and a luminosity function. 'GCascadeDiffuse' produces a cascaded diffuse differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)] \*SuperscriptBox[\(cm\), \(-2\)]\\\ \*SuperscriptBox[\(s\), \(-1\)]\\\ \*SuperscriptBox[\(sr\), \(-1\)])\))\),
TraditionalForm]\) , taking into account cosmological and electromagnetic energy losses. It is recommended to use a step size of \!\(\*FormBox[SuperscriptBox[\(10\), \(-6\)],
TraditionalForm]\). The injected differential spectrum \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)])\),
TraditionalForm]\) should be formatted as an array of values of power per unit area \!\(\*FormBox[\((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\),
TraditionalForm]\) \!\(\*FormBox[SuperscriptBox[\(s\), \(-1\)],
TraditionalForm]\)) evaluated at the energies given by the array, 'energiesGamma'.  The Luminosity distribution function is assumed to be a separable function of the form: \!\(\*FormBox[\(\*FractionBox[\(dN\), \(\*SubscriptBox[\(dV\), \(c\)] dE\)] \((\*SubscriptBox[\(E\), \(\[Gamma]\)], z)\) = \*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SubscriptBox[\(E\), \(\[Gamma]\)])\) \(\[Rho](z)\)\),
TraditionalForm]\) where the left-hand factor is the intrinsic differential spectrum of the class of sources \!\(\*FormBox[\((\*FractionBox[\(dN\), SubscriptBox[\(dE\), \(\[Gamma]\)]] \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\) and \[Rho](z) is the spacial distribution of the sources in units of density per comoving volume, so called 'luminosity function'. 'GCascadeDiffuse' uses the luminosity function as an input which must be formatted as list of ordered pairs (z,\!\(\*FormBox[\(\[Rho](z)\),
TraditionalForm]\)), of values of z and comoving density \[Rho](z). The values of z can be any finite list of z from 0 to 10.";

GCascadePoint::usage = "CCascadePoint[injected spectrum, source distance (Z),step size (dZ)]: CCascadePoint is a function which takes in a gamma-ray injected spectrum, source distance, and step size in Z. 'GCascadePoint' produces a cascaded differential spectrum, \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)]\\\ \((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(s\), \(-1\)])\))\),
TraditionalForm]\), taking into account cosmological and electromagnetic energy losses. It is recommended to use a step size of \!\(\*FormBox[SuperscriptBox[\(10\), \(-6\)],
TraditionalForm]\). The injected differential spectrum \!\(\*FormBox[\((\*FractionBox[\(dN\), \(dE\)])\),
TraditionalForm]\) should be formatted as an array of values of power per unit area \!\(\*FormBox[\((\*SuperscriptBox[\(GeV\), \(-1\)]\\\ \*SuperscriptBox[\(cm\), \(-2\)]\),
TraditionalForm]\) \!\(\*FormBox[SuperscriptBox[\(s\), \(-1\)],
TraditionalForm]\)) evaluated at the energies given by the array, 'energiesGamma'.";

changeMagneticField::usage = "changeMagneticField[BField_Rational,name_String] is a function which changes the IG magnetic field of CCascade by exporting new libraries. See tutorial for complete instructions.";



Begin["`Private`"]
Print["You Have Imported GCascade. Setup will take about 2 minutes."]; 

libraryLocation = NotebookDirectory[]<>"Libraries\\"; 

libraryNames="DiffCascadeSpectralTable";

kpc = 3.086*10^19; h = 4.135667/10.0^15; joule2eV = 6.242*10.0^18.0; c = 299792458.0; 
k = 8.6173324/10^5; 
t = 2.725; electronMass = 510998.; 
classicalRadius = 2.81794/10^15; 
thomp = (8./3)*Pi*classicalRadius^2; 
\[Mu] = 1.2566/10^6; interGalacticMagneticField = 1./10^13/10000; 


(*Cosmological parameters from Planck '15*)

packageName="GCascadeV1";
H0=67.7 (*(km/s)/Mpc*);
\[CapitalOmega]\[CapitalLambda]=0.691;
\[CapitalOmega]m=0.308;
 
(*Function which calculates the distance light travels from zCurrent to some redshift zCurrent-dz*)
lightDistanceTravelled[zCurrent_,dz_]:=NIntegrate[c/(H0(1+zVal)Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,zCurrent-dz,zCurrent},Method->{"Trapezoidal"},PrecisionGoal->2];

(*Black Body differential Number Density*)

dndEtar[Etar_,t_] := ((8.0*Pi)/(h*c)^3)*(Etar^2/(Exp[Etar/(k*t)]-1.))



(*Generates an array of n logarithmically spaced values from 10^a to 10^b*)
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)]

(*Generates a logarithmic array of 300 energies from 10^8 to 10^20 eV, these serve as the energies of gamma-rays *)
energiesGamma = logspace[8,21,300];

(*Function to calculate the spacing between energies in the gamma-ray energy array*)
deltaGamma[enerGamma_] := enerGamma*((10.0^((21.-8.)/(300.-1.))) - 1);
(*Generates an array which holds the logarithmic spacing in the gamma-ray energy array, i.e. an array of values Subscript[dE, \[Gamma]]*)
dEnergiesGamma= Table[deltaGamma[energiesGamma[[i]]],{i,1,Length[energiesGamma]}];

(*Function to calculate the energy integral given an array, "spectra" which must be values of dN/dE evaluated at values of energiesGamma. Returns scalar*)
enerInt[spectra_]:=Total[energiesGamma*spectra*dEnergiesGamma]

(*Function to generate bare plots of E^2dN/dE*)
(*specPlot[spec_]:=ListLogLogPlot[Thread[{energiesGamma*10^-9,energiesGamma^2 spec}],PlotRange->All]*)

specPlot[spec_]:=LogLogPlot[Interpolation[Thread[{energiesGamma*10^-9,(energiesGamma^2 spec)*10^-9}]][x],{x,0.1,10^12}];

(*Function that generates broken powerlaws for a given energy, spectral index, cutoff, and normalization*)
cutoffPowerLaw[enerGamma_,gamma_,cutoff_,amp_] := (amp*enerGamma^-gamma)*Exp[-enerGamma/cutoff]; 
(*Generates an array which holds the values of z which will be treated as the regions of constant EBL. These range from 0 to 3.9 in steps of 0.01*)
zReg=Range[0.,10.,0.01];

(*Import Differential spetral tables & inverse mean free paths*)
(*Declare containers for the diff spectra and inverse mean free path*)

(*Import inverse mean free path*)
IMFP=Import[libraryLocation<>"IMFP.csv"];

(*declare container for specctral table*)
diffSpectra=ConstantArray[0,Length[zReg]];

(*Initiate iterator*)
counter=1;

(*Loops of zVal and imports all the relevant tables storing them into the above containers*)
Do[
diffSpectra[[counter]]= Import[libraryLocation<>libraryNames<>"_z_"<>ToString[zVal]<>".csv"];
counter++;
,{zVal,zReg}]

(*Generates a 3 dimensional array of extinction coefficients as a function of z and gamma energy using values from the above arrays and a distance of 1 kpc, \[Tau][zVal,energiesGamma]*)
extinctionCoeffs=Table[Exp[-kpc*IMFP[[j,i]]],{j,1,Length[zReg]},{i,1,Length[energiesGamma]}];


GAttenuate[injSpectra_,zStart_,stepSize_]:=Quiet[Module[
{
zCur=Range[zStart,0,-stepSize],
luminosityLength=(1+zStart)*NIntegrate[c/(H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2],
zCurReg={},
stepArray={},
finalResult={},
singleCycle
},

zCurReg=Table[Nearest[zReg,zCur[[i]],1][[1]],{i,1,Length[zCur]}];
stepArray= Table[lightDistanceTravelled[zCurReg[[i]],stepSize],{i,1,Length[zCur]}];

singleCycle[spectrum_,zCurIndex_]:=Module[{
func={},
zRegionIndex=Position[zReg,zCurReg[[zCurIndex]]][[1,1]],
attenuatedSpec={},
result={},
strechedResult={},
stretchedEnergies=energiesGamma*(1+zCurReg[[zCurIndex]])/(1+zCurReg[[zCurIndex]]+stepSize)
},
attenuatedSpec=extinctionCoeffs[[zRegionIndex]]^stepArray[[zCurIndex]]*spectrum;
func=Interpolation[Thread[{energiesGamma,attenuatedSpec}],InterpolationOrder->1];
strechedResult=Table[func[i],{i,stretchedEnergies}];
Return[strechedResult];
];
finalResult=Last[FoldList[singleCycle,injSpectra,Range[1,Length[zCur]-1]]]/(4*\[Pi]*(luminosityLength*kpc*100)^2) (*GeV^-1 cm^-2 s^-1*);
 Return[finalResult];
]];

(*Develop functions to calculate Klein-Nishina step ICS supression*)

(*CMB energy Density*)
CMBdensity[temp_]:=Integrate[dndEtar[enerTar,temp]*enerTar,{enerTar,0,Infinity}];

(*CMB Gamma1*)
gamma1[temp_]:= Evaluate[((3.*Sqrt[5.]*electronMass)/(8*Pi*k*temp))];

(*Power going into ICS*)
comptonPower[enerElec_,temp_]:=Evaluate[(4./3.)thomp*c*(enerElec^2/electronMass^2)*(CMBdensity[temp]/(1.+(enerElec^2/(electronMass^2*gamma1[temp]^2))))];

(*Magnetic Field Density*)
magneticEnergyDensity[magneticField_]:=Evaluate[(magneticField^2/(2.*\[Mu]))*((6.242*(10^18)(*eV*))/1.(*Joule*)) ];

(*Power going into Synchrotron*)
synchPower[enerElec_,magneticField_]:= Evaluate[(4./3.)thomp*c*(enerElec^2/electronMass^2)*magneticEnergyDensity[magneticField]];

(*Proportion of power going into ICS*)
powerProportionGen[enerElec_,temp_,magneticField_] := comptonPower[enerElec,temp]/(comptonPower[enerElec,temp]+synchPower[enerElec,magneticField])  ;

(*series of arrays that parametrize this in oder to account for it as a function of redshift during propagation*)

powerProp=ParallelTable[powerProportionGen[energiesGamma[[i]]*(1+zVal),t(1+zVal),interGalacticMagneticField*(1+zVal)^2],{zVal,1,10,0.01},{i,1,Length[energiesGamma]}];



$HistoryLength=0;
(*Setup diffuse background distance arrays*)
firstStep=Table[If[i==0,Table[10^-6*10^i,{j,1,10}],Table[10^-6*10^i,{j,1,9}]]//N,{i,0,4}]//Flatten ;
diffuseSteps={firstStep,ConstantArray[0.1,99]}//Flatten;
diffuseDistances=Table[Total[diffuseSteps[[;;i]]],{i,1,Length[diffuseSteps]}];
maxXIndex[x_]:=Position[diffuseDistances,Nearest[diffuseDistances,x][[1]]][[1,1]];
zMaxFunc[x_]:=diffuseDistances[[maxXIndex[x]]];



CosmoCascade[injSpectra_,zStart_,stepSize_]:=Module[
{zCur = Range[zStart, 0, -stepSize],
 zCurReg,
 stepArray,
 finalResult,
 params,
 singleCycle,
 stepArrayZReg},

singleCycle[spectrum_,thiszCurReg_, thisStepSize_]:=
Module[
{func,
 zRegionIndex = Position[zReg,thiszCurReg][[1,1]], 
 attenuatedSpec,
 result,
 strechedResult,
 stretchedEnergies = energiesGamma*((1. + thiszCurReg)/(1. + thiszCurReg + stepSize))},
 attenuatedSpec = extinctionCoeffs[[zRegionIndex]]^thisStepSize*spectrum; 
 result = Total[(spectrum - attenuatedSpec)*diffSpectra[[zRegionIndex]]*dEnergiesGamma] + attenuatedSpec; 
 func = Interpolation[Thread[{energiesGamma, result}], InterpolationOrder -> 1]; 
 strechedResult = Table[func[i],{i, stretchedEnergies}]; 
 Return[strechedResult]
];
stepArrayZReg=Table[lightDistanceTravelled[zReg[[i]], stepSize],{i,Length[zReg]}];
zCurReg = Table[Nearest[zReg, zCur[[i]], 1][[1]], {i, 1, Length[zCur]}];
stepArray=Table[stepArrayZReg[[Position[zReg,zCurReg[[i]]][[1,1]]]],{i,Length[zCurReg]}];
params=Table[{zCurReg[[i]],stepArray[[i]]},{i,1,Length[zCur] - 1}];
finalResult = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, params];
Return[Evaluate[finalResult]]
];

GCascadePoint[injSpectra_,zStart_,stepSize_]:=Module[
{
finalResult={},
luminosityLength=(1+zStart)*NIntegrate[c/(H0*Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m((1+zVal)^3)]),{zVal,0,zStart},Method->{"Trapezoidal"},PrecisionGoal->2]
},

(*Check formating of Injected spectrum*)
If[ToString[Head[injSpectra]]==ToString[List] && Dimensions[injSpectra]== Dimensions[energiesGamma],
Print["Injected spectrum properly formatted"],
Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of luminosity as a function of energy for your source"];
Abort[];
 ];


finalResult= CosmoCascade[injSpectra,zStart,stepSize]/(4*\[Pi]*(luminosityLength*kpc*100)^2);

 Return[finalResult];

];


changeMagneticField[BField_Rational,name_String]:=Module[{
newIntergGalacticMagneticField=BField/10000,
powerProp={} ,
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
Export[libraryLocation<>name<>"DiffCascadeSpectralTable_z_"<>ToString[zReg[[zValIndex]]]<>".csv",ParallelTable[Total[NormalizedPairProductionSpectralTable[[zValIndex,j]]*synchCorrectedSecGammaTable[[zValIndex]]*dEnergiesGamma],{j,1,Length[energiesGamma]},DistributedContexts -> {packageName<>"`Private`"}]];
,{zValIndex,1,Length[zReg]}];

Return["New magnetic field successfully implemented"];
];




CosmoCascadeV2[injSpectra_,zStart_,zEnd_,stepSize_]:=Module[
{zCur = Range[zStart, zEnd, -stepSize],
 zCurReg,
           stepArray,
 finalResult,
 params,
singleCycle,
stepArrayZReg},

singleCycle[spectrum_,thiszCurReg_, thisStepSize_]:=
Module[
{func,
 zRegionIndex = Position[zReg,thiszCurReg][[1,1]], 
              attenuatedSpec,
 result,
 strechedResult,
 stretchedEnergies = energiesGamma*((1. + thiszCurReg)/
            (1. + thiszCurReg + stepSize))},
 
        attenuatedSpec = extinctionCoeffs[[zRegionIndex]]^thisStepSize*spectrum; 
         result = Total[(spectrum - attenuatedSpec)*diffSpectra[[zRegionIndex]]*dEnergiesGamma] + attenuatedSpec; 
         func = Interpolation[Thread[{energiesGamma, result}], InterpolationOrder -> 1]; 
         strechedResult = Table[func[i],{i, stretchedEnergies}]; 
Return[strechedResult]
];

stepArrayZReg=Table[lightDistanceTravelled[zReg[[i]], stepSize],{i,Length[zReg]}];
zCurReg = Table[Nearest[zReg, zCur[[i]], 1][[1]], {i, 1, Length[zCur]}];
stepArray=Table[stepArrayZReg[[Position[zReg,zCurReg[[i]]][[1,1]]]],{i,Length[zCurReg]}];
params=Table[{zCurReg[[i]],stepArray[[i]]},{i,1,Length[zCur] - 1}];
      finalResult = Fold[singleCycle[#1,#2[[1]],#2[[2]]]&, injSpectra, params];

      Return[Evaluate[finalResult]]
];


GCascadeDiffuse[injSpectra_,zStart_,stepSize_,zFunction_]:=Module[
{finalResult={},zFunc={},volumeNorms={},params={}},

If[ToString[Head[injSpectra]]==ToString[List]&&Dimensions[injSpectra]==Dimensions[energiesGamma],Print["Injected spectrum properly formatted"],Print["Injected spectrum is not formatted correctly. Please use a list of length "<>ToString[Length[energiesGamma]]<>" containing values of luminosity as a function of energy for your source"];Abort[];];

If[ToString[Head[zFunction]]==ToString[List]&&Dimensions[zFunction]==Dimensions[Thread[{diffuseDistances,diffuseDistances}]],Print["Successfully incorporated luminosity function. CCascade will produce spectra for diffuse background generated by your injected spectra."];,Print["Luminosity function is not formatted correctly. Please use a list of length "<>ToString[Length[zReg]]<>" containing values of luminosity density as a function of z from z=0 to z=10"];Abort[];];

zFunc:=Interpolation[zFunction,InterpolationOrder->1];

finalResult=FoldList[CosmoCascadeV2[]];

volumeNorms=Table[((c*kpc* 100 *zFunc[diffuseDistances[[i]]]) diffuseSteps[[maxXIndex[diffuseDistances[[i]]]]])/(H0 Sqrt[\[CapitalOmega]\[CapitalLambda]+\[CapitalOmega]m (1+diffuseDistances[[i]])^3] (1+diffuseDistances[[i]])^2),{i,1,Length[diffuseDistances[[1;;maxXIndex[zStart]]]]}] ;
params = Table[{volumeNorms[[i]],diffuseDistances[[i]],diffuseDistances[[i-1]]},{i,2,Length[volumeNorms]}];
params = Prepend[params,{zFunc[diffuseDistances[[1]]],diffuseDistances[[1]],0}];
finalResult=Fold[CosmoCascadeV2[#1+(#2[[1]]*injSpectra),#2[[2]],#2[[3]],stepSize]&,Table[0,{i,1,Length[injSpectra]}],Reverse[params]];
Return[finalResult];
Clear["params$*",
"zFunc$*",
"finalResult$*",
"volumeNorms$*",
"zCur$*",
"zCurReg$*",
"stepArray$*",
 "params$*",
"singleCycle$*",
"stepArrayZReg$*"
];
];


End[ ]

EndPackage[ ]
