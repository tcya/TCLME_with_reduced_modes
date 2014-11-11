(* ::Package:: *)

BeginPackage["ModeAnalysis`"];

(* 

Copyright (C)2014  Eric R. Bittner& Xunmo  
This program is free software:you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation,either version 3 of the License,or
any later version. This program is distributed in the 
hope that it will be useful,but WITHOUT ANY WARRANTY;without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.You should have received a 
copy of the GNU General Public License
along with this program.If not,see http://www.gnu.org/licenses/.

*)

ReadRotM::usage = "Reads a QChem output file to find the diabatic RotMatrix";

ReadEne::usage = "Reads a QChem output file to find thee adiabaticH(0,0) and adiabaticH(1,1) ";

ReadNormalModes::usage = "Reads QCHem file to load normal modes";

ReadGradient::usage="Reads  QChem file to find gradient  ";

Project1by1::usage = "Project1by1 sorts phonon modes into a 1 by 1 hierarchical form. Only works for one mode case, for general 3 by 3 case it should be modified.  "



NMode::usage = "  "


 


ReadRotM[filename_]:=Module[{stream,RotMat},
stream=OpenRead[filename];
(*Exactly same as Eq .12 in Subotnik's paper. Subscript[H, D]=U^T Subscript[H, A]U*)
RotMat=Transpose[Partition[ToExpression[Last[Transpose[StringSplit[FindList[filename,"diabatic RotMatrix"]]]]],2]];
Close[stream];
Return[RotMat]]

ReadEne[filename_]:=Module[{stream,e1,e2},
stream=OpenRead[filename];
e1=ToExpression[Last[StringSplit[Find[stream,"showmatrix adiabatH(0,0)"]]]]*27.2;
e2=ToExpression[Last[StringSplit[Find[stream,"showmatrix adiabatH(1,1)"]]]]*27.2(*in eV*);
Close[stream];
Return[{e1,e2}](* in hartree *)
]

ReadNormalModes[filename_](* for QChem*):= Module[{stream, NAtom, NMode, NBlock, Freq, ForceConst, RedMass, IRActive, IRIntens, RamanActive, NormalModeCoord, IAn, AtomicMass, tmp},
  Freq = {};
  ForceConst = {};
  RedMass = {};
  IRActive = {};
  IRIntens = {};
  RamanActive = {};
  NormalModeCoord = {};
  IAn = {};
  AtomicMass = {};
  
  stream = OpenRead[filename];
  NMode = ToExpression[Last[Last[StringSplit[FindList[stream, "Mode"]]]]];
  NBlock = Floor[(NMode - 0.1)/3] + 1;
  NAtom = (NMode + 6)/3;
  Close[stream];
  
  stream = OpenRead[filename];
  Find[stream, "Raman Active"];
  Read[stream];
  tmp = Table[First[StringSplit[Read[stream, String]]], {NAtom}];
  AtomicMass = ElementData[#, "AtomicWeight"] & /@ tmp;
  IAn = ElementData[#, "AtomicNumber"] & /@ tmp;
  Close[stream];
  
  stream = OpenRead[filename];
  Do[
   Find[stream, "Mode"];
   AppendTo[Freq, ToExpression[Drop[StringSplit[Read[stream, String]], 1]]];
   AppendTo[ForceConst, ToExpression[Drop[StringSplit[Read[stream, String]], 2]]];
   AppendTo[RedMass, ToExpression[Drop[StringSplit[Read[stream, String]], 2]]];
   AppendTo[IRActive, Drop[StringSplit[Read[stream, String]], 2]];
   AppendTo[IRIntens, ToExpression[Drop[StringSplit[Read[stream, String]], 2]]];
   AppendTo[RamanActive, Drop[StringSplit[Read[stream, String]], 2]];
   Read[stream];
   NormalModeCoord = Join[NormalModeCoord, ToExpression[Transpose[Table[Partition[Rest[StringSplit[Read[stream, String]]], 3], {NAtom}]]]];
   , {i, 1, NBlock}];
  
  Close[stream];
  {Freq, ForceConst, RedMass, IRActive, IRIntens, RamanActive} = Flatten /@ {Freq, ForceConst, RedMass, IRActive, IRIntens, RamanActive};
  Return[{NAtom, NMode, AtomicMass, Freq, RedMass, NormalModeCoord}]
]

ReadGradient[filename_, NAtom_] := Module[{stream, grad, NBlock, tmp},
  grad = {};
  NBlock = Floor[NAtom/6] + 1;
  
  stream = OpenRead[filename];
  Do[Find[stream, "Gradient of the state energy (including CIS Excitation Energy)"];,{1}];
  Find[stream, "Gradient of the state energy (including CIS Excitation Energy)"];
   Do[
    Read[stream];
    grad = Join[grad, Transpose[Table[ToExpression[Rest[StringSplit[Read[stream, String]]]], {3}]]];
    , {i, 1, NBlock}];
  Return[grad]
]



(* Project1by1 sorts phono modes into a 1 by 1 hierarchical form. Only works for one mode case, for general 3 by 3 case it should be modified. *)
Project1by1[gg_,Hess_]:=Module[
{S,P,Q,PHP,QHQ,L,\[Omega]c,\[Omega]ck,\[Omega]b,Mc,Mck,Mb,np,ns,Sinv,G,p,k,g},
(* g should be 1-D, namely, ns=1 *)
{ns,np}=Dimensions[gg];
g=gg;
\[Omega]c={};
Mc={};
p=ConstantArray[0,{np,np}];

Do[{
S = Table[g[[i]].g[[j]],{i,ns},{j,ns}];
Sinv = Inverse[S] ;
(* define projection operators *)

P = Sum[Sinv[[i,j]]  g[[i]]\[CircleTimes]g[[j]],{i,ns},{j,ns}] ; 
p=p+P; (* P=Subscript[P, 1]+Subscript[P, 2]+... *)
Q = IdentityMatrix[np]-p;
PHP = P.Hess.P ;
QHQ = Q.Hess.Q ; 
{\[Omega]ck,Mck}=Chop[Transpose[Take[Sort[Transpose[Eigensystem[PHP]]],-ns]],10^-8];
{\[Omega]b,Mb}=Chop[Transpose[Take[Sort[Transpose[Eigensystem[QHQ]]],ns+k-1-np]],10^-8];
\[Omega]c=Chop[Join[\[Omega]c,\[Omega]ck],10^-8];
Mc=Chop[Join[Mc,Mck],10^-8];
\[Omega]all= Join[\[Omega]c,\[Omega]b]//Chop;
M = Join[Mc,Mb]//Chop ;
g=Chop[{Join[ConstantArray[0,k],( M.Hess.Transpose[M])[[k]][[k+1;;np]]].M},10^-8]; (* Define nth mode=(0,0,0,...,Subscript[\[Gamma], k+1],Subscript[\[Gamma], k+2],...,Subscript[\[Gamma], N]) *)
Gnew=M.gg[[1]];
},{k,1,np-1}];
](* The kth step gives k*k tridiagonal submatrix on lefttop, diagonal N-k on rightbottom and only the kth mode is coupled to the N-k bath modes *)


Begin["Private`"]

A_\[CircleTimes]B_:= Outer[Times,A,B]


NMode[n_,hess_,M_,g_]:=Module[{h,gn,Mnm,\[Omega]prime,Gn},
h=Take[(M.hess.Transpose[M]),{1,n},{1,n}];
{\[Omega]prime,Mn}=Transpose[Sort[Transpose[Eigensystem[h]]]];
gn=Join[Take[g,n],ConstantArray[0,n-Length[Take[g,n]]]];
Gn=Mn.gn;
Return[{Sqrt[\[Omega]prime],Gn}]
]



End[]
EndPackage[]

