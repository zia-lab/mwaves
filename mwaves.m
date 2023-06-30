(* -----------------------------------------------------------------

+------------------------------------------------------------------+
|                                                                  |
|                                                                  |
|             ____ ___     _      ______ __   _____  _____         |
|            / __ `__ \   | | /| / / __ `/ | / / _ \/ ___/         |
|           / / / / / /   | |/ |/ / /_/ /| |/ /  __(__  )          |
|          /_/ /_/ /_/    |__/|__/\__,_/ |___/\___/____/           |
|                                                                  |
|                                                                  |
|                                                                  |
+------------------------------------------------------------------+

This code was authored by David Lizarazo in 2023.  It  implements  a
few useful functions for the numerical analysis of waves propagating
in homogeneous media.  The code is written in the Wolfram Language.

----------------------------------------------------------------- *)

BeginPackage["mwaves`"]

moduleDir = DirectoryName[$InputFileName];

(*Load colormap data*)

colorMapFile    = FileNameJoin[{moduleDir, "data", "colormaps.txt"}];
colorMapRawData = ReadList[colorMapFile, Record, RecordSeparators -> {"\n\n", "\n"}];
colormapNames   = colorMapRawData[[1 ;; ;; 257]];
colormapData    = #[[2 ;;]] & /@ Partition[colorMapRawData, 257];
colormapData    = ToExpression /@ colormapData;
colormaps       = AssociationThread[colormapNames, colormapData];

GetColorMap::usage = "GetColorMap[name] returns a color function for the colormap with the given name. Data for these colormaps was ported from the Python package cmasher.
This function returns a color function that can be used with the ColorFunction option of Plot, DensityPlot, etc.

The available colormaps are: amber, amber_r, amethyst, amethyst_r, apple, apple_r, arctic, arctic_r, bubblegum, bubblegum_r, chroma, chroma_r, copper, copper_r, copper_s, copper_s_r, cosmic, cosmic_r, dusk, dusk_r, eclipse, eclipse_r, ember, ember_r, emerald, emerald_r, emergency, emergency_r, emergency_s, emergency_s_r, fall, fall_r, flamingo, flamingo_r, freeze, freeze_r, fusion, fusion_r, gem, gem_r, ghostlight, ghostlight_r, gothic, gothic_r, guppy, guppy_r, holly, holly_r, horizon, horizon_r, iceburn, iceburn_r, infinity, infinity_r, infinity_s, infinity_s_r, jungle, jungle_r, lavender, lavender_r, lilac, lilac_r, neon, neon_r, neutral, neutral_r, nuclear, nuclear_r, ocean, ocean_r, pepper, pepper_r, pride, pride_r, prinsenvlag, prinsenvlag_r, rainforest, rainforest_r, redshift, redshift_r, sapphire, sapphire_r, savanna, savanna_r, seasons, seasons_r, seasons_s, seasons_s_r, seaweed, seaweed_r, sepia, sepia_r, sunburst, sunburst_r, swamp, swamp_r, torch, torch_r, toxic, toxic_r, tree, tree_r, tropical, tropical_r, viola, viola_r, voltage, voltage_r, waterlily, waterlily_r, watermelon, watermelon_r, wildfire, wildfire_r, heat, heat_r.
"

GetColorMap[name_String] := 
  Module[{rgbList}, rgbList = colormaps[name];
   With[{colorData = rgbList}, 
    RGBColor[colorData[[Round[1 + #*(Length[rgbList] - 1)]]]] &]];


ScalarFieldProp::usage = "ScalarFieldProp takes a field component in an aperture plane and propagates that to an observation plane by using an implementation of the direct integration of the Rayleigh-Sommerfeld diffraction integral. This implementation is based on the method described in Shen and Wang (2006). The field is sampled in the aperture plane using a uniform grid, and the observation plane is also sampled using a uniform grid. The field is assumed to be zero outside of the aperture plane. Contrary to the Fresnel or Kirchhoff integrals, the Rayleigh-Sommerfeld integral does not require the use of the paraxial approximation, and is therefore valid for arbitrary distances between the aperture plane and the observation plane. Here these two planes are assumed to be parallel.

Parameters
----------
+ Lobs (float): spatial width of the obsevation region, in μm. The observation region is assumed to be a squared centered on (x,y) = (0,0), and extending from -Lobs/2 to Lobs/2 in both the x and y directions.
+ z (float): distance between the aperture plane and the observation plane, given in μm. The aperture plane is assumed to be at z=0.
+ apertureFunction (function): a bi-variate function that returns the complex amplitude of the field in the aperture plane. Input to the function is assumed to be in cartesian coordinates x,y.
+ λfree (float): wavelength in vaccum of field, given in μm.
+ nref  (float): refractive index of the propagating medium.

Options
-------
+ \"numSamples\" (int or Automatic): number of samples to use in the aperture plane and the observation plane. The aperture plane is sampled using a uniform grid, and the observation plane is also sampled using a uniform grid. The default is Automatic in which case numSamples is calculated so that the sample size is equal to half the wavelength of the wave inside of the propagating medium.

Returns
-------
{numSamples, xCoords, yCoords, field}
+ xCoords (List): x coordinates of the observation plane, given in μm.
+ yCoords (List): y coordinates of the observation plane, given in μm.
+ field   (List): complex amplitude of the field in the observation plane. The top left corner of the array corresponds to the lower left corner of the observation plane. The coordinates associated with each element in the given list should be taken from xCoords and yCoords.

References
----------
+ Shen, Fabin, and Anbo Wang. \"Fast-Fourier-transform based numerical integration method for the Rayleigh-Sommerfeld diffraction formula.\" Applied optics 45, no. 6 (2006): 1102-1110.

Example (double slit diffraction)
---------------------------------
(*geometry of the double slits*)
{slitSep, slitWidth, slitLen} = {3, 1, 4};

(*full width of the observation plane*)
Lobs = 10;
numSamples = 100;
z    = 10;
\[Lambda] = 0.532;

doubleSlit = 
  Compile[{{x, _Real}, {y, _Real}}, (If[
      And[Abs[x - slitSep/2] < slitWidth/2, Abs[y] < slitLen/2], 1., 0.] 
      + If[And[Abs[x + slitSep/2] < slitWidth/2, Abs[y] < slitLen/2],
      1., 
      0.])];

({xCoord, yCoord, field} = ScalarFieldProp[Lobs, z, doubleSlit, \[Lambda]];
 ArrayPlot[Abs[field], FrameTicks -> Automatic, DataReversed -> True, 
  ImageSize -> 800, 
  DataRange -> {{-Lobs/2, Lobs/2}, {-Lobs/2, Lobs/2}}])

"
Options[ScalarFieldProp] = {"numSamples" -> Automatic};
ScalarFieldProp[Lobs_, z_, apertureFunction_, \[Lambda]free_, nref_:1, OptionsPattern[]] := (
  \[Lambda] = \[Lambda]free / nref;
  If[OptionValue["numSamples"] == Automatic,
   numSamples = Round[2*Lobs/\[Lambda]];
   ,
   numSamples = OptionValue["numSamples"];
   ];
  (*1D samples in the observation plane, must be equal*)
  apertureSamples = numSamples;
  (*width of the aperture region, in um, 
  needs to be the same as the width of the observation region*)
  Lap = Lobs;
  (*wave number*)
  k = 2. \[Pi]/\[Lambda];
  (*sample size in aperture*)
  \[CapitalDelta]\[Zeta] = Lap / apertureSamples;
  \[CapitalDelta]\[Eta] = \[CapitalDelta]\[Zeta];

  g = Compile[
    {{x, _Real}, {y, _Real}},
    (E^(I k Sqrt[x^2 + y^2 + z^2])
      z (-I k + 1/Sqrt[x^2 + y^2 + z^2]))/(
    2. \[Pi] (x^2 + y^2 + z^2))];
  gr = Compile[
    {{r, _Real}},
    z/(2.*\[Pi]) (E^(I k  r) *(-I k + 1./r))/r^2];
  
  (*\[Zeta] and \[Eta] are coordinates in the source plane*)
  \[Zeta]    = Range[-Lap/2, Lap/2, Lap/(apertureSamples - 1)] // N;
  \[Eta]     = Range[-Lap/2, Lap/2, Lap/(apertureSamples - 1)] // N;

  (*x and y are coordinates in the observation plane*)
  x          = Range[-Lobs/2, Lobs/2, Lobs/(numSamples - 1)] // N;
  y          = Range[-Lobs/2, Lobs/2, Lobs/(numSamples - 1)] // N;
  coordArray = Outer[List, \[Zeta], \[Eta]];
  
  (*Put Together the U matrix*)
  U = Map[apertureFunction @@ # &, coordArray, {2}];
  U = PadRight[U, {2*apertureSamples - 1, 2*apertureSamples - 1}];
  
  Hx1       = ConstantArray[x[[1]], {apertureSamples - 1, 2 apertureSamples - 1}];
  Hx2       = Transpose[ConstantArray[x, 2 apertureSamples - 1]];
  Hx        = ArrayFlatten[{{Hx1}, {Hx2}}];
  H\[Zeta]2 = ConstantArray[\[Zeta][[1]], {apertureSamples - 1, 2 apertureSamples - 1}];
  H\[Zeta]1 = Transpose[ConstantArray[Reverse[\[Zeta]], 2 apertureSamples - 1]];
  H\[Zeta]  = ArrayFlatten[{{H\[Zeta]1}, {H\[Zeta]2}}];
  Hx\[Zeta] = Hx - H\[Zeta];
  
  Hy1      = ConstantArray[y[[1]], {apertureSamples - 1, 2 apertureSamples - 1}];
  Hy2      = Transpose[ConstantArray[y, 2 apertureSamples - 1]];
  Hy       = Transpose[ArrayFlatten[{{Hy1}, {Hy2}}]];
  H\[Eta]2 = ConstantArray[\[Eta][[1]], {apertureSamples - 1, 2 apertureSamples - 1}];
  H\[Eta]1 = Transpose[ConstantArray[Reverse[\[Eta]], 2 apertureSamples - 1]];
  H\[Eta]  = Transpose[ArrayFlatten[{{H\[Eta]1}, {H\[Eta]2}}]];
  Hy\[Eta] = Hy - H\[Eta];
  
  x\[Zeta]y\[Eta]coordArray = {Hx\[Zeta], Hy\[Eta]};
  x\[Zeta]y\[Eta]coordArray = Transpose[x\[Zeta]y\[Eta]coordArray, {3, 1, 2}];
  ClearAll[Hx1, Hx2, Hx, H\[Zeta]1, H\[Zeta]2,
           H\[Zeta]1, Hx\[Zeta], Hy1, Hy2, Hy,
           H\[Eta]1, H\[Eta]2, H\[Eta], Hy\[Eta]];
  
  (*Evaluate r across the observation plane*)
  rFun = Compile[{{x, _Real}, {y, _Real}}, Sqrt[x^2 + y^2 + z^2]];
  rEva = Map[rFun @@ # &, x\[Zeta]y\[Eta]coordArray, {2}];
  
  (*Evaluate gr*)
  H = Map[gr, rEva, {2}];
  
  (*Evaluate the fast Fourier transforms*)
  FFU  = Fourier[U];
  FFH  = Fourier[H];
  (*Convolution*)
  FFUH = FFU*FFH;
  (*Come back to real space*)
  S    = InverseFourier[FFUH];
  
  (*Nick just the good part*)
  field = \[CapitalDelta]\[Eta] * \[CapitalDelta]\[Zeta] * S[[apertureSamples ;;, apertureSamples ;;]];
  field = Transpose[field];

  (*What goes around comes around*)
  Return[{numSamples,x,y,field}];
  )

Options[VectorFieldProp] = {"numSamples" -> Automatic};
VectorFieldProp::usage = "VectorFieldProp takes the three cartesian components of a field in an aperture plane and propagates them to an observation plane by using an implementation of the direct integration of the Rayleigh-Sommerfeld diffraction integral. This implementation is based on the method described in Shen and Wang (2006). The field is sampled in the aperture plane using a uniform grid, and the observation plane is also sampled using a uniform grid. The field is assumed to be zero outside of the aperture plane. Contrary to the Fresnel or Kirchhoff integrals, the Rayleigh-Sommerfeld integral does not require the use of the paraxial approximation, and is therefore valid for arbitrary distances between the aperture plane and the observation plane. Here these two planes are assumed to be parallel. What is given here is only valid when the three components of the given fields can be propagated individually as if they were scalar fields.

Parameters
----------
+ Lobs (float): spatial width of the obsevation region, in μm. The observation region is assumed to be a squared centered on (x,y) = (0,0), and extending from -Lobs/2 to Lobs/2 in both the x and y directions.
+ z (float): distance between the aperture plane and the observation plane, given in μm. The aperture plane is assumed to be at z=0.
+ {apertureFunctionx, apertureFunctiony, apertureFunctionz} (List with three functions): a list with three two-variable functions which return the complex amplitude of the field in the aperture plane. Input to the functions is assumed to be in cartesian coordinates x,y. No attempt is made to make sure that these are valid electromagnetic fields.
+ λfree (float): wavelength in vaccum of field, given in μm.
+ nref  (float): refractive index of the propagating medium.

Options
-------
+ \"numSamples\" (int or Automatic): number of samples to use in the aperture plane and the observation plane. The aperture plane is sampled using a uniform grid, and the observation plane is also sampled using a uniform grid. The default is Automatic in which case numSamples is calculated so that the sample size is equal to half the wavelength of the wave inside of the propagating medium.

Returns
-------
{numSamples, xCoords, yCoords, {fieldx, fieldy, fieldz}}
+ xCoords (List): x coordinates of the observation plane, given in μm.
+ yCoords (List): y coordinates of the observation plane, given in μm.
+ {fieldx, fieldy, fieldz}   (List of List): each of the given lists are 2-dimensional lists which contain the complex amplitude of the cooresponding component of the field in the observation plane. The top left corner of each list corresponds to the lower left corner of the observation plane. The coordinates associated with each element in the given lists should be taken from xCoords and yCoords.

References
----------
+ Shen, Fabin, and Anbo Wang. \"Fast-Fourier-transform based numerical integration method for the Rayleigh-Sommerfeld diffraction formula.\" Applied optics 45, no. 6 (2006): 1102-1110.

"
VectorFieldProp[Lobs_, z_, {apertureFunctionx_, apertureFunctiony_, apertureFunctionz_}, \[Lambda]free_, nref_:1, OptionsPattern[]] := (
  {numSamplesx,xx,yx,fieldx} = ScalarFieldProp[Lobs, z, apertureFunctionx, \[Lambda]free, nref, "numSamples" -> OptionValue["numSamples"]];
  {numSamplesy,xy,yy,fieldy} = ScalarFieldProp[Lobs, z, apertureFunctiony, \[Lambda]free, nref, "numSamples" -> OptionValue["numSamples"]];
  {numSamplesz,xz,yz,fieldz} = ScalarFieldProp[Lobs, z, apertureFunctionz, \[Lambda]free, nref, "numSamples" -> OptionValue["numSamples"]];
  Return[{numSsamplesx, xx, yx, {fieldx, fieldy, fieldz}}];
  )

(* Begin["`Private`"]

End[]; *)

EndPackage[]