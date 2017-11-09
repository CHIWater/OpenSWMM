//-----------------------------------------------------------------------------
//   Seasonal.h
//
//   Project: EPA SWMM5 extension
//   Version: 5.1
//   Date:    07/18/15   (Build 5.1.909)
//   Author:  T. Xiao (CHI)
//
//   Public interface for Seasonal variations functions.
//
//-----------------------------------------------------------------------------

#ifndef Seasonal_H
#define Seasonal_H

#ifndef EXTERN
#define  EXTERN                        // defined as 'extern' in headers.h
#endif

// Extended attributes for Subcatchment
typedef struct
{
	int			NPat;		// Time pattern index for pervious area Manning's N			(OPENSWMM 5.1.909)
	int			DSPat;		// Time pattern index for pervious area	depression storage	(OPENSWMM 5.1.909)
	//----------Patterns for Green-Ampt infiltration-------------------														   
	int			SPat;		// Time pattern	index for Suction head						(OPENSWMM 5.1.909)
	int			KsPat;		// Time pattern	index for Conductivity						(OPENSWMM 5.1.909)
	int			IMDmaxPat;	// Time pattern	index for Initial deficit					(OPENSWMM 5.1.909)
	//----------Patterns for Horton infiltration-------------------														   
	int			MaxIRPat;	// Time pattern	index for Max infiltration rate				(OPENSWMM 5.1.910.1)
	int			MinIRPat;	// Time pattern	index for Min infiltration rate				(OPENSWMM 5.1.910.1)
	int			DecayPat;	// Time pattern	index for decay constant 					(OPENSWMM 5.1.910.1)
}  TSubEx;																				   
																						   
EXTERN TSubEx* SubEx;		// Array of extended attributes for subcatchments			(OPENSWMM 5.1.909)
EXTERN int currSub;			// Index for current subcatchment							(OPENSWMM 5.1.909)
							// Set in subcatch_getRunoff() in subcatch.c

// public methods
int    seasonal_readParams(char* tok[], int ntoks);
double seasonal_getOneParam(double defValue, int pattern);
double seasonal_getalpha(TSubarea *subarea);
double seasonal_getL(TGrnAmpt *infil);
void seasonal_writeFile(TGrnAmpt *infil);
void seasonal_Close();
void seasonal_createObjects();

#endif