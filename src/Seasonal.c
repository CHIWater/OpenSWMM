//-----------------------------------------------------------------------------
//   Seasonal.c
//
//   Project:  EPA SWMM5 Extension
//   Version:  5.1
//   Date:     07/18/15   (Build 5.0.909)
//             03/22/16   (Build 5.1.910.1)
//   Author:   T. Xiao (CHI)
//
//   This module handles processing related to dynamicly changing 5 parameters of
// a subcatchment during simulation by multiplying it with a time pattern.
//
//	 Add support to three new time patterns for Horton infiltration in 5.1.910.1
//		(Max/min infiltration rate and decay constant)
//   and change "SHeadP" to "SectHeadP" for Green-Ampt infiltration.
//
//   The time pattern should be a MONTHLY time pattern.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "Seasonal.h"

char* SeasonalTypeWords[] =
    {
	"NPervP",
	"DSPervP",
	"SuctHeadP",				  //(OPENSWMM 5.1.910.1)
	"ConductP",
	"IDeficitP",
	"MAXINFRP",					  //(OPENSWMM 5.1.910.1)
	"MININFRP",					  //(OPENSWMM 5.1.910.1)
	"DECAYP",					  //(OPENSWMM 5.1.910.1)
	NULL};

TFile FSeasonal;

//=============================================================================
int seasonal_readParams(char* toks[], int ntoks)
//
//  Purpose: reads seasonal time pattern information from line of input data file
//  Input:   toks = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//
//  Format for first line that defines a Seasonal variatons is:
//    Sub_ID  Seasonal_Type	Time_Pattern
{
    int j, k, m;

    // --- check for minimum number of tokens
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that subcatchment exists in database
	j = project_findObject(SUBCATCH, toks[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, toks[0]);

    // --- check if second token is valid type of seasonal variation
    m = findmatch(toks[1], SeasonalTypeWords);
    if ( m < 0 )
        return error_setInpError(ERR_KEYWORD, toks[1]);
	
	// --- check that time pattern exists in database
	k = project_findObject(TIMEPATTERN, toks[2]);
    if ( k < 0 ) return error_setInpError(ERR_NAME, toks[2]);

	// --- check that time pattern is MONTHLY pattern
	if(Pattern[k].type != MONTHLY_PATTERN)
		return error_setInpError(ERR_SEASONAL, toks[2]);

	switch (m)
	{
		case 0:
			SubEx[j].NPat = k;
			break;
		case 1:
			SubEx[j].DSPat = k;
			break;
		case 2:
			if(InfilModel == GREEN_AMPT || InfilModel == MOD_GREEN_AMPT)
				SubEx[j].SPat = k;
			break;
		case 3:
			if(InfilModel == GREEN_AMPT || InfilModel == MOD_GREEN_AMPT)
				SubEx[j].KsPat = k;
			break;
		case 4:
			if(InfilModel == GREEN_AMPT || InfilModel == MOD_GREEN_AMPT)
				SubEx[j].IMDmaxPat = k;
			break;
		case 5:
			if(InfilModel == HORTON || InfilModel == MOD_HORTON)
				SubEx[j].MaxIRPat = k;
			break;
		case 6:
			if(InfilModel == HORTON || InfilModel == MOD_HORTON)
				SubEx[j].MinIRPat = k;
			break;
		case 7:
			if(InfilModel == HORTON || InfilModel == MOD_HORTON)
				SubEx[j].DecayPat = k;
			break;
		default:
			break;
	}
	
	return 0;
}

int getMonth()
//
//  Input: None
//  Output:  month of the current simulation time
//
{
	int month;
	DateTime currentDate;
	currentDate = getDateTime(NewRunoffTime);
	month = datetime_monthOfYear(currentDate) - 1;
	return month;
}

double seasonal_getOneParam(double defValue, int pattern)
//
//  Input:   defValue = default value of the Green-Ampt infiltration parameter (S, Ks or IMDmax)
//						or subarea dStore
//						or Horton infiltration parameter (Max infil rate, Min infil rate or decay) (5.1.910.1)
//			 pattern = index for pattern multiplier of the parameter
//  Output:  modified value
//  Purpose: apply time pattern factor if exists
//
{
	if(pattern >= 0)
	{
		double f;
		int month = getMonth();

		// Get pattern multiplier for the parameter in this month
		f = inflow_getPatternFactor(pattern, month, -1,-1);
		return defValue * f;
	}
	else
	{
		// Return default value
		return defValue;
	}
}

double seasonal_getalpha(TSubarea *subarea)
//
//  Input:   subarea = ptr. to subarea object (pervious area)
//  Output:  alpha
//  Purpose: apply N time pattern factor if exists
//
{
	if(SubEx[currSub].NPat >= 0)
	{
		int month;
		double nFactor;
		month = getMonth();

		// Get pattern multiplier for the N parameter in this month
		nFactor = inflow_getPatternFactor(SubEx[currSub].NPat, month, -1, -1);
		if(nFactor == 0)
			return subarea->alpha;
		else
			return subarea->alpha / nFactor;
	}
	else
	{
		return subarea->alpha;
	}
}

double seasonal_getL(TGrnAmpt *infil)
//
//  Input:   infil = ptr. to Green-Ampt infiltration object
//  Output:  Lu (depth of upper soil zone in ft)
//  Purpose: apply Ks time pattern factor
//
{
	if(SubEx[currSub].KsPat >= 0)
	{
		double f;
		int month;
		month = getMonth();

		// Get pattern multiplier for Ks (suction head)
		f = inflow_getPatternFactor(SubEx[currSub].KsPat, month, -1,-1);
		if(f > 0)
			return infil-> Lu * sqrt(f);
		else
			return infil->Lu;	// Return 0 caused devided by zero in grnampt_setT()
	}
	else
	{
		// Reture default L
		return infil->Lu;
	}
}

void seasonal_writeFile(TGrnAmpt *infil)
{
	DateTime currentDate;
	int y,m,d;
	char str[9];
    char dateStr[11 + 8 + 1];
	char * str2;
	// Check Green-Ampt infiltration
	if(InfilModel != GREEN_AMPT && InfilModel != MOD_GREEN_AMPT)
		return;
	// Check subcatchment
	if(Nobjects[SUBCATCH]==0)
		return;

    currentDate = getDateTime(NewRunoffTime);
	datetime_timeToStr(currentDate, str);
	datetime_decodeDate(currentDate,&y,&m,&d);
    sprintf(dateStr, "%02d/%02d/%4d %s", m,d, y, str);
	
	if ( !FSeasonal.file )
	{
		str2=Finp.name;
		str2=strcat(str2,".tsf");
		FSeasonal.file = fopen(str2,"wt");
		// Write header info		
		fprintf(FSeasonal.file, "IDs:	%s	%s	%s", Subcatch[0].ID, Subcatch[0].ID, Subcatch[0].ID);
		fprintf(FSeasonal.file, "\nDate/Time	FU	IMD	F");
		fprintf(FSeasonal.file, "\nMM/dd/yyyy	ft	frac	ft");
	}
	
    fprintf(FSeasonal.file, "\n%s", dateStr);
	fprintf(FSeasonal.file, "	%12.5f", infil->Fu);
	fprintf(FSeasonal.file, "	%12.5f", infil->IMD);
	fprintf(FSeasonal.file, "	%12.5f", infil->F);
}

void seasonal_Close()
{
	// Purpose: close seasonal time series file
	// called in swmm_end()

	if(FSeasonal.file)
	{
		fclose(FSeasonal.file);
		FSeasonal.file = NULL;
	}
}

void seasonal_createObjects()
{
	// Purpose: creates seasonal objects
	// called in project.createObjects()

	int j;
	SubEx = (TSubEx*) calloc(Nobjects[SUBCATCH], sizeof(TSubEx));
	currSub = 0;
	FSeasonal.file = NULL;
	for(j=0;j<Nobjects[SUBCATCH];j++)
	{
		// Initialize time pattern index for N and dStore
		SubEx[j].NPat = -1;					//(OPENSWMM 5.1.909)
		SubEx[j].DSPat = -1;				//(OPENSWMM 5.1.909)
												 
		// Initialize time pattern index for Green-Apmt infiltration			 
		SubEx[j].SPat = -1;					//(OPENSWMM 5.1.909)
		SubEx[j].KsPat = -1;				//(OPENSWMM 5.1.909)
		SubEx[j].IMDmaxPat = -1;			//(OPENSWMM 5.1.909)

		// Initialize time pattern index for Horton infiltration
		SubEx[j].MaxIRPat = -1;				//(OPENSWMM 5.1.910.1)
		SubEx[j].MinIRPat = -1;				//(OPENSWMM 5.1.910.1)
		SubEx[j].DecayPat = -1;				//(OPENSWMM 5.1.910.1)
	}
}
