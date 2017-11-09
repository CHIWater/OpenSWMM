//-----------------------------------------------------------------------------
//   WaterAge.c
//
//   Project:  Matric SWMM5
//   Version:  5.1
//   Date:     11/08/16  (Build 5.1.011.v1)
//   Author:   T. Xiao (CHI)
//
//   Water age modeling functions.
//
//   Build 5.1.011.v1:
//   - Model water age as a pollutant for each node and link in the network.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const double ZeroVolume = 0.0353147; // 1 liter in ft3

static int subWaterAgeIndex;		 // Water age index in subcatchment results
static int nodeWaterAgeIndex;        // Water age index in node results
static int linkWaterAgeIndex;        // Water age index in link results

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void setNodeOldState(int);
static void setLinkOldState(int);
static void addLinkWaterAge(int);
static void findLinkWaterAge(int, double);
static void findStorageWaterAge(int, double);
static void findNodeWaterAge(int);
static double getMixedWaterAge(double, double, double, double, double);

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
// void WaterAge_init(int, int, int); (called by output_open)
// void WaterAge_setOldState(void); (called by routing_execute in routing.c)
// void WaterAge_execute(double step); (called by routing_execute in routing.c)
// void WaterAge_getSubcatchWaterAge(float x[]); (called by subcatch_getResults)
// void WaterAge_getNodeWaterAge(int j, double f, double f1 float x[]); (called by node_getResults)
// void WaterAge_getLinkWaterAge(int j, double f, double f1 float x[]); (called by link_getResults)
// int WaterAge_readResetNodes(char* tok[], int ntoks) (called by input.parseLine())

void WaterAge_init(int subIndex, int nodeIndex, int linkIndex)
//  Purpose: set water age index for node/link results
{
	subWaterAgeIndex = subIndex - 1;
	nodeWaterAgeIndex = nodeIndex - 1;
	linkWaterAgeIndex = linkIndex - 1;
}

int WaterAge_readResetNodes(char* tok[], int ntoks)
{
	int j, t;
	for (t = 0; t < ntoks; t++)
	{
		j = project_findObject(NODE, tok[t]);
		if (j < 0)
			return error_setInpError(ERR_NAME, tok[t]);
//		Node[j].resetWaterAge = TRUE;
	}
	return 0;
}

void WaterAge_setOldState()
//  Purpose: replaces a node/link's old water age values with new ones.
{
	int j;
	for (j = 0; j < Nobjects[NODE]; j++)
		setNodeOldState(j);
	for (j = 0; j < Nobjects[LINK]; j++)
		setLinkOldState(j);
}

void setNodeOldState(int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: replaces a node's old water age values with new ones.
//
{
	Node[j].oldWaterAge = Node[j].newWaterAge;
	Node[j].newWaterAge = 0.0;
	Node[j].linkWaterAgeLoading = 0.0;
}

void setLinkOldState(int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: replaces link's old water age values with current ones.
//
{
	Link[j].oldWaterAge = Link[j].newWaterAge;
	Link[j].newWaterAge = 0.0;
}

void WaterAge_execute(double tStep)
//
//  Input:   tStep = routing time step (sec)
//  Output:  none
//  Purpose: routes water quality constituents through the drainage
//           network over the current time step.
//
{
	int    i, j;

	// --- add water age each link contributes to its downstream node
	for (i = 0; i < Nobjects[LINK]; i++)
		addLinkWaterAge(i);

	// --- find new water age at each node  
	for (j = 0; j < Nobjects[NODE]; j++)
	{
		// Check to reset incoming water age to zero
//		if (Node[j].resetWaterAge)
//			Node[j].linkWaterAgeLoading = 0;

		// --- find new quality at the node 
		if (Node[j].type == STORAGE || Node[j].oldVolume > FUDGE)
			findStorageWaterAge(j, tStep);
		else
			findNodeWaterAge(j);
	}

	// --- find new water age in each link
	for (i = 0; i < Nobjects[LINK]; i++)
		findLinkWaterAge(i, tStep);
}

void addLinkWaterAge(int i)
//
//  Input:   i = link index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: adds flow out of link to the water age at the link's downstream node.
//
{
	int    j;
	double qLink, w;

	// --- find inflow to downstream node
	qLink = Link[i].newFlow;

	// --- identify index of downstream node
	j = Link[i].node2;
	if (qLink < 0.0) j = Link[i].node1;
	qLink = fabs(qLink);

	// --- set node water age from link
	w = qLink * Link[i].oldWaterAge;
	Node[j].linkWaterAgeLoading += w;
}

void findNodeWaterAge(int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: finds new water age in a node with no storage volume.
//
{
	double qNode;

	// --- if there is flow into node then water age = sum(upstream link water age loading) / node flow
	qNode = Node[j].inflow;
	if (qNode > ZERO)
	{
		Node[j].newWaterAge = Node[j].linkWaterAgeLoading / qNode;
	}

	// --- otherwise water age is 0
	else
	{
		Node[j].newWaterAge = 0;
	}
}

void  findStorageWaterAge(int j, double tStep)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new water age in a node with storage volume.
//  
{
	double qIn,           // inflow rate at the storage node (cfs)
		wIn,              // water age loading from upstream links (hrs*cfs)
		v1,               // volume at start of time step (ft3)
		wa1,              // initial water age (hrs)
		wa2;              // final water age (hrs)

	// --- get inflow rate & initial volume
	qIn = Node[j].inflow;
	v1 = Node[j].oldVolume;

	// --- start with water age at start of time step 
	wa1 = Node[j].oldWaterAge;

	// --- mix resulting contents with inflow from all sources
	//     (temporarily accumulated in Node[j].newQual)
	wIn = Node[j].linkWaterAgeLoading;
	wa2 = getMixedWaterAge(wa1, v1, wIn, qIn, tStep);

	// --- set water age to zero if remaining volume is negligible
	if (Node[j].newVolume <= ZeroVolume)
		wa2 = 0.0;

	// --- assign new concen. to node
	Node[j].newWaterAge = wa2;
}

void findLinkWaterAge(int i, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new water age in a link at end of the current time step.
//
{
	int j,                // upstream node index
		k;                // conduit index
	double wIn,           // water age (hrs)
		qIn,              // inflow rate (cfs)
		v1,               // link volume at start of time step (ft3)
		v2,               // link volume at end of time step (ft3)
		wa1,              // current water age within link (hrs)
		wa2,              // new water age within link (hrs)
		qSeep,            // rate of seepage loss (cfs)
		vEvap,            // volume lost to evaporation (ft3)
		vLosses,          // evap. + seepage volume loss (ft3)
		barrels;          // number of barrels in conduit
	double a1;

	// --- identify index of upstream node
	j = Link[i].node1;
	if (Link[i].newFlow < 0.0)
		j = Link[i].node2;

	// --- link water age is that of upstream node when
	//     link is not a conduit or is a dummy link
	if (Link[i].type != CONDUIT || Link[i].xsect.type == DUMMY)
	{
		Link[i].newWaterAge = Node[j].newWaterAge;
		return;
	}

	// --- Steady Flow routing requires special treatment
	if (RouteModel == SF)
	{
		// set water age to that of upstream node
		Link[i].newWaterAge = Node[j].newWaterAge;
		return;
	}

	// --- get flow rates
	k = Link[i].subIndex;
	barrels = Conduit[k].barrels;
	qIn = fabs(Conduit[k].q1) * barrels;
	qSeep = Conduit[k].seepLossRate * barrels;
	vEvap = Conduit[k].evapLossRate * barrels * tStep;

	// --- get starting and ending volumes
	v1 = Link[i].oldVolume;
	v2 = Link[i].newVolume;
	vLosses = qSeep*tStep + vEvap;

	// --- adjust inflow to compensate for volume change under Dynamic
	//     Wave routing (which produces just a single (out)flow rate
	//     for a conduit)
	if (RouteModel == DW)
	{
		qIn = qIn + (v2 + vLosses - v1) / tStep;
		qIn = MAX(qIn, 0.0);
	}

	// --- start with water age at start of time step
	wa1 = Link[i].oldWaterAge;

	// --- mix resulting water age with inflow from upstream node
	wIn = Node[j].newWaterAge * qIn;

	// Check empty conduits
	// FUDGE is 0.0001 ft, which is too small
	// Use 0.5 mm instead
	if (Link[i].newDepth <= 0.0016)// && v2 < 0.035)
	{
		// Set new water age to 0 since SWMM5 doesn't allow empty conduit
		Link[i].newWaterAge = 0;
		// Set volume to 0 to allow quick response to water inflows
		Link[i].newVolume = 0;
		return;
	}

	wa2 = getMixedWaterAge(wa1, v1, wIn, qIn, tStep);

	// --- assign new water age to link
	Link[i].newWaterAge = wa2;
}

double getMixedWaterAge(double wa, double v1, double wIn, double qIn, double tStep)
//
//  Input:   wa = water age in reactor at start of time step (hrs)
//           v1 = volume in reactor at start of time step (ft3)
//           wIn = water age from upstream node/links (hrs)
//           qIn = flow inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  returns water age at end of time step (hrs)
//  Purpose: finds water age within a completely mixed reactor.
//
{
	double vIn, cIn, cMax;

	wa = wa + tStep / 3600.0;

	// --- if no inflow then reactor water age is added with time step
	if (qIn <= ZERO) return wa;

	// --- compute water age of any inflow
	vIn = qIn * tStep;
	cIn = wIn * tStep / vIn;

	// --- mixture water age can't exceed either original or inflow water age
	cMax = MAX(wa, cIn);

	// --- mix inflow with current reactor contents
	wa = (wa*v1 + wIn*tStep) / (v1 + vIn);
	wa = MIN(wa, cMax);
	wa = MAX(wa, 0.0);
	return wa;
}

void WaterAge_getSubcatchWaterAge(int j, float x[])
//
//  Input:   j = subcatchment index
//           x[] = array of nodal reporting variables
//  Output:  none
//  Purpose: computes weighted average of old and new water age at a subcatchment.
//
{
	float z;

	// Water age is zero for subcatchment
	z = 0;
	x[subWaterAgeIndex] = z;
}

void WaterAge_getNodeWaterAge(int j, double f, double f1, float x[])
//
//  Input:   j = node index
//           f, f1 = weighting factor
//           x[] = array of nodal reporting variables
//  Output:  none
//  Purpose: computes weighted average of old and new water age at a node.
//
{
	double z;

	z = f1*Node[j].oldWaterAge + f*Node[j].newWaterAge;
	x[nodeWaterAgeIndex] = (float)z;
}

void WaterAge_getLinkWaterAge(int j, double f, double f1, float x[])
//
//  Input:   j = link index
//           f, f1 = weighting factor
//           x[] = array of link reporting variables
//  Output:  none
//  Purpose: computes weighted average of old and new water age at a link.
//
{
	double z;

	z = f1*Link[j].oldWaterAge + f*Link[j].newWaterAge;
	x[linkWaterAgeIndex] = (float)z;
}