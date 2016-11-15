/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-cyclocopter/module-cyclocopter.cc,v 1.5 2016/02/17 15:41:37 masarati Exp $ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/*
 * Author: Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>
 *
 * Reworked as runtime loadable module by
 * Pierangelo Masarati <pierangelo.masarati@polimi.it>
 */

/* Rotor elements */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>
#include <map>

#include "dataman.h"
#include "userelem.h"
#include "indvel.h"

// For CyclocopterDMST
#include <vector>
#include "aeroelem.h"
#include <complex>

/* CyclocopterInflow - begin */

/* Base class for cycloidal rotor  inflow models:
 * all models a derived from this class
*/

class CyclocopterInflow
: virtual public Elem, public UserDefinedElem, public InducedVelocity {
protected:
	const StructNode* pRotor;

	bool bFlagAverage;
	// == true: uses the mean of the forces over a cycle
	//	to calculate the induced velocity on the next rotation
	// == false: uses an instant value (possibly filtrated)

	doublereal dRadius;		// Rotor radius
	doublereal dSpan;		// Blade length
	doublereal dArea;		// Cylinder longitudinal area

	doublereal dKappa;		// Hover correction coefficient
	doublereal dOmega;		// Rotor Angular velocity

	DriveOwner Weight;
	doublereal dWeight;

	Mat3x3 RRot;

	/* data used for output calculation */
	Mat3x3 RRotorTranspose;
	doublereal dUindMean;

	// Filter coefficients for the second order butteworth filter
	doublereal a1, a2, b0, b1, b2;

	// Filter coefficients for the second order butteworth filter
	void SetFilterCoefficients(doublereal dOmegaFilter, doublereal dDeltaT);

public:
	CyclocopterInflow(unsigned int uL, const DofOwner* pDO);
	virtual ~CyclocopterInflow(void);

	virtual Elem::Type GetElemType(void) const;
	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborate the internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// Contribution to the restart file
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);

	// Relative to the ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// *******FOR PARALLEL SOLVER********
	// Provides the type and label of the nodes connected to the element
	// useful for the assembly of the DOF connexion matrix
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

CyclocopterInflow::CyclocopterInflow(unsigned int uL, const DofOwner* pDO)
: Elem(uL, flag(0)),
UserDefinedElem(uL, pDO),
InducedVelocity(uL, 0, 0, flag(0)),
pRotor(0),
bFlagAverage(false),
dRadius(0.),
dSpan(0.),
dArea(0.),
Weight(0),
dWeight(0.),
RRot(::Eye3),
RRotorTranspose(::Eye3),
dUindMean(0.),
a1(0.), a2(0.), b0(1.), b1(0.), b2(0.)
{
	NO_OP;
}

CyclocopterInflow::~CyclocopterInflow(void)
{
	NO_OP;
}

Elem::Type
CyclocopterInflow::GetElemType(void) const
{
	return Elem::LOADABLE;
}

InducedVelocity::Type
CyclocopterInflow::GetInducedVelocityType(void) const
{
	return InducedVelocity::USER_DEFINED;
}

void
CyclocopterInflow::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterInflow::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << "0."                	 /* 9 */
                        << " " << "0."                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << "0."                	 /* 14 */
                        << " " << "0."   		 /* 15 */
                        << " " << "0."           	 /* 16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

std::ostream&
CyclocopterInflow::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

void
CyclocopterInflow::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

void
CyclocopterInflow::SetInitialValue(VectorHandler& X)
{
	NO_OP;
}

void
CyclocopterInflow::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pCraft;
	connectedNodes[1] = pRotor;
}

unsigned int
CyclocopterInflow::iGetInitialNumDof(void) const
{
	return 0;
}

void 
CyclocopterInflow::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
CyclocopterInflow::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
CyclocopterInflow::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

void
CyclocopterInflow::SetFilterCoefficients(doublereal dOmegaFilter,
	doublereal dDeltaT)
{
	/* Butterworth discrete low-pass filter coefficients */
	if (dDeltaT > 0. && dOmegaFilter > 0.) {
		doublereal dTmp = 4. + 2.*sqrt(2.)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter;
		a1 = (-8. + 2.*dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;
		a2 = (4. - 2.*sqrt(2.)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;

		dTmp = dOmegaFilter*dOmegaFilter*dDeltaT*dDeltaT/dTmp;
		b0 = dTmp;
		b1 = 2.*dTmp;
		b2 = dTmp;
	}
}

static bool
ReadRotorData(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel,
	const StructNode *& pCraft,
	Mat3x3& rrot,
	const StructNode *& pRotor)
{
     	/* aircraft node */
	pCraft = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	/* rotor orientation with respect to aircraft */
     	rrot = ::Eye3;
     	if (HP.IsKeyWord("orientation")) {
     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);

     	} else if (HP.IsKeyWord("hinge")) {
		silent_cerr("InducedVelocity(" << uLabel << "): deprecated keyword \"hinge\"; use \"orientation\" instead at line " << HP.GetLineData() << std::endl);

     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);
     	}

     	/* rotor node */
     	pRotor = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	return true;
}

static bool
ReadUniform(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel,
	bool& bFlagAve,
	doublereal& dR,
	doublereal& dL,
	DriveCaller *& pdW,
	doublereal& dOmegaFilter,
	doublereal& dKappa,
	doublereal& dDeltaT,
	doublereal& dOmega)
{
	bFlagAve = HP.GetYesNoOrBool();

	dR = HP.GetReal();
	if (dR <= 0.) {
		silent_cerr("ReadUniform(" << uLabel << "): "
			"illegal null or negative radius"
			"for rotor" << uLabel << " at line " << HP.GetLineData()
			<< std::endl);
		return false;
	}

	dL = HP.GetReal();
	if (dL <= 0.) {
		silent_cerr("ReadUniform(" << uLabel << "): "
			"illegal null or negative blade"
			"length for rotor" << uLabel
			<< " at line " << HP.GetLineData()
			<< std::endl);
		return false;
	}

	pdW = 0;
	if (HP.IsKeyWord("delay")) {
		pdW = HP.GetDriveCaller();

	} else {
		SAFENEW(pdW, NullDriveCaller);
	}

	dOmegaFilter = 0.;
	if (HP.IsKeyWord("omegacut")) {
		dOmegaFilter = HP.GetReal();
		if (dOmegaFilter <= 0) {
			silent_cerr("Illegal null or negative filter"
				"cut frequency for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}
	} else {
		dOmegaFilter = 0.;
	}

	dKappa = 1.;
	if (HP.IsKeyWord("kappa")) {
		dKappa = HP.GetReal();
		if (dKappa <= 0) {
			silent_cerr("Illegal null or negative hover"
				"correction coefficient" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}
	}

	dDeltaT = 0.;
	if (HP.IsKeyWord("timestep")) {
		dDeltaT = HP.GetReal();
		if (dDeltaT <= 0) {
			silent_cerr("Illegal null or negative time"
				"step for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}

	} else {
		dDeltaT = 0.;
	}

	dOmega = 1.;
	if (HP.IsKeyWord("omega")) {
		dOmega = HP.GetReal();
		if (dOmega <= 0) {
			silent_cerr("Illegal null or negative omega"
			    << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}
	}

	return true;
}

/* CyclocopterInflow - end */


/* CyclocopterNoInflow - begin */

class CyclocopterNoInflow
: virtual public Elem, public CyclocopterInflow {
public:
	CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterNoInflow(void);

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restores the induced velocity to an element
	// based on the azimuth position (used to
	// iterare between induced velocity calculation and
	// calculation of aerodynamic forces
	// for the KARI induced velocity model of the cycloidal rotor)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restores the induced velocity dalla metà superiore del rotore
	// ( only for KARI cycloidal rotor)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	}
	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	}
	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	}

};

CyclocopterNoInflow::CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal no inflow ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));
}

CyclocopterNoInflow::~CyclocopterNoInflow(void)
{
	NO_OP;
}

SubVectorHandler&
CyclocopterNoInflow::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	if (bToBeOutput()) {
		RRotorTranspose = pCraft->GetRCurr()*RRot;
		RRotorTranspose = RRotorTranspose.Transpose();
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterNoInflow::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Calculates the moment only if output is required */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}
}

Vec3
CyclocopterNoInflow::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{

	return Zero3;
}

/* CyclocopterNoInflow - end */

/* CyclocopterUniform1D - begin */

/*

From Moble: the induced velocity is opposite to the
the force generated by the rotor in the direction 3
of the rotor reference (the direction 1 of this
reference must be align with the rotation axis of
the rotor). It should make sense just for a multiblade
rotor (not for the one-blade rotor) when the generated
force is mainly in one direction.

*/

class CyclocopterUniform1D
: virtual public Elem, public CyclocopterInflow {
protected:
	Vec3 RRot3;
	Mat3x3 RRotor;

	mutable doublereal dUindMeanPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	doublereal dTz, dTzMean;
	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* data for force filtering */
	doublereal Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterUniform1D(void);

	// Elaborate internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restores the induced velocity to an element
	// based on the azimuth position (used to
	// iterare between induced velocity calculation and
	// calculation of aerodynamic forces
	// for the KARI induced velocity model of the cycloidal rotor)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restores the induced velocity from the higher half of the rotor
	// ( only for KARI cycloidal rotor)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterUniform1D::CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRot3(::Zero3),
RRotor(::Eye3),
dUindMeanPrev(0.),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
dTz(0.), dTzMean(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(0.), Uk_1(0.), Uk_2(0.), Yk(0.), Yk_1(0.), Yk_2(0.)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal uniform 1D ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT, dOmega)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterUniform1D::~CyclocopterUniform1D(void)
{
	NO_OP;
}

void
CyclocopterUniform1D::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << FMeanOut                	 /* 14 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

void
CyclocopterUniform1D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
	/* calculates the mean of the forces generated by the rotor over a cycle*/
	dTzMean += dTz;
	FMean += F;
	iStepCounter++;
	// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean /= iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			dTzMean = dTzMean/iStepCounter;
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*copysign(std::sqrt(std::abs(dTzMean)/(2*dRho*dArea)), dTzMean);
			dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
			dTzMean = 0.;
		}
		FMean = ::Zero3;
		iStepCounter = 0;
	}
	dAzimuthPrev = dAzimuth;

	/* update the inputs and outputs of the filter */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dUindMeanPrev = dUindMean;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;
	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

SubVectorHandler&
CyclocopterUniform1D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity (Moble version)*/
	/* Transpose of the rotor rotation matrix */
	RRotor = pCraft->GetRCurr()*RRot;
	RRot3 = RRotor.GetVec(3);
	RRotorTranspose = RRotor.Transpose();
	/* Force in the rotor coordinate system */
	F = RRotorTranspose*Res.Force();
	dTz = RRot3*Res.Force();
	if (!bFlagAverage) {
		/* filter the forces */
		Uk = dTz;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		dTz = Yk;
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*copysign(std::sqrt(std::abs(dTz)/(2*dRho*dArea)), dTz);

		dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterUniform1D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{

	/* Calculates the azimuth position of the first blade */
	if (bFlagIsFirstBlade == true) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}


	/* Calculates the moment only if output is required */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform1D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return RRot3*dUindMean;
}

/* CyclocopterUniform1D - end */

/* CyclocopterUniform2D - begin */

/*

Uniform inflow: the induced velocity is opposite to
the force generated by the rotor in the plane
perdendicular to the rotor rotation axis. The
rotor reference must have direction 1 aligned
with the rotor rotation axis!

*/

class CyclocopterUniform2D
: virtual public Elem, public CyclocopterInflow {
protected:
	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* data for force filtering */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterUniform2D(void);

	// Elaborate internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restores the induced velocity to an element
	// based on the azimuth position (used to
	// iterare between induced velocity calculation and
	// calculation of aerodynamic forces
	// for the KARI induced velocity model of the cycloidal rotor)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restores the induced velocity dalla metà superiore del rotore
	// ( only for KARI cycloidal rotor)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterUniform2D::CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRotor(::Eye3),
dUind(::Zero3), dUindPrev(::Zero3),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(::Zero3), Uk_1(::Zero3), Uk_2(::Zero3), Yk(::Zero3), Yk_1(::Zero3), Yk_2(::Zero3)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal uniform 2D ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT, dOmega)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterUniform2D::~CyclocopterUniform2D(void)
{
	NO_OP;
}

void
CyclocopterUniform2D::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << dUind                	 /* 11-13 */
                        << " " << FMeanOut                	 /* 14-16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}

}

void
CyclocopterUniform2D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
#if 0
	if (bFlagAverage) {
		/* calculates the mean of the forces generated by the rotor over a cycle*/
		FMean = FMean + F;
		iStepCounter++;
		// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
		if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
			FMean = FMean/iStepCounter;
			/* Force in the plane normal to the rotation axis */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Induced velocity: calculated according to dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = sqrt(dT/(2*dRho*dArea));
			/* Induced velocity components in the coordinate
	 		* system of the rotor */
			dUind = 0.;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

			FMean = 0.;
			iStepCounter = 0;
		}
	}
#endif

	/* calculates the mean of the forces generated by the rotor over a cycle*/
	FMean = FMean + F;
	iStepCounter++;
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			/* Force in the plane normal to the rotation axis */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Induced velocity: calculated according to dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
			/* Induced velocity components in the coordinate
	 		* system of the rotor */
			dUind = ::Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
		}

		FMean = ::Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* update the inputs and outputs of the filter */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;

	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

SubVectorHandler&
CyclocopterUniform2D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Transpose of the rotor rotation matrix */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Force in the rotor coordinate system */
	F = RRotorTranspose*Res.Force();
	if (!bFlagAverage) {
		/* filter the forces */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;
		/* Force in the plane normal to the rotation axis */
		doublereal dT = sqrt(F(2)*F(2) + F(3)*F(3));
		/* Induced velocity: calculated according to dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
		/* Induced velocity components in the coordinate
	 	* system of the rotor */
		dUind = ::Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterUniform2D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Calculates the azimuth position of the first blade */
	// if (bFlagIsFirstBlade && bFlagAverage) {
	if (bFlagIsFirstBlade) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}

	/* Calculates the moment only if output is required */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform2D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return RRotor*dUind;
}
/* CyclocopterUniform2D - end */


/* CyclocopterPolimi - begin */
/*
The induced velocity is opposite to
the force generated by the rotor in the plane
perdendicular to the rotor rotation axis. The
rotor reference must have direction 1 aligned
with the rotor rotation axis!
The inflow velocity distribution is constant along
the cylinder span and is function of r:
Vi(r) = Kc*cos((pi/2)*(r/R)) + Ks*sin(pi*(r/R))
where the cofficients Kc and Ks are based on the
momentum theory
*/

class CyclocopterPolimi
: virtual public Elem, public CyclocopterInflow {
protected:
	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	doublereal dXi;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* data for force filtering */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

	unsigned int iCounter;
	unsigned int iRotationCounter;

	doublereal dpPrev, dp;

	bool flag_print;

public:
	CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterPolimi(void);

	// Elaborate internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restores the induced velocity to an element
	// based on the azimuth position (used to
	// iterare between induced velocity calculation and
	// calculation of aerodynamic forces
	// for the KARI induced velocity model of the cycloidal rotor)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restores the induced velocity dalla metà superiore del rotore
	// ( only for per KARI cycloidal rotor)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterPolimi::CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRotor(::Eye3),
dUind(::Zero3), dUindPrev(::Zero3),
dXi(0.),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(::Zero3), Uk_1(::Zero3), Uk_2(::Zero3), Yk(::Zero3), Yk_1(::Zero3), Yk_2(::Zero3),
iCounter(0), iRotationCounter(0),
dpPrev(0.), dp(0.),
flag_print(true)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal Polimi ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT, dOmega)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterPolimi::~CyclocopterPolimi(void)
{
	NO_OP;
}

void
CyclocopterPolimi::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
	/* calculates the mean of the forces generated by the rotor over a cycle*/
	FMean = FMean + F;;
	iStepCounter++;
	// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			/* Force in the plane normal to the rotation axis */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Induced velocity: calculated according to dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
			/* Induced velocity components in the coordinate
	 		* system of the rotor */
			dUind = Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
			/* angle by which the tension is rotated */
			dXi = atan2(FMean(3), FMean(2)) - M_PI/2.;
		}

		FMean = Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* update the inputs and outputs of the filter */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;

	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

void
CyclocopterPolimi::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean              /* 8 */
                        << " " << dUind                	 /* 9 -11*/
                        << " " << dXi                	 /* 12 */
                        << " " << dAzimuth               /* 13 */
                        << " " << FMeanOut                	 /* 14-16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

SubVectorHandler&
CyclocopterPolimi::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Transpose of the rotor rotation matrix */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Force in the rotor coordinate system */
	F = RRotorTranspose*Res.Force();
	if (!bFlagAverage) {
		/* filter the forces */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;
		/* Force in the plane normal to the rotation axis */
		doublereal dT = sqrt(F(2)*F(2) + F(3)*F(3));
		/* Induced velocity: calculated according to dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
		/* Induced velocity components in the coordinate
	 	* system of the rotor */
		dUind = Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
		/* angle by which the tension is rotated */
		dXi = atan2(F(3), F(2)) - M_PI/2.;
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterPolimi::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Calculates the azimuth position of the first blade */
	if (bFlagIsFirstBlade) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}

	/* Calculates the moment only if output is required */
	if (bToBeOutput()) {
		Res.AddForces(F,M,X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterPolimi::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));

	doublereal d1 = XRel(2);
	doublereal d2 = XRel(3);

	/* dPsi0 is useless because the relative
	 * angle is used: (dp-dXi)!!! */
	doublereal dpp = atan2(d2, d1);

	doublereal r = sqrt(d1*d1 + d2*d2)*cos(dpp - dXi);

	return RRotor*((dUind*(M_PI/2.))*cos((M_PI/2.)*(r/dRadius)));
}

/* CyclocopterPolimi - end */


/*AUTHOR: Kuldeep Singh <kuldeepsingh050895@gmail.com>
Copyright (C) 2016(-2016) all rights reserved.
The copyright of this patch(For the class CyclocopterDMST or DMST Module) is transferred
to Pierangelo Masarati and Paolo Mantegazza
for use in the software MBDyn as described 
in the GNU Public License version 2.1 */

/* CyclocopterDMST - begin */

/*
Double Multiple Stream Tube inflow Model: 
The induced velocity at a particular azimuth location is opposite to 
the "blade force at that azimuth location" and in the plane
perpendicular to the rotor rotation axis. The
rotor reference must have direction 1 aligned
with the rotor rotation axis!
*/

class CyclocopterDMST
: virtual public Elem, public CyclocopterInflow {
protected:
	Mat3x3 RRotor;
	
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	int NBlades;
	std::vector<unsigned> uBlade; // to store the blade labels
	std::vector<const Elem*> pBlade; 
	std::vector<const AerodynamicBody*> pBladeAero; // stores the blade pointers

	std::vector<Vec3> dUindBlade;
	mutable std::vector<Vec3> dUindPrevBlade;

	std::vector<Vec3> dBladeForce; // Blade force in global reference frame

	// This stores the induce velocity magnitude at each azimuth location
	// this will be updated after every iteration to keep the it updated
	// Length of this array will be (2*PI/(dOmega*dDeltaT)), which is equal to twice of stream tubes  
	std::vector<doublereal> dIndVelMagAll; 
	
	std::vector<Vec3> dBladePos; // Blade position in global reference frame
	std::vector<doublereal> dAngleBtBladeAndNetforce;  // angle between Blade position with respect to Net rotor force vector 
	std::vector<doublereal> dBladeUindMag; // Blade induced velocity magnitude 

	// induced velocity vectors at each blade 
	std::vector<Vec3> dBladeUind; 
	std::vector<Vec3> dBladeUindPrev; 

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;
	doublereal dDeltaT;

	// Will make sure there are dIndVelMagAll do not have empty elements
	bool EmptyCheck; // Initially dIndVelMagAll is empty

	/* data for force filtering */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterDMST(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterDMST(void);

	// Elaborate internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);
	
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);


	virtual doublereal IndVelMag(int i, doublereal dRho);

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restores the induced velocity to an element
	// based on the azimuth position (used to
	// iterate between induced velocity calculation and
	// calculation of aerodynamic forces
	// for the KARI induced velocity model of the cycloidal rotor)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// ( only for KARI cycloidal rotor)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterDMST::CyclocopterDMST(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRotor(::Eye3),
dUind(::Zero3), dUindPrev(::Zero3),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(::Zero3), Uk_1(::Zero3), Uk_2(::Zero3), Yk(::Zero3), Yk_1(::Zero3), Yk_2(::Zero3),
EmptyCheck(true)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Kuldeep Singh <kuldeepsingh050895@gmail.com>	\n"
"Mentored by: Louis Gagnon <louis.gagnon@polimi.it> \n" "based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models	\n"
"		for cycloidal rotors using Double Multiple Stream Tube method. \n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal DMST ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT, dOmega)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	int BladeLabel = 1 ;
	while(BladeLabel)
	{
		BladeLabel = HP.GetInt();
		if (BladeLabel != 0)
		{
			uBlade.push_back(BladeLabel);
		}
	}
	NBlades = uBlade.size();

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterDMST::~CyclocopterDMST(void)
{
	NO_OP;
}

void
CyclocopterDMST::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	// Fixing the size of pointer storage vectors
	pBlade.resize(NBlades);
	pBladeAero.resize(NBlades);
	
	// Assigning the pointers
	for (int i = 0; i < NBlades; ++i)
	{
		pBlade[i] = pDM->pFindElem(Elem::AERODYNAMIC,uBlade[i]);
		pBladeAero[i] = dynamic_cast< const AerodynamicBody *>(pBlade[i]);
	}

	//Defining the size to store the corresponding values
	dBladeForce.resize(NBlades); 
	dBladePos.resize(NBlades); 
	dAngleBtBladeAndNetforce.resize(NBlades); 
	dBladeUind.resize(NBlades);
	dBladeUindPrev.resize(NBlades);
	dBladeUindMag.resize(NBlades);
	dIndVelMagAll.resize(floor(2*M_PI/(dOmega*dDeltaT)));
}

void
CyclocopterDMST::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   		/* 1 */
                        << " " << RRotorTranspose*Res.Force()   /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()  /* 5-7 */
                        << " " << dUindMean                	 	/* 8 */
                        << " " << iStepCounter             	 	/* 9 */
                        << " " << dAzimuth               		/* 10 */ // First blade azimuth location
                        << " " << dBladeUind[0]            	 	/* 11-13 */ // Induced velocity at first blade
                        << " " << FMeanOut                	 	/* 14-16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}

}

doublereal 
CyclocopterDMST::IndVelMag(int i, doublereal dRho)
{
	return sqrt((dKappa*dBladeForce[i].Norm()*NBlades*pow(sin((M_PI/2)-dAngleBtBladeAndNetforce[i]),2))/(4*M_PI*dRho*dRadius)); 
}

void
CyclocopterDMST::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;

	/* calculates the mean of the forces generated by the rotor over a cycle*/
	FMean = FMean + F;
	iStepCounter++;

	// Once the array dIndVelMagAll is full switch to full DMST model 
	if (iStepCounter > floor(2*M_PI/(NBlades*dOmega*dDeltaT)) )
	{	
		EmptyCheck = false;
	}

	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) 
	{
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		FMean = ::Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	for (int i = 0; i < NBlades; ++i)
	{
		// After convergence storing the converged value of induced velocity
		dBladeUindPrev[i] = dBladeUind[i]; 
	}

	/* update the inputs and outputs of the filter */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;

	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

SubVectorHandler&
CyclocopterDMST::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Transpose of the rotor rotation matrix */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();

	/* Force in the rotor coordinate system */
	F = RRotorTranspose*Res.Force();
	/* filter the forces */
	Uk = F;
	Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
	F = Yk;

	///////////////////////////////////////////////////////////////
	// 2D Model
	/* Force in the plane normal to the rotation axis */
	doublereal dT = sqrt(F(2)*F(2) + F(3)*F(3));
	/* Induced velocity: calculated according to dT */
	doublereal dRho = dGetAirDensity(GetXCurr());
	dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
	/* Induced velocity components in the coordinate
 	* system of the rotor */
	dUind = ::Zero3;
	if (dT > std::numeric_limits<doublereal>::epsilon()) 
	{
		dUind(2) = dUindMean*F(2)/dT;
		dUind(3) = dUindMean*F(3)/dT;
	}
	dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
	dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
	dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

	dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));


	//////////////////////////////////////////////////////////////////////
	// DMST Model
	Vec3 NetForce = RRotorTranspose*Res.Force(); //(0, Fx, Fy) Net rotor force

	for (int i = 0; i < NBlades; ++i)
	{
		dBladeForce[i] = RRotorTranspose*pBladeAero[i]->dGetForces(); // storing the all blade force 
		dBladePos[i] = RRotorTranspose*pBladeAero[i]->dGetPosition(); // storind the position of blades
		dAngleBtBladeAndNetforce[i] = acos(dBladePos[i].Dot(NetForce)/sqrt(dBladePos[i].Dot(dBladePos[i] * NetForce.Dot(NetForce))));
		
		// At first time step NetForce is (0,0,0) but it should be nonzero vector because angle with respect to null vector is undefined 
		if (NetForce.Norm() == 0)
		{	
			dAngleBtBladeAndNetforce[i] = 2*i*M_PI/NBlades;
		} 
	}

	// Updating the dIndVelMagAll  
	for (int i = 0; i < NBlades; ++i)
	{
		 doublereal psi = atan2(dBladePos[i](3),dBladePos[i](2)); // Azimuth location of the blade i 
		 if (psi < 0)
		 {
		 	psi = psi + 2*M_PI; // atan2() range (-PI, PI) so to make it (0, 2*PI) => -psi = psi + 2*PI 
		 }
		 int index = floor(psi/(dOmega*dDeltaT)); // where to change the induced velocity
		 dIndVelMagAll[index] = IndVelMag(i,dRho);
	}

	////////////////////////////////
	// Upper half and Lower half both 
	// Because until we have nonzero U_up, we can not get nonzero "w"  
	// and without nonzero "w" the root method will not give U_d 
	////////////////////////////////	
	#if 1
	if(EmptyCheck == true)
	{
		// Calculating magnitude of Induced velocity of individual blade
		for (int i = 0; i < NBlades; ++i)
		{	
			// Induced velocity magnitude is calculated using the blade force of the individual blade and the Azimuth location of the blade
			dBladeUindMag[i] = IndVelMag(i,dRho);

			// Calculating the Induced velocity component for Upper half
			dBladeUind[i] = ::Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon())
			{
				if (dAngleBtBladeAndNetforce[i] > 0) // upper half 
				{
					// Option 1: Blade Azimuth location
					dBladeUind[i](2) = -dBladeUindMag[i]*dBladePos[i](2)/dBladePos[i].Norm();
					dBladeUind[i](3) = -dBladeUindMag[i]*dBladePos[i](3)/dBladePos[i].Norm();					
				}
				else /*dAngleBtBladeAndNetforce[i] < 0*/ // Lower half
				{
					// Option 1: Blade Azimuth location
					dBladeUind[i](2) = dBladeUindMag[i]*dBladePos[i](2)/dBladePos[i].Norm();
					dBladeUind[i](3) = dBladeUindMag[i]*dBladePos[i](3)/dBladePos[i].Norm();
				}

				// Option 2: Local blade net force vector direction
				// dBladeUind[i](2) = dBladeUindMag[i]*dBladeForce[i](2)/dBladeForce[i].Norm();
				// dBladeUind[i](3) = dBladeUindMag[i]*dBladeForce[i](3)/dBladeForce[i].Norm();
				
				// Option 3: Rotor, Net thrust vector direction
				// dBladeUind[i](2) = dBladeUindMag[i]*F(2)/dT;
				// dBladeUind[i](3) = dBladeUindMag[i]*F(3)/dT;

				// Option 4: Given in the Article BRNDICT ET AL   
				// dBladeUind[i](3) = -dBladeUindMag[i];
			
			}
			dBladeUind[i](1) = (1 - dWeight)*dBladeUind[i](1) + dWeight*dBladeUindPrev[i](1);
			dBladeUind[i](2) = (1 - dWeight)*dBladeUind[i](2) + dWeight*dBladeUindPrev[i](2);
			dBladeUind[i](3) = (1 - dWeight)*dBladeUind[i](3) + dWeight*dBladeUindPrev[i](3);
		}
	}
	#endif

	if (EmptyCheck == false)
	{
		// Calculating magnitude of Induced velocity of individual blade
		for (int i = 0; i < NBlades; ++i)
		{	
			////////////////////////////////
			// Upper half
			////////////////////////////////	
			if (cos(dAngleBtBladeAndNetforce[i]) > 0) // Angle between NetForce and Blade position is in range (-PI/2, PI/2) 
			{
				// Induced velocity magnitude is calculated using the blade force of the individual blade and the Azimuth location of the blade
				dBladeUindMag[i] = IndVelMag(i,dRho);

				// Calculating the Induced velocity component for Upper half
				dBladeUind[i] = ::Zero3;
				if (dT > std::numeric_limits<doublereal>::epsilon())
				{
					// Option 1: Blade Azimuth location
					dBladeUind[i](2) = -dBladeUindMag[i]*dBladePos[i](2)/dBladePos[i].Norm();
					dBladeUind[i](3) = -dBladeUindMag[i]*dBladePos[i](3)/dBladePos[i].Norm();

					// Option 2: Local blade net force vector direction
					// dBladeUind[i](2) = -dBladeUindMag[i]*dBladeForce[i](2)/dBladeForce[i].Norm();
					// dBladeUind[i](3) = -dBladeUindMag[i]*dBladeForce[i](3)/dBladeForce[i].Norm();
					
					// Option 3: Rotor, Net thrust vector direction
					// dBladeUind[i](2) = dBladeUindMag[i]*F(2)/dT;
					// dBladeUind[i](3) = dBladeUindMag[i]*F(3)/dT;

					// Option 4: Given in the Article BRNDICT ET AL   
					// dBladeUind[i](3) = -dBladeUindMag[i];	
				}

				dBladeUind[i](1) = (1 - dWeight)*dBladeUind[i](1) + dWeight*dBladeUindPrev[i](1);
				dBladeUind[i](2) = (1 - dWeight)*dBladeUind[i](2) + dWeight*dBladeUindPrev[i](2);
				dBladeUind[i](3) = (1 - dWeight)*dBladeUind[i](3) + dWeight*dBladeUindPrev[i](3);
			}
			/////////////////////////////////
			// Lower half
			/////////////////////////////////
			else // (dAngleBtBladeAndNetforce[i]) < 0)
			{
				doublereal w;
				int MirrorPointIndex; // Lower half blade position's mirror image point in the upper half 
				// VECTOR ABOUT WHICH MIRROR IMAGE IS TAKEN, IS VECTOR ORTHOGNAL TO THE NetForce Vector which is = (0, -NetForce(3), NetForce(2)) 

				Vec3 NetForceOrtho = ::Zero3; // Vector Orthogonal to Net force Vector 
				NetForceOrtho(2) = NetForce(3);  
				NetForceOrtho(3) = -NetForce(2);  
				
				doublereal m = (NetForceOrtho(3)/NetForceOrtho(2)); // Slop of the line

				doublereal a = dBladePos[i](2) ; 
				doublereal b = dBladePos[i](3) ;
				
				doublereal c = ((2*(m*b+a))/(m*m+1)) - a ; // Mirror image point first coordinate  
				doublereal d = ((2*m*(m*b+a))/(m*m+1)) - b ; // Mirror image point second coordinate 


				doublereal MirrorAzimuth = atan2(d,c); // Mirror image point's Azimuth angle in the global frame 
				if (MirrorAzimuth < 0)
				{
					MirrorAzimuth = MirrorAzimuth + 2*M_PI ; 
				}
				MirrorPointIndex = floor(MirrorAzimuth / (dOmega*dDeltaT));
				

				w = 2*dIndVelMagAll[MirrorPointIndex]/sin(-dAngleBtBladeAndNetforce[i] + M_PI/2);

				//////////////////////////////////////////
				//Find induced velocity using roots   
				#if 1
				// Direct root calculation method to find the induced velocity in lower-half of the rotor
				float sinPsi = sin(-dAngleBtBladeAndNetforce[i] + M_PI/2);
				float Gamma = pow(((dKappa*dBladeForce[i].Norm()*NBlades)/(4*M_PI*dRho*dRadius)),2);
				/*x1, x2, x3 and x4 are the roots of the 4th order polynomial*/
				std::complex<double> x1 = 0+0j; 	
				std::complex<double> x2 = 0+0j; 	
				std::complex<double> x3 = 0+0j; 	
				std::complex<double> x4 = 0+0j; 

				std::complex<double> c1 = 1 + 0j; 	
				std::complex<double> c2 = 2 + 0j; 	
				std::complex<double> c3 = 3 + 0j; 	
				std::complex<double> c4 = 4 + 0j; 	
				std::complex<double> c6 = 6 + 0j; 	
				std::complex<double> c8 = 8 + 0j; 	
				std::complex<double> c9 = 9 + 0j; 	
				std::complex<double> c12 = 12 + 0j; 	
				std::complex<double> c16 = 16 + 0j; 	
				std::complex<double> c27 = 27 + 0j; 	
				std::complex<double> c36 = 36 + 0j; 	

				std::complex<double> sinPsi_c = sinPsi + 0j; 	
				std::complex<double> w_c = w + 0j; 	
				std::complex<double> Gamma_c = Gamma + 0j; 	
				
				x1 = -c1/c2*sinPsi_c*w_c - c1/c2*sqrt(c2*pow(sinPsi_c,c2)*pow(w_c,c2) - c4/c3*pow(w_c,c2) - c1/c9*(pow(w_c,c4) -
				c12*Gamma_c)/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c +
				c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) -
				c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) - c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)) +
				c6*(pow(sinPsi_c,c3)*pow(w_c,c3) - sinPsi_c*pow(w_c,c3))/sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - c1/c6*sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)));

				x2 = -c1/c2*sinPsi_c*w_c + c1/c2*sqrt(c2*pow(sinPsi_c,c2)*pow(w_c,c2) - c4/c3*pow(w_c,c2) - c1/c9*(pow(w_c,c4) -
				c12*Gamma_c)/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c +
				c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) -
				c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) - c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)) +
				c6*(pow(sinPsi_c,c3)*pow(w_c,c3) - sinPsi_c*pow(w_c,c3))/sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - c1/c6*sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)));


				x3 = -c1/c2*sinPsi_c*w_c - c1/c2*sqrt(c2*pow(sinPsi_c,c2)*pow(w_c,c2) - c4/c3*pow(w_c,c2) - c1/c9*(pow(w_c,c4) -
				c12*Gamma_c)/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c +
				c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) -
				c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) - c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)) -
				c6*(pow(sinPsi_c,c3)*pow(w_c,c3) - sinPsi_c*pow(w_c,c3))/sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) + c1/c6*sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)));
			
				x4 = -c1/c2*sinPsi_c*w_c + c1/c2*sqrt(c2*pow(sinPsi_c,c2)*pow(w_c,c2) - c4/c3*pow(w_c,c2) - c1/c9*(pow(w_c,c4) -
				c12*Gamma_c)/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c +
				c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) -
				c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) - c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)) -
				c6*(pow(sinPsi_c,c3)*pow(w_c,c3) - sinPsi_c*pow(w_c,c3))/sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) - pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))) + c1/c6*sqrt((pow(w_c,c4) + c3*pow((c1/c27*pow(w_c,c6) -
				c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3))*(c3*pow(sinPsi_c,c2) - c2)*pow(w_c,c2) - c12*Gamma_c + c9*pow((c1/c27*pow(w_c,c6)
				- c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) - c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) -
				c1)*pow(w_c,c8) - (c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c2/c3)))/pow((c1/c27*pow(w_c,c6) - c2/c3*(c3*pow(sinPsi_c,c2)*pow(w_c,c2) -
				c2*pow(w_c,c2))*Gamma_c + c2/c3*sqrt(c1/c3)*sqrt(-((pow(sinPsi_c,c2) - c1)*pow(w_c,c8) -
				(c27*Gamma_c*pow(sinPsi_c,c4) - c36*Gamma_c*pow(sinPsi_c,c2) + c8*Gamma_c)*pow(w_c,c4) -
				c16*pow(Gamma_c,c2))*Gamma_c)),(c1/c3)));



				// Choose the real positive root 
				if (x1.real() > 0 && fabs(x1.imag()) < 0.000001)
					{
						dBladeUindMag[i] = x1.real() ;
					}else if(x2.real() > 0 && fabs(x2.imag()) < 0.000001)
					{
						dBladeUindMag[i] = x2.real() ;
					}else if(x3.real() > 0 && fabs(x3.imag()) < 0.000001)
					{
						dBladeUindMag[i] = x3.real();
					}else if(x4.real() > 0 && fabs(x4.imag()) < 0.000001)
					{
						dBladeUindMag[i] = x4.real();
					}		
				#endif

				// Calculating the Induced velocity component Lower Half
				dBladeUind[i] = ::Zero3;
				if (dT > std::numeric_limits<doublereal>::epsilon())
				{
					// Option 1: Blade Azimuth location
					// dBladeUind[i](2) = dBladeUindMag[i]*dBladePos[i](2)/dBladePos[i].Norm() - w*NetForce(2)/NetForce.Norm();
					// dBladeUind[i](3) = dBladeUindMag[i]*dBladePos[i](3)/dBladePos[i].Norm() - w*NetForce(3)/NetForce.Norm();
					
					dBladeUind[i](2) = dBladeUindMag[i]*dBladePos[i](2)/dBladePos[i].Norm();
					dBladeUind[i](3) = dBladeUindMag[i]*dBladePos[i](3)/dBladePos[i].Norm();

					// Option 2: Local blade net force vector direction
					// dBladeUind[i](2) = dBladeUindMag[i]*dBladeForce[i](2)/dBladeForce[i].Norm() - w*NetForce(2)/NetForce.Norm();
					// dBladeUind[i](3) = dBladeUindMag[i]*dBladeForce[i](3)/dBladeForce[i].Norm() - w*NetForce(3)/NetForce.Norm();
					
					// Option 3: Rotor, Net thrust vector direction
					// dBladeUind[i](2) = fabs(dBladeUindMag[i]-w)*F(2)/dT;
					// dBladeUind[i](3) = fabs(dBladeUindMag[i]-w)*F(3)/dT;

					// Option 3: According to the Article BRNDICT ET AL   
					// dBladeUind[i](2) = dBladeUindMag[i] + w*dBladePos[i](3)/dBladePos[i].Norm();		
					// dBladeUind[i](3) = w*dBladePos[i](2)/dBladePos[i].Norm();		
				}

				dBladeUind[i](1) = (1 - dWeight)*dBladeUind[i](1) + dWeight*dBladeUindPrev[i](1);
				dBladeUind[i](2) = (1 - dWeight)*dBladeUind[i](2) + dWeight*dBladeUindPrev[i](2);
				dBladeUind[i](3) = (1 - dWeight)*dBladeUind[i](3) + dWeight*dBladeUindPrev[i](3);
				///////////////////////////////////////////
			}   
		}
	}
	
	ResetForce();
	WorkVec.Resize(0);
	return WorkVec;
}

void
CyclocopterDMST::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Calculates the azimuth position of the first blade */
	if (bFlagIsFirstBlade) 
	{
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}

	/* Calculates the moment only if output is required */
	if (bToBeOutput()) 
	{
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else 
	{
		Res.AddForce(F);
	}
}
int GlobalCount=0;
Vec3
CyclocopterDMST::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	//printf("%f %f %f\n",dUind(1),dUind(2),dUind(3));
	GlobalCount = GlobalCount + 1;
	Vec3 dUindSet;
	Vec3 dUindDMSTAvg = ::Zero3;
	for (int i = 0; i < NBlades; ++i)
	{
		if (uLabel == uBlade[i])
		{
			dUindSet = dBladeUind[i];
		}
		dUindDMSTAvg +=(dBladeUind[i]/NBlades);
	}
	std::cout << "dBlade 1 position " << dBladePos[0] << "  angle   " << atan2(dBladePos[0](2),dBladePos[0](3))*180/M_PI<< std::endl;
	// std::cout << uLabel <<"\t" << GlobalCount <<"\t" << "2D" << "\t" << dUind << "\t" <<dUind.Norm() << "\t" << "DMST" << "\t" << dUindSet << "\t" << dUindSet.Norm() << std::endl;
	return RRotor*(dUindSet); // DMST model
	// return RRotor*(dUind); // 2D model
}

/* CyclocopterDMST - end */



#if 0 // KARI NOT IMPLEMENTED YET!

/* CyclocopterKARI - begin */

/*

From:

"A New VTOL UAV Cyclocopter with Cycloidal Blades System",

Chul Yong Yun, Illkyung Park,
Ho Yong Lee, Jai Sang Jung, In Seong Hwang,
Seung Jo Kim,
Sung Nam Jung

Presented at the American Helicopter Society 60th Annual Forum,
Baltimore, MD, June 7-10, 2004.
*/

class CyclocopterKARI
: virtual public Elem, public CyclocopterInflow {
public:
	CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterKARI(void);

	// Elaborate internal state after convergence
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; it is assumed that each element knows,
	// by the OutputHandler, where to write its own output
	virtual void Output(OutputHandler& OH) const;

	// Contribution to the restart file
	virtual std::ostream& Restart(std::ostream& out) const;

	// residual assembly
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Adds to the forces the contribution from an element
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);
#endif

	// Restores the induced velocity to an element
	// based on the azimuth position
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// *******FOR PARALLEL SOLVER********
	// Provides the type and label of the nodes connected to the element
	// useful for the assembly of the DOF connexion matrix
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

CyclocopterKARI::CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, 0),
CyclocopterInflow(uL, pDO),
pRotor(0),
RRot(::Eye3)
{
	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));
}

CyclocopterKARI::~CyclocopterKARI(void)
{
	NO_OP;
}

void
CyclocopterKARI::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterKARI::Output(OutputHandler& OH) const
{
	NO_OP;
}

std::ostream&
CyclocopterKARI::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

SubVectorHandler&
CyclocopterKARI::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	return WorkVec;
}

#if 0
void
CyclocopterKARI::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	NO_OP;
}
#endif

Vec3
CyclocopterKARI::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return Zero3;
}

void
CyclocopterKARI::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

/* CyclocopterKARI - end */

#endif // KARI NOT IMPLEMENTED YET!

bool
mbdyn_cyclocopter_set(void)
{
	UserDefinedElemRead *rf;
	
	rf = new UDERead<CyclocopterNoInflow>;
	if (!SetUDE("cyclocopter" "no" "inflow", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter no inflow\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterUniform1D>;
	if (!SetUDE("cyclocopter" "uniform" "1D", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter uniform 1D\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterUniform2D>;
	if (!SetUDE("cyclocopter" "uniform" "2D", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter uniform 2D\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterPolimi>;
	if (!SetUDE("cyclocopter" "Polimi", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter Polimi\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterDMST>;
	if (!SetUDE("cyclocopter" "DMST", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter DMST\""
			<< std::endl);

		return false;
	}

#if 0
	rf = new UDERead<CyclocopterKARI>;
	if (!SetUDE("cyclocopter" "KARI", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter KARI\""
			<< std::endl);

		return false;
	}
#endif

	return true;
}

#ifdef MBDYN_MODULE

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!mbdyn_cyclocopter_set()) {
		silent_cerr("cyclocopter: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // MBDYN_MODULE
