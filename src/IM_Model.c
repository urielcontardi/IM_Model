/// \file		IM_Model.c
///
/// \brief	
///
/// \author		Uriel Abe Contardi (urielcontardi@hotmail.com)
/// \date		26-01-2025
///
/// \version	1.0
///
/// \note		Revisions:
/// 			26-01-2025 <urielcontardi@hotmail.com>
/// 			First revision.
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INCLUDES                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "IM_Model.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           DEFINES AND MACROS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      LOCAL TYPEDEFS AND STRUCTURES                       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
typedef struct {
    double dis_alpha;    // derivative
    double dis_beta;     // derivative
    double dfluxR_alpha; // derivative
    double dfluxR_beta; // derivative
    double dwmec;       // derivative
    double is_alpha;
    double is_beta;
    double fluxR_alpha;
    double fluxR_beta;
    double wmec;
    double Te;
} IM_States_t;

typedef struct {
    double valpha;
    double vbeta;
    double v0;
} IM_InternalInputs_t;

typedef struct {
    IM_InternalInputs_t inp;
    IM_States_t states;
} IM_PrivateData_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        LOCAL FUNCTIONS PROTOTYPES                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static void _updateInternalParameters(IM_Model_t *self);
static void _vabc2AphaBeta(IM_Model_t *self);
static void _updateStates(IM_Model_t *self);
static void _updateOutputs(IM_Model_t *self);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      STATIC VARIABLES AND CONSTANTS                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

static IM_PrivateData_t _privateData;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void IM_Init(IM_Model_t *self) {
    if (self == NULL) {
        return; 
    }

    // Clear the structure
    memset(self, 0, sizeof(IM_Model_t));

    // Allocation of the for the priv data
    self->priv = (void*)&_privateData;
    if (self->priv == NULL) {
        return;
    }

    // Clear priv data
    memset(self->priv, 0, sizeof(IM_PrivateData_t));
}

void IM_SetParams(IM_Model_t *self, const IMParams *params) {
    if (self == NULL || params == NULL) {
        return;
    }

    memcpy(&self->params, params, sizeof(IMParams));
    _updateInternalParameters(self);
}

void IM_SetInputs(IM_Model_t *self, const IMInputs *inputs) {
    if (self == NULL || inputs == NULL) {
        return;
    }

    memcpy(&self->inp, inputs, sizeof(IMInputs));
    _vabc2AphaBeta(self);
}

void IM_SimulateStep(IM_Model_t *self) {
    _updateStates(self); 
    _updateOutputs(self);
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              LOCAL FUNCTIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void _updateInternalParameters(IM_Model_t *self) {

}

void _vabc2AphaBeta(IM_Model_t *self) {

    // Get the input values
    double Va = self->inp.Va;
    double Vb = self->inp.Vb;
    double Vc = self->inp.Vc;
    
    // Clark Transform (abc -> αβ)
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_InternalInputs_t *intInputs = &privateData->inp;
    intInputs->valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
    intInputs->vbeta =  (1.0 / sqrt(3.0)) * (Vb - Vc);
    intInputs->v0 = (1.0/3.0) * (Va+Vb+Vc);

}

void _updateStates(IM_Model_t *self) {
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *states = &privateData->states;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls;
    double Lr  = self->params.Lr;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Last States
    double is_alpha = states->is_alpha;
    double is_beta = states->is_beta;
    double dfluxR_alpha = states->dfluxR_alpha;
    double dfluxR_beta = states->dfluxR_beta;
    double fluxR_alpha = states->fluxR_alpha;
    double fluxR_beta = states->fluxR_beta;
    double wmec = states->wmec;
    double w = wmec * npp;
    double Te = states->Te;

    // Inputs (u_alpha, u_beta)
    double u_alpha = intInputs->valpha;
    double u_beta = intInputs->vbeta;
    double Tload = self->inp.Tload;

    // Derivative: Currents / Flux / Speed
    states->dis_alpha = (1/Ls)*(u_alpha - Rs*is_alpha - Lm*dfluxR_alpha);
    states->dis_beta = (1/Ls)*(u_beta - Rs*is_beta - Lm*dfluxR_beta);
    states->dfluxR_alpha = (Rr*Lm/Lr)*is_alpha + w*fluxR_beta - (Rr/Lr)*fluxR_alpha;
    states->dfluxR_beta = (Rr*Lm/Lr)*is_beta + w*fluxR_alpha - (Rr/Lr)*fluxR_beta;
    states->dwmec = (Te-Tload)/J;

    // Update Torque
    states->Te = ((3*npp*Lm)/(2*Lr))*(fluxR_alpha*is_beta - fluxR_beta*is_alpha);

    // Update States
    states->is_alpha    = states->is_alpha + states->dis_alpha * Ts;
    states->is_beta     = states->is_beta + states->dis_beta * Ts;
    states->fluxR_alpha = states->fluxR_alpha + states->dfluxR_alpha * Ts;
    states->fluxR_beta  = states->fluxR_beta + states->dfluxR_beta * Ts;
    states->wmec        = states->wmec + states->dwmec * Ts;
}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *states = &privateData->states;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // Inverse Clarke Transform (αβ -> abc)
    self->out.ia = states->is_alpha + intInputs->v0;
    self->out.ib = -0.5 * states->is_alpha + (sqrt(3.0) / 2.0) * states->is_beta + intInputs->v0;
    self->out.ic = -0.5 * states->is_alpha - (sqrt(3.0) / 2.0) * states->is_beta + intInputs->v0;
    self->out.wmec = states->wmec;
    self->out.wr = states->wmec * self->params.npp;

}
