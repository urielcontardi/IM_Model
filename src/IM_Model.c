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
    double is_alpha;
    double is_beta;
    double dfluxR_alpha; // derivative
    double dfluxR_beta; // derivative
    double fluxR_alpha;
    double fluxR_beta;
    double Te;
    double w;
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
    double w = states->w;

    // Inputs (u_alpha, u_beta)
    double u_alpha = intInputs->valpha;
    double u_beta = intInputs->vbeta;
    double Tload = self->inp.Tload;

    // Update Currents
    states->is_alpha = is_alpha + (Lr/(Lm*Lm-Ls))*(Rs*is_alpha+(Lm*dfluxR_alpha/Lr)-u_alpha) * Ts;
    states->is_beta = is_beta + (Lr/(Lm*Lm-Ls))*(Rs*is_beta+(Lm*dfluxR_alpha/Lr)-u_beta) * Ts;

    // Update Flux
    states->dfluxR_alpha = (states->is_alpha*Rr*Lm/Lr)-(npp*w*fluxR_beta)-(Rr*Lm*fluxR_alpha/Lr);
    states->dfluxR_beta = (states->is_beta*Rr*Lm/Lr)+(npp*w*dfluxR_beta)-(Rr*Lm*fluxR_beta/Lr);
    states->fluxR_alpha = states->fluxR_alpha + states->dfluxR_alpha * Ts;
    states->fluxR_beta = states->fluxR_beta + states->dfluxR_beta * Ts;

    // Update Torque
    states->Te = ((3*npp*Lm)/(2*Lr))*(states->fluxR_alpha*states->is_beta-states->fluxR_beta*states->is_alpha);
    
    // Update Speed
    states->w = states->w + ((states->Te-Tload)/J) * Ts;
}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *states = &privateData->states;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // Inverse Clarke Transform (αβ -> abc)
    self->out.ia = states->is_alpha + intInputs->v0;
    self->out.ib = -0.5 * states->is_alpha + (sqrt(3.0) / 2.0) * states->is_beta + intInputs->v0;
    self->out.ic = -0.5 * states->is_alpha - (sqrt(3.0) / 2.0) * states->is_beta + intInputs->v0;
    self->out.wr = states->w;

}