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
    double ir_alpha;
    double ir_beta;
    double wmec;
} IM_States_t;

typedef struct {
    double valpha;
    double vbeta;
    double v0;
} IM_InternalInputs_t;

typedef struct {
    IM_InternalInputs_t inp;
    IM_States_t states;
    IM_States_t out;
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

static IM_PrivateData_t _privateData = {0};

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
    intInputs->vbeta =  (sqrt(3.0) / 3.0) * (Vb - Vc);
    intInputs->v0 = (1.0/3.0) * (Va+Vb+Vc);

}

void _updateStates(IM_Model_t *self) {
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *states = &privateData->states;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = (3.0/2.0)*self->params.Lm;
    double Ls  = (3.0/2.0)*self->params.Ls;
    double Lr  = (3.0/2.0)*self->params.Lr;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inputs (u_alpha, u_beta)
    double u_alpha = intInputs->valpha;
    double u_beta = intInputs->vbeta;
    double Tload = self->inp.Tload;

    double sigma = (double) 1.0f/(Ls*Lr-Lm*Lm);

    // Update States
    states->is_alpha = sigma * (
                                (Lr * (u_alpha - Rs*out->is_alpha)) +
                                (-Lm * (-npp*out->wmec*(Lm*out->is_beta + Lr * out->ir_beta) - Rr * out->ir_alpha))
                                );

    states->ir_alpha = sigma * (
                               (-Lm * (u_alpha - Rs*out->is_alpha)) + 
                               (Ls * (-npp*out->wmec*(Lm*out->is_beta + Lr * out->ir_beta) - Rr * out->ir_alpha))
                               );

    states->is_beta = sigma * (
                                (Lr * (u_beta - Rs*out->is_beta)) +
                                (-Lm * (npp*out->wmec*(Lm*out->is_alpha + Lr * out->ir_alpha) - Rr * out->ir_beta))
                                );

    states->ir_beta = sigma * (
                            (-Lm * (u_beta - Rs*out->is_beta)) +
                            (Ls * (npp*out->wmec*(Lm*out->is_alpha + Lr * out->ir_alpha) - Rr * out->ir_beta))
                            );

    // Torque
    double Te = (3.0/2.0) * (npp * Lm / Lr) * 
                ((Lm * out->is_alpha + Lr * out->ir_alpha) * out->is_beta -
                (Lm * out->is_beta  + Lr * out->ir_beta)  * out->is_alpha);
  
    states->wmec = (Te - Tload) / J;

    // Euler Discretization
    out->is_alpha    = out->is_alpha + states->is_alpha * Ts;
    out->is_beta     = out->is_beta + states->is_beta * Ts;
    out->ir_alpha    = out->ir_alpha + states->ir_alpha * Ts;
    out->ir_beta     = out->ir_beta + states->ir_beta * Ts;
    out->wmec        = out->wmec + states->wmec * Ts;

}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // Inverse Clarke Transform (αβ -> abc)
    self->out.ia = out->is_alpha + intInputs->v0;
    self->out.ib = -0.5 * out->is_alpha + (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.ic = -0.5 * out->is_alpha - (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.wmec = out->wmec;
    self->out.wr = out->wmec * self->params.npp;

}
