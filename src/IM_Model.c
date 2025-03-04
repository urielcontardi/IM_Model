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
    double fluxR_alpha;
    double fluxR_beta;
    double wr;
    double wm;
    double Te;
} IM_States_t;

typedef struct {
    double valpha;
    double vbeta;
    double v0;
} IM_InternalInputs_t;

typedef struct {
    IM_InternalInputs_t inp;
    IM_States_t out;
} IM_PrivateData_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        LOCAL FUNCTIONS PROTOTYPES                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static void _updateInternalParameters(IM_Model_t *self);
static void _vabc2AphaBeta(IM_Model_t *self);
static void _updateOutputs(IM_Model_t *self);

static void _updateStatesA(IM_Model_t *self);
static void _updateStatesB(IM_Model_t *self);
static void _updateStatesC(IM_Model_t *self);

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

void IM_TypeModel(IM_Model_t *self, IMType model)
{
    self->type = model;
    memset(self->priv, 0, sizeof(IM_PrivateData_t));
}

void IM_SimulateStep(IM_Model_t *self) {

    // Compute selected model
    switch (self->type )
    {
        default:
        case MODEL_A:
            _updateStatesA(self);
            break;

        case MODEL_B:
            _updateStatesB(self);
            break;

        case MODEL_C:
            _updateStatesC(self);
            break;
    
    }

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

void _updateStatesA(IM_Model_t *self) {
    /// Article: FPGA based real-time model of three-phase induction machine
    /// https://ieeexplore.ieee.org/document/8395534
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls;
    double Lr  = self->params.Lr;
    double J   = self->params.J; 
    double npp = self->params.npp;
    // double poles = 2*npp;
    double Ts  = self->params.Ts;

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Article Equations :
    double dfluxR_alpha = (Rr*Lm/Lr)*out->is_alpha - out->wr * out->fluxR_beta - (Rr/Lr) * out->fluxR_alpha;
    double dfluxR_beta = (Rr*Lm/Lr)*out->is_beta + out->wr * out->fluxR_alpha - (Rr/Lr) * out->fluxR_beta;
    double dis_alpha = (Lr/(Lm*Lm - Lr*Ls))*(Rs * out->is_alpha + (Lm/Lr)*dfluxR_alpha - v_alpha_s);
    double dis_beta = (Lr/(Lm*Lm - Lr*Ls))*(Rs * out->is_beta + (Lm/Lr)*dfluxR_beta - v_beta_s);

    double Te = ((3*npp*Lm)/(2*Lr))*(out->fluxR_alpha*out->is_beta - out->fluxR_beta*out->is_alpha);
    double dwm = (Te - Tload)/J;

    // Update States
    out->is_alpha    = out->is_alpha + dis_alpha * Ts;
    out->is_beta     = out->is_beta + dis_beta * Ts;
    out->fluxR_alpha = out->fluxR_alpha + dfluxR_alpha * Ts;
    out->fluxR_beta  = out->fluxR_beta + dfluxR_beta * Ts;
    out->wm          = out->wm + dwm * Ts;
    out->wr          = out->wm * npp;
    out->Te          = Te;

}

void _updateStatesB(IM_Model_t *self) {
    /// Article: C++ based dynamic model of AC induction motor in discrete time domain
    /// https://ieeexplore.ieee.org/document/7512961
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
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

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double v_alpha_r = 0.0f; // Induction Motor
    double v_beta_r  = 0.0f;
    double Tload     = self->inp.Tload;
    double sigma = 1.0 / (Lm*Lm - Lr*Ls);
    double wr = out->wr;

    // A_d / B_d Matrix
    double Ad[4][4] = {
        {1.0 + sigma * Ts * (Lr * Rs),  -sigma * Ts * (Lm * Lm * wr), -sigma * Ts * (Lm * Rr), -sigma * Ts * (Lm * Lr * wr)},
        {sigma * Ts * (Lm * Lm * wr),    1.0 + sigma * Ts * (Lr * Rs), sigma * Ts * (Lm * Lr * wr), -sigma * Ts * (Lm * Rr)},
        {-sigma * Ts * (Lm * Rs),       sigma * Ts * (Lm * Ls * wr), 1.0 + sigma * Ts * (Ls * Rr), sigma * Ts * (Lr * Ls * wr)},
        {-sigma * Ts * (Lm * Ls * wr), -sigma * Ts * (Lm * Rs), -sigma * Ts * (Lr * Ls * wr), 1.0 + sigma * Ts * (Ls * Rr)}
    };

    double Bd[4][4] = {
        {-sigma * Lr * Ts,  0,                sigma * Lm * Ts,   0},
        {0,                 -sigma * Lr * Ts, 0,                 sigma * Lm * Ts},
        {sigma * Lm * Ts,   0,                -sigma * Ls * Ts,  0},
        {0,                 sigma * Lm * Ts,  0,                 -sigma * Ls * Ts}
    };

    // States
    double Ik[4] = {out->is_alpha, out->is_beta, out->ir_alpha, out->ir_beta};
    double Vk[4] = {v_alpha_s, v_beta_s, v_alpha_r, v_beta_r};
    double Iknext[4] = {0};

    // Update States
    for (int i = 0; i < 4; i++) {
        Iknext[i] = 0;
        for (int j = 0; j < 4; j++) {
            Iknext[i] += Ad[i][j] * Ik[j] + Bd[i][j] * Vk[j];
        }
    }

    // Update Torque / Speed
    double Te = (npp * Lm / 3.0) * (Iknext[2] * Iknext[1] - Iknext[3] * Iknext[0]);
    double dwm = (Te - Tload)/J;

    // Update Output
    out->is_alpha = Iknext[0];
    out->is_beta  = Iknext[1];
    out->ir_alpha = Iknext[2];
    out->ir_beta  = Iknext[3];
    out->wm       = out->wm + Ts * dwm;
    out->wr       = out->wm * npp;
    out->Te       = Te;
}

void _updateStatesC(IM_Model_t *self) {
    /// Article: Real-Time Emulator of an Induction Motor: FPGA-based Implementation
    /// https://ieeexplore.ieee.org/document/6421152/
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
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

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Constanst
    double Tr = Lr/Rr;
    double sigma = 1.0 - ((Lm*Lm)/(Ls*Lr));
    double K = Lm/(sigma * Ls * Lr);
    double gama = (Rs/(sigma * Ls)) + ((Rr*Lm*Lm)/(sigma*Ls*Lr*Lr));
    
    // Derivative
    double dis_alpha = 0.0f;
    double dis_beta = 0.0f;
    double dfluxR_alpha = 0.0f;
    double dfluxR_beta = 0.0f;
    double dwm = 0.0f;

    dis_alpha = -gama*out->is_alpha + (out->fluxR_alpha*K/Tr) + npp*out->wm*K*out->fluxR_beta + (v_alpha_s/(sigma*Ls));
    dis_beta = -gama*out->is_beta + (out->fluxR_beta*K/Tr) - npp*out->wm*K*out->fluxR_alpha + (v_beta_s/(sigma*Ls));
    dfluxR_alpha = (out->is_alpha*Lm/Tr) - out->fluxR_alpha/Tr - npp*out->wm*out->fluxR_beta;
    dfluxR_beta = (out->is_beta*Lm/Tr) - out->fluxR_beta/Tr + npp*out->wm*out->fluxR_alpha;
    dwm = (npp*Lm/(J*Lr))*(out->fluxR_alpha*out->is_beta - out->fluxR_beta*out->is_alpha) - (Tload/J);

    // Compute next value
    out->is_alpha = out->is_alpha + dis_alpha * Ts;
    out->is_beta  = out->is_beta + dis_beta * Ts;
    out->fluxR_alpha = out->fluxR_alpha + dfluxR_alpha * Ts;
    out->fluxR_beta  = out->fluxR_beta + dfluxR_beta * Ts;
    out->wm = out->wm + dwm * Ts;

}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // Inverse Clarke Transform (αβ -> abc)
    self->out.ia = out->is_alpha + intInputs->v0;
    self->out.ib = -0.5 * out->is_alpha + (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.ic = -0.5 * out->is_alpha - (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.wmec = out->wm;
    self->out.wr = out->wr;
    self->out.Te = out->Te;

}
