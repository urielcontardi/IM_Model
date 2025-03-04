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
    double wr;
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
    intInputs->v0 = 0.0f; // (1.0/3.0) * (Va+Vb+Vc);

}

void _updateStates(IM_Model_t *self) {
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
        {sigma * Ts * (-Lr*Rs) + 1.0,     sigma * Ts * (-Lm*Lm*wr),    sigma * Ts * (-Lm*Rr),     sigma * Ts * (-Lm*Lr*wr)    },
        {sigma * Ts * (Lm*Lm*wr),   sigma * Ts * (Lr*Rs) + 1.0,        sigma * Ts * (Lm*Lr*wr),   sigma * Ts * (-Lm*Rr)       },
        {sigma * Ts * (-Lm*Rs),     sigma * Ts * (Lm*Ls*wr),     sigma * Ts * (Ls*Rr) + 1.0,      sigma * Ts * (Lr*Ls*wr)     },
        {sigma * Ts * (-Lm*Ls*wr),  sigma * Ts * (-Lm*Rs),       sigma * Ts * (-Lr*Ls*wr),  sigma * Ts * (Ls*Rr) + 1.0       }
    };
    
    double Bd[4][4] = {
        {-Lr * Ts * sigma,      0,                          Lm * Ts * sigma,            0                   },
        {0,                     -Lr * Ts * sigma,           0,                          Lm * Ts * sigma     },
        {Lm * Ts * sigma,       0,                          -Ls * Ts * sigma ,          0                   },
        {0,                     Lm * Ts * sigma,            0,                          -Ls * Ts * sigma    }
    };

    // States
    double Ik[4] = {out->is_alpha, out->is_beta, out->ir_alpha, out->ir_beta};
    double Vk[4] = {v_alpha_s, v_beta_s, v_alpha_r, v_beta_r};
    double Iknext[4] = {0};

    // Update States
    for (int i = 0; i < 4; i++) {
        Iknext[i] = 0;
        for (int j = 0; j < 4; j++) {
            Iknext[i] += Ad[i][j] * Ik[j];
            Iknext[i] += Bd[i][j] * Vk[j];
        }
    }

    // Update Torque
    double Te = (npp * Lm / 3.0) * (Iknext[2] * Iknext[1] - Iknext[3] * Iknext[0]);
    
    // Update Speed
    double wr_next = out->wr + Ts * (npp / (2.0 * J)) * (Te - Tload);

    // Update Output
    out->is_alpha = Iknext[0];
    out->is_beta  = Iknext[1];
    out->ir_alpha = Iknext[2];
    out->ir_beta  = Iknext[3];
    out->wr       = wr_next;
}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // Inverse Clarke Transform (αβ -> abc)
    self->out.ia = out->is_alpha + intInputs->v0;
    self->out.ib = -0.5 * out->is_alpha + (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.ic = -0.5 * out->is_alpha - (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    self->out.wmec = out->wr * self->params.npp;
    self->out.wr = out->wr;

}
