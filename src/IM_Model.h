/// \file		IM_Model.h
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

#ifndef IM_MODEL_H
#define IM_MODEL_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         TYPEDEFS AND STRUCTURES                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

typedef struct {
    double Rs;        // Resistência estatórica
    double Rr;        // Resistência rotórica
    double Lm;        // Indutância mutua
    double Ls;        // Indutância estatórica
    double Lr;        // Indutância rotórica
    double J;         // Momento de inércia
    double npp;       // Número de pares de polos
    double Ts;        // Tempo de amostragem
} IMParams;

typedef struct {
    double Va;
    double Vb;
    double Vc;
    double Tload;
} IMInputs;

typedef struct {
    double ia;
    double ib;
    double ic;
    double wr;
    double wmec;
} IMOutputs;

typedef struct {
    IMParams params;
    IMInputs inp;
    IMOutputs out;
    void *priv;
} IM_Model_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
void IM_Init(IM_Model_t *self);
void IM_SetParams(IM_Model_t *self, const IMParams *params);
void IM_SetInputs(IM_Model_t *self, const IMInputs *inputs);
void IM_SimulateStep(IM_Model_t *self);

#endif // IM_MODEL_H