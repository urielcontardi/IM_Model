/// \file		psim.c
///
/// \brief	
///
/// \author		Uriel Abe Contardi (urielcontardi@hotmail.com)
/// \date		27-02-2025
///
/// \version	1.0
///
/// \note		Revisions:
/// 			27-02-2025 <urielcontardi@hotmail.com>
/// 			First revision.
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INCLUDES                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <windows.h>
#include <IM_Model.h>

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
typedef enum
{
	VA=0,
    VB,
    VC,
    TLOAD,
    RS,
    RR,
    LM,
    LS,
    LR,
    J,
    NPP,
    TS,
	TOTAL_INPUTS
} Inputs_t;

typedef enum
{
    IA=0,
    IB,
    IC,
    WR,
	TOTAL_OUTPUTS
} Outputs_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        LOCAL FUNCTIONS PROTOTYPES                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static IM_Model_t* _getIMInstance(void);
static void _setParameters(IM_Model_t* model, double *simInputs);
static void _setInputs(IM_Model_t* model, double *simInputs);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      STATIC VARIABLES AND CONSTANTS                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) void simuser(double t, double delt, double *simInputs, double *simOut);

__declspec(dllexport) void simuser(t, delt, simInputs, simOut)
	double t, delt;
	double *simInputs, *simOut;
{

    // Induction Machine Model
    IM_Model_t* pIM = _getIMInstance();

    // Only runs the _setup when it's the first simulation cycle
	if (t <= delt){
        IM_Init(pIM);
        _setParameters(pIM, simInputs);
    }

    _setInputs(pIM, simInputs);
    IM_SimulateStep(pIM);

    simOut[IA] = (double)pIM->out.ia;
    simOut[IB] = (double)pIM->out.ib;
    simOut[IC] = (double)pIM->out.ic;
    simOut[WR] = (double)pIM->out.wr;

}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              LOCAL FUNCTIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static IM_Model_t* _getIMInstance(void)
{
	static IM_Model_t model;
	return &model;
}

static void _setParameters(IM_Model_t* model, double *simInputs)
{
    IMParams params;
    params.Rs = simInputs[RS];
    params.Rr = simInputs[RR];
    params.Lm = simInputs[LM];
    params.Ls = simInputs[LS];
    params.Lr = simInputs[LR];
    params.J = simInputs[J];
    params.npp = simInputs[NPP];
    params.Ts = simInputs[TS];
    IM_SetParams(model, &params);
}

static void _setInputs(IM_Model_t* model, double *simInputs)
{
    IMInputs inp;
    inp.Va = simInputs[VA];
    inp.Vb = simInputs[VB];
    inp.Vc = simInputs[VC];
    inp.Tload = simInputs[TLOAD];
    IM_SetInputs(model, &inp);
}