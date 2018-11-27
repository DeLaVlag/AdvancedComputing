#include "variables.h"

#ifdef DOUBLE_PRECISION
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

void CompDend(global mod_prec *cellCompParamsPtr, int k, int j);///
void DendHCurr(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_newComp1);///
void DendCaCurr(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_newComp1);///
void DendKCurr(global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_prevComp2, global mod_prec *chPrms_newComp1);///	
void DendCal(global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_prevComp2, global mod_prec *chPrms_newComp1);///
void DendCurrVolt(mod_prec chComps_iC, global mod_prec *chComps_iApp, global mod_prec *chComps_vDend, global mod_prec *chComps_newVDend, 
    global mod_prec *chComps_vSoma, global mod_prec *chComps_q, global mod_prec *chComps_r, global mod_prec *chComps_s, 
    global mod_prec *chComps_newI_CaH);///
mod_prec IcNeighbors(global mod_prec *neighVdend, mod_prec prevV_dend);///
void CompSoma(global mod_prec *cellCompParamsPtr, int k, int j);///
void SomaCalcium(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_prevComp2, global mod_prec *chPrms_newComp1, 
    global mod_prec *chPrms_newComp2);///
void SomaSodium(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_prevComp2, global mod_prec *chPrms_newComp1, 
    global mod_prec *chPrms_newComp2);///
void SomaPotassium(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_prevComp2, 
    global mod_prec *chPrms_newComp1, global mod_prec *chPrms_newComp2);///
void SomaPotassiumX(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_newComp1);///
void SomaCurrVolt(global mod_prec *chComps_g_CaL, global mod_prec *chComps_vDend, global mod_prec *chComps_vSoma, 
    global mod_prec *chComps_newVSoma, global mod_prec *chComps_vAxon, global mod_prec *chComps_k, global mod_prec *chComps_l, 
    global mod_prec *chComps_m, global mod_prec *chComps_h, global mod_prec *chComps_n, global mod_prec *chComps_x_s);///
void CompAxon(global mod_prec *cellCompParamsPtr, int k, int j);///
void AxonSodium(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_newComp1, global mod_prec *chPrms_newComp2);///
void AxonPotassium(global mod_prec *chPrms_v, global mod_prec *chPrms_prevComp1, global mod_prec *chPrms_newComp1);///
void AxonCurrVolt(global mod_prec *chComps_vSoma, global mod_prec *chComps_vAxon, global mod_prec *chComps_newVAxon, global mod_prec *chComps_m_a, 
    global mod_prec *chComps_h_a, global mod_prec *chComps_x_a);///