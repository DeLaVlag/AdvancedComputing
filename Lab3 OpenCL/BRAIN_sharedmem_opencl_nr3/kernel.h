#include "variables.h"

#ifdef DOUBLE_PRECISION
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

void CompDend(local mod_prec *cellCompParamsPtr);///
void DendHCurr(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_newComp1);///
void DendCaCurr(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_newComp1);///
void DendKCurr(local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_prevComp2, local mod_prec *chPrms_newComp1);///	
void DendCal(local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_prevComp2, local mod_prec *chPrms_newComp1);///
void DendCurrVolt(mod_prec chComps_iC, local mod_prec *chComps_iApp, local mod_prec *chComps_vDend, local mod_prec *chComps_newVDend, 
    local mod_prec *chComps_vSoma, local mod_prec *chComps_q, local mod_prec *chComps_r, local mod_prec *chComps_s, 
    local mod_prec *chComps_newI_CaH);///
mod_prec IcNeighbors(local mod_prec *neighVdend, mod_prec prevV_dend);///
void CompSoma(local mod_prec *cellCompParamsPtr);///
void SomaCalcium(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_prevComp2, local mod_prec *chPrms_newComp1, 
    local mod_prec *chPrms_newComp2);///
void SomaSodium(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_prevComp2, local mod_prec *chPrms_newComp1, 
    local mod_prec *chPrms_newComp2);///
void SomaPotassium(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_prevComp2, 
    local mod_prec *chPrms_newComp1, local mod_prec *chPrms_newComp2);///
void SomaPotassiumX(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_newComp1);///
void SomaCurrVolt(local mod_prec *chComps_g_CaL, local mod_prec *chComps_vDend, local mod_prec *chComps_vSoma, 
    local mod_prec *chComps_newVSoma, local mod_prec *chComps_vAxon, local mod_prec *chComps_k, local mod_prec *chComps_l, 
    local mod_prec *chComps_m, local mod_prec *chComps_h, local mod_prec *chComps_n, local mod_prec *chComps_x_s);///
void CompAxon(local mod_prec *cellCompParamsPtr);///
void AxonSodium(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_newComp1, local mod_prec *chPrms_newComp2);///
void AxonPotassium(local mod_prec *chPrms_v, local mod_prec *chPrms_prevComp1, local mod_prec *chPrms_newComp1);///
void AxonCurrVolt(local mod_prec *chComps_vSoma, local mod_prec *chComps_vAxon, local mod_prec *chComps_newVAxon, local mod_prec *chComps_m_a, 
    local mod_prec *chComps_h_a, local mod_prec *chComps_x_a);///