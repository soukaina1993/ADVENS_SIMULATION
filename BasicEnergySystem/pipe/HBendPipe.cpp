//
// Created by lucile.schulthe on 24.11.2022.
//

#include "HBendPipe.h"

using namespace std;

HBendPipe::HBendPipe(const char *iPipeName, TInputFile &InputFile){
    cerr << "Not implemented" << endl;
    exit(1);
}

HBendPipe::HBendPipe(HPipe* pipe){
    flow=pipe->flow;
    PipeNumber=pipe->PipeNumber;
    Length=pipe->Length;
    Diam=pipe->Diam;
    RugAbs=pipe->RugAbs;
    U_pipe=pipe->U_pipe;
    U_total=pipe->U_total;
    tau=pipe->tau;
    DeltaPr=pipe->DeltaPr;
    DeltaT=pipe->DeltaT;
    Re=pipe->Re;
    Speed=pipe->Speed;
    F_Moody=pipe->F_Moody;
    pStateIn=pipe->pStateIn;
}

HBendPipe::HBendPipe(HPipe& pipe){
    flow=pipe.flow;
    PipeNumber=pipe.PipeNumber;
    Length=pipe.Length;
    Diam=pipe.Diam;
    RugAbs=pipe.RugAbs;
    U_pipe=pipe.U_pipe;
    U_total=pipe.U_total;
    tau=pipe.tau;
    DeltaPr=pipe.DeltaPr;
    DeltaT=pipe.DeltaT;
    Re=pipe.Re;
    Speed=pipe.Speed;
    F_Moody=pipe.F_Moody;
    pStateIn=pipe.pStateIn;
}

HBendPipe::HBendPipe() {
    flow=new Flow;
    pStateIn=new TPhysicState;
    flow_created = state_created = true;
}

HBendPipe::HBendPipe(Flow *iflow) {
    flow=iflow;
    pStateIn=new TPhysicState;          // allocate new physical state
    state_created = true;
}

HBendPipe::HBendPipe(Flow_physic *iflow) {
    flow=iflow;
    pStateIn=iflow->pStateIn;           // pointer to physical state
}

HBendPipe::~HBendPipe() {
    if (flow_created) {
        delete flow;
        flow_created=false;
    }
    if (state_created) {
        delete pStateIn;
        state_created=false;
    }
}



void HBendPipe::PressureDropCalculation() {
    double zeta=0.1; // TODO: depends from angle, bendradius, Re, roughness,... find general formula
/*
 https://tecciness.de/hilfe/bogendruckverlust.php       (if angle=90)
 https://doc.modelica.org/om/Modelica.Fluid.Dissipation.Utilities.SharedDocumentation.PressureLoss.Bend.dp_curvedOverall.html
 */
    if (Speed==0.0)
        Speed = *flow->Massflow / pStateIn->Ro / (M_PI * sqr(Diam) / 4);
    DeltaPr = zeta * pStateIn->Ro * sqr(Speed)/2.0;
}

double HBendPipe::Heattransfer_conv(){
    // TODO
}