#include "simulation.h"

Simulation::Simulation()
{
//    s =new CellSolver();
//    s.solve(100, 100);
//    s.solve(300, none);
   //simulate_first_porous_system_fromFile();
    try_config();
}


void Simulation::prepareSystem(){

}

void Simulation::prepareSystem(string filename){

}

void Simulation::thermalize(double T){

}

void Simulation::makePores(int num_pores, double r_min, double r_max){

}

void Simulation::solve(int timesteps){

}






void Simulation::simulate_first_porous_system_fromFile(){
    string filename = "thermalized_T_1_05_nano_20_pores.xyz";
    string path = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/";
    int num_time = 300;
    MakeMatrix* nano_porous_matrix = new MakeMatrix();
    CellSolver* solver = nano_porous_matrix->makeFromFile(path+filename, 20, 5.72, num_time);
    solver->thermostatON = false;
    double dt = 0.01;
    string outfile = "after_thermalized_porous";
    bool writeVMD = true;
    bool writeMeasurements = false;
    bool berendsen = false, andersen = false;
    solver->solve(0, num_time, dt, outfile,writeVMD, writeMeasurements,0, berendsen, andersen);

}

void Simulation::simulate_first_porous_system(){
    int max_time = 300;
    MakeMatrix* nano_porous_matrix = new MakeMatrix();
    nano_porous_matrix->makePorousMatrix(20,30,max_time);
    CellSolver* solver = nano_porous_matrix->mySolver;
    double T = 1.05*solver->T0;
    double dt = 0.01;
    int num_time = 100;
    string filename = "first_nano_porous";
    bool writeVMD = true;
    bool writeMeasurements = false;
    bool berendsen= true, andersen = false;
    //thermalize
    solver->thermostatON = true;
    solver->solve(0,num_time,dt,filename,writeVMD, writeMeasurements, T, berendsen, andersen);
    solver->thermostatON = false;

    //without thermostat
    solver->solve(num_time, max_time-num_time, dt,filename, writeVMD, writeMeasurements, T, berendsen, andersen);

}


void Simulation::try_config(){
    CellSolver* s = new CellSolver("/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/simulations_config/first_config_file.cfg");
    s->solve();

}
