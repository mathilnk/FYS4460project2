#import sys,os
from Simulation import *
def run_that_works():
	config_folder = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Simulations/make_thermalize_run/"
	configFile = "make_thermalize_and_run_config1.cfg"
	configFile = "make_and_run_config.cfg";
	sim = Simulation("fil",config_folder+configFile)
	sim.run()
config_folder = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Simulations/make_thermalize_run/"
fil =  config_folder+ "make_and_run";
sim = Simulation(fil);
sim.make_crystal(8,100,5.72,"ar")
sim.thermalize(100, 100);
sim.solve(100, 'true')

