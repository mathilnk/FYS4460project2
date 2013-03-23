import sys, os;


class Simulation:
    def __init__(self,v):
        self.v = v;
        self.start_state_file = "state_file.xyz";
        self.command_line_file = "cmd_line_arg.txt";

    def prepare_system(self, filename):
        self.filename = filename;
        if(filename == "None"):
            #from scratch
            self.v=0;
        else:
            #put the the content of filename into start_state_file.
            self.v = 0;
    def thermalize(self, T, N):
        #run the program with thermostat for N timesteps
        self.v = 0;
    def makePores(self, num, r_min, r_max):
        #make random spheres of random size
        self.v = 0;
    def run(self):
        #run the program with the command_line_file.
        project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release/";
        executable = "Project2";
        cmd_arg = " ./results/cut_cyllinder.xyz 8 3 5.72 100 100 0 Ar 200 4 0 1";
        terminal_run = project_dir+executable + cmd_arg;
        os.system(terminal_run);
        print terminal_run

    def solve(self, N):
        #run the program with N timesteps
        #the xyz files get filename_ending solve to separate from files made during thermalize
        self.v =0;
    def glueVMDfiles(self):
        #makes two files, one glued thermalize..., and one glued from running solve
        self.v = 0;
    def compile_program(self):
        #compiles. Should only happen once
        terminal_compile = "qmake-qt4 /home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2/Project2.pro -r -spec linux-g++-64"
        terminal_compile2 = "make -w /home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release"
        os.system(terminal_compile);
        os.system(terminal_compile2);


if __name__=="__main__":
    simulator = Simulation(6);
    simulator.compile_program();
    simulator.run();

