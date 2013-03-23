import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
def plot_radial():
	plt.figure()
	x = [0.1822*i for i in range(len(radial_dist))];
	plt.plot(x,radial_dist);
	plt.xlabel('r (AAngstrom)');
	plt.ylabel('Distribution');
	temp,pressure =  find_measurements(current_time_step);
	plt.title("g(r) for T=%s and P=%s"%(temp,pressure));
	#plt.title("g(r) for a perfect crystal");
	#plt.axis([0,200,0,120000]);
	return 0;
	
def find_measurements(time_step):
	myline = measure_lines[time_step+2].split();
	
	temp = myline[3];
	pressure = myline[4];
	return temp,pressure;
	
project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release"
if(len(sys.argv)==2):
	filename = sys.argv[1];
	infile = open(dir_p  + filename+ "0199_radial.txt", 'r');
	#measure_file = open(project_dir + filename+"0_energy.txt",'r');
	lines = infile.readlines();
	radial_dist = np.zeros(len(lines)-1);
	measure_file = open(project_dir + filename+"0_energy.txt",'r');
	measure_lines = measure_file.readlines();
	current_time_step = 198;
	for i in range(1,len(lines)-1):
		radial_dist[i] = eval(lines[i]);
	
	plot_radial();
	plt.show()
		
		
if(len(sys.argv)>3):
	filename = sys.argv[1];
	start_time_step = int(sys.argv[2]);
	stop_time_step = int(sys.argv[3]);
	infile = open(project_dir + "glued_" + filename+ "_radial.txt", 'r');
	measure_file = open(project_dir + filename+"_energy.txt",'r');
	measure_lines = measure_file.readlines();
	
	lines = infile.readlines();
	#find start
	start =0;
	N = len(lines);
	current_time_step =0;
	#find num of bins
	numOfBins = 0;
	for i in range(1,N):
		words = lines[i].split();
		if(words[0]=="new"):
			break;
		else:
			numOfBins+=1;
	#plot for the right time interval
	count = 0;
	for i in range(N):
		words = lines[i].split();
		#print stop_time_step, current_time_step;
		if(words[0]=="new"):
			current_time_step = int(words[-1]);
			if(current_time_step>=start_time_step and current_time_step<=stop_time_step):
				if(current_time_step>start_time_step):
					plot_radial();
					
				radial_dist = np.zeros(numOfBins);
				count=0;
			if(current_time_step>stop_time_step):
				break
		elif(current_time_step>=start_time_step):
			#print "her"
			radial_dist[count] = eval(words[0]);
			count +=1;
			
		

	plt.show();

