from scitools.std import *
import sys, os, glob;
if(len(sys.argv)>2):
    project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Release/"
    #project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/LastVersionPro/FYS4460Project1/project1-build-desktop-Qt_4_8_1_in_PATH__System__Release/"
    #project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/Project2-build-Desktop-Qt_4_8_1_in_PATH__System__Release/"
    results_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/"
    base_name = sys.argv[1]
    num_files = eval(sys.argv[2])
    glued_name = "glued_" + base_name + ".xyz";
    glued_file = open(results_dir + glued_name, 'w');

    current = 0;
    

    for i in range(num_files):
        filename = base_name + str(current)+ ".xyz";
        infile = open(project_dir+filename, 'r');
        filestring = infile.read();
        glued_file.write(filestring);
        current+=1;
        infile.close();
    glued_file.close();
    if(len(sys.argv)>3):
	for i in glob.glob("%s*.xyz"%(project_dir +base_name)):
		os.remove(i)
	
else:
    print "you forgot the commandline arguments.\nsys.argv[1] = basename\nsys.argv[2] = number of files"
