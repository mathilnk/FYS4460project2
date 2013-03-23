import numpy as np
import matplotlib.pyplot as mpl;

#(1) keep the number of processord constant P=4 Here E_p is S_p divided by 4.
time_pro_4 =np.array([9.6402e+00,1.9215e+01,3.3312e+01,5.5370e+01,8.3471e+01,1.2396e+02,1.7064e+02,2.2809e+02,2.9443e+02]);
time_pro_1 =np.array([2.6987e+01,5.2623e+01,9.2187e+01,1.4794e+02,2.2185e+02,3.1644e+02,4.3345e+02,6.0226e+02,8.0646e+02]);
number_of_atoms = 4*np.array([8,10,12,14,16,18,20,22,24])**3;
speed = time_pro_1/time_pro_4;

mpl.plot(number_of_atoms,time_pro_4, '-o');
mpl.title("Time spent as a function of number of atoms, P=4");
mpl.xlabel("Number of atoms");
mpl.ylabel("Time spent (s)");
mpl.figure();

mpl.plot(number_of_atoms, speed,'-o');
mpl.title("Speed of the simulation as a function of number of atoms, P=4");
mpl.xlabel("Number of atoms");
mpl.ylabel("Speed");
mpl.figure();

mpl.plot(number_of_atoms, speed/4,'-o');
mpl.title("Parallel efficiency as a function of number of atoms");
mpl.xlabel("Number of atoms");
mpl.ylabel("Efficiency");


#(2) keep the number of atoms constant and increase the number of processors
time_N_10 = np.array([57.2695,37.2261,32.5168,32.3311,30.2027,28.7962,27.7622,27.7292,29.0996,28.8743]);
speed_N_10=np.array([1,1.53842,1.76123,1.77134,1.89617,1.98879,2.06286,2.06531,1.96805,1.9834]);
num_of_pros= np.linspace(1,10,10);
efficiency_N_10 = speed_N_10/num_of_pros;
print num_of_pros
print time_N_10
mpl.figure();
mpl.plot(num_of_pros, time_N_10, '-o');
mpl.title("Time as a function of number of processors. N=10^3");
mpl.xlabel("Number of processors");
mpl.ylabel("Time (s)");

mpl.figure()
mpl.plot(num_of_pros, speed_N_10, '-o')
mpl.title("Speed as a function of number of processors. N=10^3");
mpl.xlabel("Number of processors");
mpl.ylabel("Speed");

mpl.figure()
mpl.plot(num_of_pros, efficiency_N_10, '-o')
mpl.title("Efficiency as a function of number of processors. N=10^3");
mpl.xlabel("Number of processors");
mpl.ylabel("Efficiency");

#(3) keep the ratio N/P constant.
#time_NP_ratio=np.array([1.4956e+01,2.9110e+01,6.3429e+01,1.3972e+02,2.6105e+02]);
#speed_NP_ratio = np.array([1.8930,1.9383,1.5740,1.1399,0.9135]);

time_NP_ratio = np.array([1.6555e+01,2.7368e+01,4.8319e+01,8.3535e+01,1.3780e+02]);
speed_NP_ratio = np.array([1.6885,2.0823,2.0616,1.9047,1.7240])
N = 4*np.array([8,10,12,14,16])**3
#P = np.array([16,31,54,85,128]);
P = np.array([4,7,13,21,32]);
eff = speed_NP_ratio/P;

mpl.figure()
mpl.plot(N,time_NP_ratio, '-o');
mpl.title("Time as a function of number of atoms (P=N/128)");
mpl.xlabel("Number of atoms");
mpl.ylabel("Time (s)");


mpl.figure()
mpl.plot(N, speed_NP_ratio, '-o');
mpl.title("Speed as a function of number of atoms, P=N/128");
mpl.xlabel("Number of atoms");
mpl.ylabel("Speed");


mpl.figure()
mpl.plot(N, eff, '-o');
mpl.title("Parallel efficiency as a function of number of atoms (P=N/128");
mpl.xlabel("Number of atoms");
mpl.ylabel("Efficiency");


mpl.show()

