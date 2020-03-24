import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import special
from scipy import integrate
from collections import defaultdict
import numpy as np
import os
import time
from matplotlib2tikz import save as tikz_save
#import seaborn as sns


###############################################################	
def read_data(R, directory, namefile):
	myfile = open(directory+namefile+".txt", 'r')

	for i in myfile:
		x1, x2 = i.split()
		R.append(float(x2))
	myfile.close() 

####################################################################################################
def plot_Rvspjump(p_jump, median_R, max_median_R, min_median_R, R_th, directory, namefile):
	
	plt.close()
	c = "r" 
	
	plt.plot(p_jump, median_R, color=c, marker="o", markersize=8, linestyle='', label="Simulations")
	plt.fill_between(p_jump, min_median_R, max_median_R, facecolor='gray', alpha=0.5)
	#plt.plot(p_jump, R_th, color="k", label="Theory")
	
	#plt.xlim(0, 1.5)
	#plt.ylim(0.7, 1.3)
	plt.ylabel(r"$R$", fontsize=20)
	plt.xlabel(r'p$_{jump}$', fontsize=20)
	#plt.xscale('log')
	#plt.yscale('log')
	#plt.legend(loc=2)
	plt.tight_layout()
	tikz_save(directory+namefile+".tex")
	plt.savefig(directory+namefile+".png")
	plt.draw()
	
####################################################################################################
def plot_timevsR(median_R, R_th, directory, namefile):
	
	plt.close()
	c = "r" 
	
	plt.plot(range(len(median_R)), median_R, color=c, marker="o", markersize=8, linestyle='', label="Simulations")
	plt.plot(range(len(R_th)), R_th, color="k", label="Theory")
	
	#plt.xlim(0, 1.5)
	#plt.ylim(0.7, 1.3)
	plt.ylabel(r"$R$", fontsize=20)
	plt.xlabel(r'$t$', fontsize=20)
	#plt.xscale('log')
	#plt.yscale('log')
	#plt.legend(loc=2)
	plt.tight_layout()
	tikz_save(directory+namefile+".tex")
	plt.savefig(directory+namefile+".png")
	plt.draw()

######################################################################################
def p_in(c, vexp, p_jump):
	return (1.-np.exp(-c*vexp))/(1.-np.exp(-c*vexp)+p_jump)
	
######################################################################################
def tau_out(c, vexp):
	return 1./(1.-np.exp(-c*vexp))


'''
#########################################################################################
def Epidemic_curve(I_old, mu, Lambda, r, c, vexp, p_jump, Rc, D, N):
	#p_jump = 0.3
	#print p_in(c, vexp, p_jump), 1.-p_in(c, vexp, p_jump)
	#print  mu*I_old, Lambda*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(pow(D,2.)-np.pi*pow(Rc,2.)))
	#print np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(pow(D,2.)-np.pi*pow(Rc,2.)))
	#print Lambda*I_old*(1.-I_old)
	#return I_old - mu*I_old + Lambda*I_old*(1.-I_old)*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(pow(D,2.)-np.pi*pow(Rc,2.)))
	#print p_in(c, vexp, p_jump)
	#print np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))
	#kinf = N*I_old*2.5*pow(10, -6)#*np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))
	#print (pow(Rc,2.)-2.*pow(r,2.))/(pow(Rc,2.))
	#kinf = N*I_old*np.pi*pow(r,2.)/(np.pi*pow(Rc,2.))#*(pow(Rc,2.)-2.*pow(r,2.))/(pow(Rc,2.))#*np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))
	#kinf = N*I_old*np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)+2.*np.pi*Rc*r)
	#print kinf
	#kinf = N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(pow(D,2.)-np.pi*pow(Rc,2.)))
	#print N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(pow(D,2.)-np.pi*pow(Rc,2.)))
	#kinf = N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)))#+pow(1.-p_in(c, vexp, p_jump),2.)/(1./c))
	#print np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)))
	#print np.pi*pow(r,2.)*(pow(1.-p_in(c, vexp, p_jump),2.)*p_jump/(np.pi*pow(1./(c)+Rc,2.)-np.pi*pow(Rc,2.)))
	#print np.pi*pow(r,2.)*(pow(1.-p_in(c, vexp, p_jump),2.)*p_jump/(np.pi*(pow(1./(c),2.)+Rc/c)))
	
	#print N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)))
	kinf =  N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(np.pi*pow(1./(c)+Rc,2.)-np.pi*pow(Rc,2.)))
	kinf =  N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(2.*np.pi*(pow(1./(c),2.)+Rc/c)))
	kinf =  N*I_old*(np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+2.*r*pow(1.-p_in(c, vexp, p_jump),2.)/(2.*np.pi*(1./(c)+Rc)*c))
	kinf =  N*I_old*np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+pow(1.-p_in(c, vexp, p_jump),2.)/(np.pi*pow(1./(c)+Rc,2.)-np.pi*pow(Rc,2.)))
	kinf =  N*I_old*(np.pi*pow(r,2.)*pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.))+2./c/vexp/c/vexp*r*r*pow(1.-p_in(c, vexp, p_jump),2.)/(2.*np.pi*(1./(c)+Rc))/(2.*np.pi*(1./(c)+Rc)))
	#print I_old - mu*I_old + (1.-I_old)*(1.-pow(1.-Lambda, kinf))
	kinfin =  np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)))#N*I_old*np.pi*pow(r,2.)*p_in(c, vexp, p_jump)/(np.pi*pow(Rc,2.))
	#kinfout =  N*I_old*2./c/vexp/c/vexp*r*r*(1.-p_in(c, vexp, p_jump))/(2.*np.pi*(1./(c)+Rc))/(2.*np.pi*(1./(c)+Rc))
	kinfout = 0.000506222# r/(vexp*Rc)*(1.-pow(p_in(c, vexp, p_jump),2.))
	#print kinfin
	#print kinfout 
	
	xjumpi = lambda x: c*c*np.exp(-2.*c*x)
	xjump = lambda x: c*np.exp(-c*x)*(np.exp(-c*(x-r))-np.exp(-c*(x+r)))/(Rc+x)
	xjump2 = lambda x: c*np.exp(-c*x)
	
	f = lambda y, x: c*np.exp(-c*x)*(np.exp(-c*(x-r))-np.exp(-c*(x+r)))*((r-y)/(Rc+x+y)+(r-y)/(Rc+x-y))#goood
	f3 = lambda z, y, x: x*c*np.exp(-c*x)*(c*np.exp(-c*y))*((r-z)/(Rc+x+z)+(r-z)/(Rc+x-z))
	#print Rc*c
	print integrate.quad(xjump, 0., np.inf)[0]*r/(vexp*np.pi)*(1.-pow(p_in(c, vexp, p_jump),2.))
	print integrate.dblquad(f, 0, np.inf, lambda x: 0, lambda x: r)[0]/(vexp*np.pi)*(1.-pow(p_in(c, vexp, p_jump),2.))#good
	print integrate.tplquad(f3, 0, np.inf, lambda x: x-r, lambda x: x+r, lambda x, y: 0, lambda x, y: r)[0]*2*c/(vexp*np.pi)*(1.-pow(p_in(c, vexp, p_jump),2.))
	kinfout = integrate.tplquad(f3, 0, np.inf, lambda x: x-r, lambda x: x+r, lambda x, y: 0, lambda x, y: r)[0]*2*c/(vexp*np.pi)*(1.-pow(p_in(c, vexp, p_jump),2.))#0.000506222# r/(vexp*Rc)*(1.-pow(p_in(c, vexp, p_jump),2.))
	#print integrate.quad(xjump2, 0., 1./c)[0]
	#print integrate.quad(xjump2, 0., np.inf)[0]
	#print (special.exp1(c*Rc)*c*Rc*np.exp(c*Rc)+1.)*r/(vexp*np.pi)
	kinf = N*I_old*kinfin + N*I_old*kinfout
	#print r/(c*vexp*np.pi*(1./(c)+Rc))#integrate.quad(xjump2, 0., 1./c)[0]*integrate.quad(xjump2, 0., 1./c)[0]/2.
	#print (1.-pow(p_in(c, vexp, p_jump),2.))
	#print r/(vexp*Rc)*(1.-pow(p_in(c, vexp, p_jump),2.))
	#print r*c*c/(np.pi*vexp)*(1.-pow(p_in(c, vexp, p_jump),2.))
	
	#print np.pi*pow(r,2.)*(pow(p_in(c, vexp, p_jump),2.)/(np.pi*pow(Rc,2.)))
	return I_old - mu*I_old + (1.-I_old)*(1.-pow(1.-Lambda, kinf))
	#return I_old - mu*I_old + (1.-I_old)*p_in(c, vexp, p_jump)*(1.-pow(1.-Lambda, kinfin))+(1.-I_old)*(1.-p_in(c, vexp, p_jump))*(1.-pow(1.-Lambda, kinfout))

#########################################################################################
def Epidemic_curve(I_old, R_old, mu, Lambda, r, c, vexp, p_jump, Rc, D, N):
	if p_jump == 0:
		kinf = N*I_old*np.pi*pow(r,2.)/(np.pi*pow(Rc,2.))
		print np.pi*pow(r,2.)/(np.pi*pow(Rc,2.))
		return I_old - mu*I_old + (1.-R_old-I_old)*(1.-pow(1.-Lambda, kinf)), R_old + mu*I_old
	elif p_jump==1:
		kinf = N*I_old*np.pi*pow(r,2.)/(pow(D,2.))
		print np.pi*pow(r,2.)/(pow(D,2.))
		return I_old - mu*I_old + (1.-R_old-I_old)*(1.-pow(1.-Lambda, kinf)), R_old + mu*I_old
	else:
		return 0
'''
		
#########################################################################################
def Epidemic_curve(I_old, R_old, mu, Lambda, r, c, vexp, p_jump, Rc, D, N):
	kinf = N*I_old*np.pi*pow(r,2.)/(np.pi*pow(Rc,2.))
	#print np.pi*pow(r,2.)/(np.pi*pow(Rc,2.))
	return I_old - mu*I_old + (1.-R_old-I_old)*(1.-pow(1.-Lambda, kinf)), R_old + mu*I_old
	

#############################################################################################################
def find_median(y, y_median, y_median_max, y_median_min):
	y = np.sort(y)
	y_median.append(np.median(y))
	n = float(len(y))
	#print n
	lower = int(n/2. - 1.96*np.sqrt(n)/2.)
	higher = int(1.+n/2. + 1.96*np.sqrt(n)/2.)
	if n ==1:
		lower=1
		higher=1
	#print lower, higher
	y_median_min.append(y[lower-1])
	y_median_max.append(y[higher-1])

###########################################################################################
def plot():
	directory = "./Data_pjump0/Output_file/Vicsek_merged/"
	directory = "./Data/Output_file/Vicsek_merged/"
	
	#	D		L	Rcmin	Rcmax	rmin	rmax	omega	gamma	vmin	vexp	Lambda	mu	perc_Infected	c	Steady_state	T	pjump	noise	delta_t	eta
	#1000	1000000000.	1	1000.	1000.	10.	10.	2.4	2.1	300	500	0.05	0.0001	0.01		5.     	200		1000	0.0	0.	1.	360
	
	n_runs = 100
	noise = 0
	N= 10000.
	p_jump = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	Initial_infected = 0.01
	Initial_recovered = 0.
	c = 5.
	vexp = 500.
	r = 10.
	Rc = 1000.
	D = 100000000.
	Lambda = 0.15
	mu = 0.1
	
	median_R, max_median_R, min_median_R, R, R_th = [], [], [], [], []
	
	
	for k in range(len(p_jump)):
		supp = 0.
		all_lists= []
		R = []
		for i in range(n_runs):
			if 10.*p_jump[k] < 10:
				directory = "./Data/Output_file/Vicsek_merged/"
			else:
				directory = "./Data/Output_file/Vicsek_merged/"
			namefile = "Lambda_over_mu_"+str(int(100.*p_jump[k]))+"_"+str(int(100.*noise))+"_"+str(i)
			read_data(R, directory, namefile)
			supp+=R[-1]
			all_lists.append(float(R[-1]))
			
		find_median(all_lists, median_R, max_median_R, min_median_R)
	
	'''
	print len(R[0])
	for i in range(len(R[0])):
		supp = 0.
		all_lists= []
		for j in range(n_runs):
			supp+=R[j][i]
			all_lists.append(float(R[j][i]))
			
		#print np.median(all_lists)
		#median_R.append(np.median(all_lists))
		median_R.append(np.mean(all_lists))
		#mean_R.append(supp/n_runs)
		if i==0:
			continue
		if i>1:
			supp_I, supp_R = Epidemic_curve(I_th[-1], R_th[-1], mu, Lambda, r, c, vexp, p_jump, Rc, D, N)
			I_th.append(supp_I), R_th.append(supp_R)
		else:
			supp_I, supp_R = Epidemic_curve(Initial_infected, Initial_recovered, mu, Lambda, r, c, vexp, p_jump, Rc, D, N)
			I_th.append(supp_I), R_th.append(supp_R)	
		print median_R[-1], R_th[-1]
	'''
	

	if p_jump[0] == 0:
		I_old = Initial_infected 
		R_old = Initial_recovered
		I = I_old  
		R = R_old 
		while (I>pow(10,-15)):
			I, R = Epidemic_curve(I_old, R_old, mu, Lambda, r, c, vexp, p_jump, Rc, D, N)
			I_old = I  
			R_old = R
			print I, R
			
		print I, R
		
	
	namefile = "Rvspjump"
	plot_Rvspjump(p_jump, median_R, max_median_R, min_median_R, R_th, directory, namefile)
		
	#namefile = "Epidemic_curve"
	#plot_timevsR(median_R, R_th, directory, namefile)
	

plot()	
	
	

