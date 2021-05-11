# -- script that read Luna data and various utilities
#
# Philippe Jean
# 19-11-2019
#
#
import sys, os
sys.path.append('/Users/philippejean/Doctorat/Python/')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.constants as cte
import scipy.io

# -- class
class luna():

	def __init__(	self,
					filename):

		self.filter_bw = 30e-12 	# Time domain filter bandwidth [s]

		self.filename = filename
		# self.param_dict = dict

		# -- load file into class
		fn, self.filetype = os.path.splitext(self.filename)
		if self.filetype == '.txt':
			self.load_file(9)
		elif self.filetype == '.mat':
			self.loadMat()
		else:
			print('File extension not supported.')
			exit()
		# self.load_bin()

		# -- create parameters list dict
		self.param_list = {
			'wvl'	: 'X Axis - Wavelength (nm)',
			'f' 	: 'X Axis - Frequency (GHz)',
			'il' 	: 'Insertion Loss (dB)',
			'gd' 	: 'Group Delay (ps)',
			'max' 	: 'Max Loss (dB)',
			'min'	: 'Min Loss (dB)',
			'jma_a' : 'JM Element a Amplitude',
			'jma_b' : 'JM Element b Amplitude',
			'jma_c' : 'JM Element c Amplitude',
			'jma_d' : 'JM Element d Amplitude',
			'jmp_a' : 'JM Element a Phase (rad)',
			'jmp_b' : 'JM Element b Phase (rad)',
			'jmp_c' : 'JM Element c Phase (rad)',
			'jmp_d' : 'JM Element d Phase (rad)'
							}

		# -- load frequency grid
		self.f 		= self.read_param(param='f')
		self.w 		= 2*np.pi*self.f
		self.wvl 	= cte.c / self.f
		self.df 	= np.abs(np.mean(np.diff(self.f)))
        
		# -- define empty jones matrix
		self.jm = np.zeros([2, 2, len(self.f),], dtype=complex)
        
		return

	def loadMat(self):
		# Load the data from a mat file containing the fields
		self.data = scipy.io.loadmat(self.filename)

	# def load_bin(self):

	# 	f = open(self.filename, 'r')

	# 	self.start_freq 		= self.fread(1, 'float64')[0]
	# 	self.samp_freq 			= self.fread(2, 'float64')[1]
	# 	self.segment 			= self.fread(12, 'uint32')[6]
	# 	self.l_dut	 			= self.fread(20, 'float64')

	# 	f.close()

	# 	return

	# def fread(self, nelements, dtype):

	#     """Equivalent to Matlab fread function"""

	#     if dtype is np.str:
	#         dt = np.uint8
	#     else:
	#         dt = dtype

	#     data_array = np.fromfile(self.filename, dt, nelements)
	#     data_array.shape = (nelements, 1)

	#     return data_array

	def load_file(self, nskip=8):
		print('Loading file %s ...\n' %self.filename)
		self.file = open(self.filename, 'r')
		for ii in range(0, nskip):
			temp = self.file.readline()
		self.head = temp.split('\t')
		self.head[-1] = self.head[-1].replace('\n', '')
		self.file.close()
		return

	def read_param(self, param='f'):
		''' read a single parameter from the data file '''
		if self.filetype == '.txt':
			param = self.param_list[param]
			if param in self.head:
				for ii in range(0,len(self.head)):
					if self.head[ii] == param:
						iparam = ii
				dat = np.loadtxt(self.filename, skiprows=9)[:,iparam]
			else:
				print(' \' param \' %s not found' %param)
				dat = 0
			return dat
		elif self.filetype == '.mat':
			return np.squeeze(self.data[param])

	def calc_il(self):
		''' calculate averaged insertion loss from jones matrix '''
		self.il = 10*np.log10(( np.abs(self.jm[0,0])**2 + np.abs(self.jm[0,1])**2 + np.abs(self.jm[1,0])**2 + np.abs(self.jm[1,1])**2 ) /2)
		fig, ax0=plt.subplots(figsize=(4,3))
		ax0.plot(self.wvl, self.il, 'k', lw=1)
		ax0.set_xlabel('Wavelength [nm]')
		ax0.set_ylabel('Mean insertion loss [dB]')
		plt.tight_layout()
		plt.show()
		return

	def calc_phase(self):
		''' calculate the phase shift across '''
		vector = np.angle((
                            self.jm[0,0,1::]*np.conj(self.jm[0,0,0:-1]) + 
							self.jm[1,0,1::]*np.conj(self.jm[1,0,0:-1]) +
							self.jm[0,1,1::]*np.conj(self.jm[0,1,0:-1]) +
							self.jm[1,1,1::]*np.conj(self.jm[1,1,0:-1]) 
                            ))
		phase = []
		for i in range(1, len(vector)):
			phase.append(np.sum(vector[0:i]))
		self.phase = phase
		return

	def calc_gd(self):
		''' calculate averaged group delay from jones matrix '''
		self.gd = np.angle((
                            self.jm[0,0,1::]*np.conj(self.jm[0,0,0:-1]) + 
							self.jm[1,0,1::]*np.conj(self.jm[1,0,0:-1]) +
							self.jm[0,1,1::]*np.conj(self.jm[0,1,0:-1]) +
							self.jm[1,1,1::]*np.conj(self.jm[1,1,0:-1]) 
                            )) / np.diff(self.w)
		return

	def get_jm(self):
		''' return Jones matrix of the device '''
		if self.jm.all() == 0:
			print('Extracting Jones matrix...\n')
			self.jm[0,0] = self.read_param(param='jma_a') * np.exp( 1j * self.read_param(param='jmp_a') )
			self.jm[0,1] = self.read_param(param='jma_b') * np.exp( 1j * self.read_param(param='jmp_b') )
			self.jm[1,0] = self.read_param(param='jma_c') * np.exp( 1j * self.read_param(param='jmp_c') )
			self.jm[1,1] = self.read_param(param='jma_d') * np.exp( 1j * self.read_param(param='jmp_d') )
		return 

	def arb_pol(self, phi=0, theta=0):
		''' 
		returns an arbitrary polarization vector for simulation 
		of the form E = [Ex; Ey]
		Ex = cos(theta) * exp(j*phi)
		Ey = sin(theta) * exp(-j*phi)
		'theta' is in [0 2*pi]
		'phi' is in [0 pi]
		'''
		ex = np.cos(theta) * np.exp(1j*phi)
		ey = np.sin(theta) * np.exp(-1j*phi)
		return np.array([[ex], [ey]])

	def simul_pol(self, phi=0, theta=0):
		self.get_jm()
		''' simulate the device response to an arbitrary polarization '''
		ein 	= self.arb_pol(phi=phi, theta=theta)
		eout_x 	= self.jm[0,0] * ein[0] + self.jm[0,1] * ein[1] # jma*ein_x + jmb*ein_y
		eout_y 	= self.jm[1,0] * ein[0] + self.jm[1,1] * ein[1] # jmc*ein_x + jmd*ein_y
		return np.squeeze(eout_x), np.squeeze(eout_y)

	def pol_gui(self):
		'''
		generate a gui to estimate the polarization axes
		'''
		self.get_jm()
		phi0 	= 0
		theta0 	= 0
		dphi 	= 0.00001
		dtheta 	= 0.00001

		self.gui_fig, ax = plt.subplots()
		plt.subplots_adjust(left=0.25, bottom=0.3)
		ex, ey = self.simul_pol(phi=phi0, theta=theta0)
		ax.plot(self.f/1000, self.applySmoothingFilter(self.read_param(param='max')), '--r', lw=0.75)
		ax.plot(self.f/1000, self.applySmoothingFilter(self.read_param(param='min')), '--b', lw=0.75)
		self.l, = ax.plot(self.f/1000, self.applySmoothingFilter(10*np.log10( np.abs(ex)**2 + np.abs(ey)**2 )), c='k', lw=1.5)
		ax.set_ylabel('Transmission [dB]')
		ax.set_xlabel('Frequency [GHz]')
		ax.margins(x=0)

		ax_phi 		= plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='silver')
		ax_theta 	= plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='silver')

		self.slider_phi 	= Slider(ax_phi, r'$\phi$ [rad]', 0.0, 2*np.pi, valinit=phi0, valstep=dphi, valfmt='%1.4f')
		self.slider_theta 	= Slider(ax_theta, r'$\theta$ [rad]', 0.0, 2*np.pi, valinit=theta0, valstep=dtheta, valfmt='%1.4f')

		self.slider_phi.on_changed(self.slider_update)
		self.slider_theta.on_changed(self.slider_update)

		plt.show()

		return

	def slider_update(self,val):
		up_phi 		= self.slider_phi.val
		up_theta 	= self.slider_theta.val
		ex, ey = self.simul_pol(phi=up_phi, theta=up_theta)
		self.l.set_ydata(self.applySmoothingFilter(10*np.log10( np.abs(ex)**2 + np.abs(ey)**2))  )
		self.gui_fig.canvas.draw_idle()

	def plot_pol_trans(self, phi=0, theta=0, scale='dB'):
		''' 
			plot the transmission response to a given polarization input 
			scale : 'dB' or 'lin'
		'''
		eout_x, eout_y = self.simul_pol(phi=phi, theta=theta)
		fig,ax0=plt.subplots(figsize=(4,3))
		if scale == 'lin':
			ax0.plot(self.f/1000, np.abs( eout_x + eout_y )**2, 'k', lw=1)
		elif scale == 'dB':
			ax0.plot(self.f/1000, 10*np.log10(np.abs( eout_x + eout_y )**2), 'k', lw=1)
		else:
			print('wrong scale, dumbass')
		ax0.set_title(r'$\phi$=%.2f rad    $\theta$=%.2f rad' %(phi,theta))
		ax0.set_xlabel('frequency [THz]')
		ax0.set_ylabel('transmission')
		ax0.set_xlim((self.f.min()/1000, self.f.max()/1000))
		ax0.grid()
		fig.tight_layout()
		plt.show()
		return

	def applySmoothingFilter(self, quantity):
		""" applies a smoothing filter to a measurement quantity . From J. Cauchon. """
		wavelength 	= self.read_param('wvl')
		time 		= self.read_param('time')

		FWHM = (self.read_param('wvl')[int(wavelength.shape[0]/2)])**2 / self.filter_bw / cte.c

		# conversion from FWHM to sigma
		sigma = FWHM/(np.sqrt(8*np.log(2)))
		gaussian = np.exp(-((time - time[int(time.shape[0]/2)])/sigma)**2)
	
		filtered_time_domain = gaussian * np.fft.fftshift(np.fft.fft(quantity))
		filtered_quantity = np.fft.ifft(np.fft.fftshift(filtered_time_domain))

		return filtered_quantity
