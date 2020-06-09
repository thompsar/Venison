import numpy as np
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from scipy import optimize

from src.background_models import bg3D

class spectrum():
    """Class that loads and parses Bruker DEER spectrum files, assumes DTA and DSC files present
    
    Parameters
    ----------
    path: Full path to .DSC or .DTA file
    
    Attributes:
    metadata: (dict) full metadata from DSC file, portions of which are explicitly saved as attributes
    real: (1D or 2D array)
    imaginary: (1D or 2D array)
    abcissa: (1D array) the x-axis for the file (time (ns) for FLTPR and wavelength (nm) for SUPR)
    ordinate: (2D array) the y-axis for the data (data,slice_number)
    
    
    """
    
    def __init__(self,path):
        """Instantiate the object
        Args:
            path (str): Full path to either the DSC or DTA file
        """
        
        self.path = path
        self.metadata = {}
        self.parse_metadata()
        self.load_data()
        self.set_defaults()
        
    def parse_metadata(self):
        
        self.metadata = dict.fromkeys(['XMIN','XWID','XPTS','YPTS'])
        self.metadata['YPTS'] = 1 #always at least one slice, but not listed in DSC if YPTS = 1

        with open(self.path.split('.')[0]+".DSC") as fd: #Extract parameters from DSC File
            for line in fd:
                if line.startswith(tuple(self.metadata.keys())):
                    [key,value] = line.split()
                    self.metadata[key] = float(value)
        

        for key in ['XPTS','YPTS']:
            self.metadata[key] = int(self.metadata[key])
    
    def load_data(self):
        
        raw_data=np.fromfile(self.path.split('.')[0]+".DTA",dtype=np.float64,sep="").byteswap() #Bruker data is 64bit big endian
        array_shape = (self.metadata['YPTS'],self.metadata['XPTS'])
        self.real2D=raw_data[::2].reshape(array_shape) #Real portion of data in ypts x xpts matrix
        self.imaginary2D=raw_data[1::2].reshape(array_shape) #Imaginary portion of data in ypts x xpts matrix
        havedata=np.sum(self.real2D[:,-1]>0) #find number of slices that have complete datasets (i.e. non zero ends)

        self.real2D=self.real2D[0:havedata] #Drop empty slices
        self.imaginary2D=self.imaginary2D[0:havedata] #Drop empty slices
        
        
        
        self.metadata['YPTS'] = havedata #correct YPTS, is this the correct thing to do or will it cause issues?
        
        self.dt = self.metadata['XWID']/(self.metadata['XPTS']-1)/1000 #delta t in µs
        self.raw_time=(np.arange(self.metadata['XPTS'])*self.dt+self.metadata['XMIN']/1000) #Time axis of raw data in µs
        self.cutoff =self.raw_time[-1] #aribitrary start of cutoff
        self.background_start = self.raw_time[int(self.raw_time.size/2)] #arbitrary choice of start of background region
        
        self.slice_mask = np.ones(self.metadata['YPTS'],dtype=bool)
        
        self.avg_slices()
        
    def set_defaults(self):
        self.dimensions = 512
        self.rmin = 1.5
        self.rmax = 10
        
        self.alphas = np.logspace(-4,6, 51)
    
        
    def delete_slice(self,index):
        self.slice_mask[index] = False
    
    def avg_slices(self):
        
        self.real = np.sum(self.real2D[self.slice_mask],axis=0)
        maximum = np.max(self.real)
        self.real = self.real/maximum
        self.imaginary = np.sum(self.imaginary2D[self.slice_mask],axis=0)/maximum
        
    def build_kernel(self):
        
        ny0 = 52.04                                          # dipolar frequency at 1 nm for g=ge
        w0  = 2*np.pi*ny0                                    #angular frequency
        self.full_t   = np.arange(0., self.tmax+2*self.dt,self.dt)         #time axis in μs, +2*dt to ensure tmax is included
        self.r   = np.linspace(self.rmin,self.rmax,self.dimensions) #distance axis in nm
    
        self.full_kernel = np.zeros((len(self.r),len(self.full_t)))

        rarray=(w0/self.r[:,np.newaxis]**3)
        tarray=self.full_t[np.newaxis,:]
        kslice=rarray*tarray
        for l in range(1001):
            self.full_kernel=self.full_kernel+np.cos(kslice*(3*(l/1000)**2-1))

        self.full_kernel=self.full_kernel/self.full_kernel[:,0][:,np.newaxis]
        self.full_kernel=self.full_kernel.T
        
        #initialize waveform
        self.waveform = np.zeros(self.full_kernel.shape[0])

        self.Lmatrix = np.zeros((self.dimensions-2,self.dimensions))
        for i in range(0, self.Lmatrix.shape[0]):
            self.Lmatrix[i, 1 + (i - 1)] = 1.
            self.Lmatrix[i, 1 + (i - 0)] = -2.
            self.Lmatrix[i, 1 + (i + 1)] = 1.
    
        
        
        def hMat(alpha):
            return np.dot(self.full_kernel,
                          np.dot(
                              np.linalg.pinv(np.dot(self.full_kernel.T,self.full_kernel)+alpha**2*np.dot(self.Lmatrix.T,self.Lmatrix)),
                      self.full_kernel.T))
        
        self.full_hMats=np.asarray([np.diag(hMat(alpha)) for alpha in self.alphas])
        
        #Kernel, t and hMats may be subsets of the full computed version.
        self.kernel = self.full_kernel
        self.t = self.full_t
        self.hMats = self.full_hMats
        
        #initialize tikhonov parameters
        self.bvector = np.concatenate([self.waveform,
                                       np.zeros(self.Lmatrix.shape[0]) ]) #insert reference
        
        self.nnls_solutions = np.zeros([self.alphas.shape[0],self.kernel.shape[1]])
    
    def fit_background(self,model,background_range):
        bgfit=optimize.curve_fit(model,self.raw_time[background_range],
                         self.real[background_range])[0]
        self.background = model(self.raw_time,*bgfit)
        return self.background
        
    
    def background_correct(self,background = None):
        if background is None:
            background = self.background[self.waveform_range]
        else:
            background = background[self.waveform_range]
            
        self.normamp=(self.zeroamp-background[0])/(background[0])
        self.waveform=(self.real[self.waveform_range]-background)/background
        self.waveform=self.waveform/self.normamp
        
        self.bvector = np.concatenate([self.waveform,
                                       np.zeros(self.Lmatrix.shape[0])])
        
        if self.kernel.shape[0]!=self.waveform.shape[0]:
            self.resize_kernel()
        return self.waveform
            
        
        
    def resize_kernel(self):
        new_size = self.waveform.shape[0]
        self.kernel=self.full_kernel[:new_size,:]
        self.t=self.full_t[:new_size]
        self.hMats=self.full_hMats[:,:new_size]
    
    def update_ranges(self):
        
        self.waveform_range = range(abs(self.raw_time-self.zeropoint).argmin(),
                       abs(self.raw_time-self.cutoff).argmin()+1) #plus one to account for non-inclusive indexing
        
        self.background_range = range(abs(self.raw_time-self.background_start).argmin(),
                         abs(self.raw_time-self.cutoff).argmin()+1)
        
        
    def alpha_regularize(self, alpha):
        return optimize.nnls(np.concatenate([self.kernel,alpha*self.Lmatrix]),self.bvector)[0]
    
    def waveform_regularize(self, waveform):
        return optimize.nnls(np.concatenate([self.kernel,self.alpha*self.Lmatrix]),
                             np.concatenate([waveform,
                                       np.zeros(self.Lmatrix.shape[0])]))[0]
    
    def generate_lcurve(self):
        
        pool = mp.Pool()
        results = [pool.apply_async(self.alpha_regularize, args=(alpha,)) for alpha in self.alphas]
        self.solutions = [p.get() for p in results]
        self.solutions=np.array(self.solutions)
        pool.close()
        
        self.tikhonovfits = np.dot(self.kernel,self.solutions.T).T #generate all the spectra from solutions
        self.tikhonovfits = self.tikhonovfits/self.tikhonovfits[:,0][:,np.newaxis] #use broadcasting to normalize solutions
        
    def validate_background(self, background_start = 0.2 , n_backgrounds = 100, end_offset = 30):
        
        """Determine optimal background fit
    
        Parameters
        ----------
        background_start: starting position for bg validation
        n_backgrounds: number of unique backgrounds to test
        end_offset: offset (in int index) from cuttoff point, want >1

        """

        cutindex = abs(self.raw_time-self.cutoff).argmin()
        bgstartindex = abs(self.raw_time-background_start).argmin()
        

        bgstepsize = np.floor((cutindex - end_offset - bgstartindex)/(n_backgrounds))
        if bgstepsize == 0:
            bgstepsize = 1

        bgpositions = range(bgstartindex, self.background_range[-1] + 1 - end_offset, int(bgstepsize))

        self.background_fit_points = self.raw_time[bgpositions]
        
        self.backgrounds = np.array([self.fit_background(bg3D,range(i,self.background_range[-1]+1)) for i in bgpositions])
        self.waveforms = np.array([self.background_correct(background) for background in self.backgrounds])

        pool = mp.Pool()
        results = [pool.apply_async(self.waveform_regularize, args=(waveform,)) for waveform in self.waveforms]
        self.validation_solutions = [p.get() for p in results]
        self.validation_solutions = np.array(self.validation_solutions)
        pool.close()
        
        self.validation_tikhonovfits = np.dot(self.kernel,self.validation_solutions.T).T #generate all the spectra from solutions
        self.validation_tikhonovfits = self.validation_tikhonovfits/self.validation_tikhonovfits[:,0][:,np.newaxis] #use broadcasting to normalize solutions
        
        #return to user selected background/waveform...there is probably a better way to handle all this.
        self.update_ranges()
        self.fit_background(bg3D,self.background_range)
        self.background_correct()