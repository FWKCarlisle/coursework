import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import io
import oct2py
from pathlib import Path


class output_data_spectra_dat():
    
    def get_file(self, pathname, filename, show_methods=0):

        #Find get the file name
        self.filename_path = pathname + "\\" + filename
        
        #Now we want to get the data and read it. 
        self.load_file(show_methods)
        
        if show_methods == 1:
            self.show_method_fun()
    
    def convert_time_to_epoch(self, time):
        dt = datetime.datetime.strptime(str(time), '%d.%m.%Y %H:%M:%S')
        return dt.timestamp()

    def load_file(self, show_methods):

        def read_metadata_info(meta_data):
            """
            Note: some files will not have this info e.g. data logger files
            Can be extended/edited to access more/different metadata info.
            """

            # Extract position information:
            try: self.x_pos = float(meta_data.split('X (m)\t')[1].split('\n')[0][0:-1])
            except: pass
            try: self.y_pos = float(meta_data.split('Y (m)\t')[1].split('\n')[0][0:-1])
            except: pass
            try: self.z_pos = float(meta_data.split('Z (m)\t')[1].split('\n')[0][0:-1])
            except: pass

            # Extract the time information:
            try: self.start_time = self.convert_time_to_epoch(
                meta_data[meta_data.index('Start time') + 11:meta_data.index('Start time') + 11 + 19])
            except: pass

            # Extract comment
            try: self.comment = meta_data.split('Comment01\t')[1].split('\n')[0][0:-1]
            except: pass

        #Open the file. The data and headings start after the '[DATA]' section.
        f = open(self.filename_path, "r")
        f.seek(0)
        file_all_data = f.read()
        f.close()
        
        #Split at the '[DATA]' section. The second in the list is our data with the headers.
        file_data = re.split(r"\[DATA\]\n",file_all_data)[1]
        meta_data = re.split(r"\[DATA\]\n", file_all_data)[0]

        read_metadata_info(meta_data)

        #Load as df to make life easy
        self.df = pd.read_csv(io.StringIO(file_data),delimiter = '\t')
        
        self.list_of_methods = list(self.df)
     
    def show_method_fun(self):
        
        print('Possible methods to use are:')
        
        for count, i in enumerate(self.list_of_methods):
            print(str(count) +') '+ i)
        
    def give_data(self, method_number):
        
        #if type(method_number) is int:
            #Do nothing, as this is the expected format
        if type(method_number) is str:
            #Convert to number but use the lookup table
            for count, i in enumerate(self.list_of_methods):
                if method_number in i:
                    method_number = count
                    break
            if method_number is not count:
                #We did not convert to an int value
                raise Exception("This input is not found in the list, check spelling is correct!")
                
        elif type(method_number) is not int: 
            raise Exception("Must be either and int or str variable")
        
        name_to_export = self.list_of_methods[method_number]
        
        exported_data = (self.df[name_to_export],self.list_of_methods[method_number])
        
        return exported_data
    
class Spectrum(output_data_spectra_dat):
    """
    The Spectrum class encapsulates 1 spectra file.
    It is a subclass of output_data_spectra_dat (Matt's spectra reading class)
    so it inherits attributes found on the file's metadata (eg. xy position).

    When performing some analysis to a spectrum, I like to:
        1. write the analysis within separate class
        (eg. KPFMSpectrumAnalysis)
        2. write a method that runs my analysis and stores only the
        objects of interest as Spectrum attributes (eg.
        Spectrum.KPFMAnalysis method)
    """

    def __init__(self, path, fileName):
        super().__init__()
        self.get_file(path, fileName) # load all data *and* metadata


    def ReadChannel(self, channel):
        """
        Read a channel's data.
        @param channel: channel name
        If the channel is not found, the available channels are
        printed. So, if you're not sure of the exact channel name, just
        type in nonsense.

        Note: channel = 'Index' is an option. May seem redundant, but may
        be useful in future to convert to a time channel, if we make
        note of sampling freq. TCP receiver sampling is limited to
        20 kHz. So, if measurement made using data logger through
        Matt's python_interface_nanonis.py, default is 20 kHz.
        Note: we may be able to play with TCP receiver to lower the 20kHz limit.
        @type channel: str
        @return: channel's data
        @rtype: arr
        """
        def CheckChannelExists(channel):
            if channel not in list(self.df):
                print('Choice of channel not found')
                self.show_method_fun()
                print('Index')
                sys.exit()

        if channel == 'Index':
            foo = self.give_data(0)[0] # load a random channel e.g. 0 to read its length
            x = list(range(len(foo)))
        else:
            if type(channel) == str: CheckChannelExists(channel)
            x = self.give_data(channel)[0]

        return x


    def KPFMAnalysis(self, xAdatomCentre=None, yAdatomCentre=None,
                     plotCalculation=False, axFit=None, axResiduals=None,
                     yChannel=None):
        """
        From KPFM spectra,  df(V), we want to calculate the Vcontact. That is a parabolic
         fit's, y=ax**2+bx+c, minima, b/(2a).

        Note: we can get a feel for the calculation's accuracy by:
            1. looking at the error on Vcontact by propagating the error quoted by lmfit
             for the fitting parameters b and a. *But* this is likely an underestimate
             because:
              - experimental variables play a role, eg. the quality of the AFM
                resonance, drift...
              - I don't know how lmfit calculates error. We have noticed that lmfit's
                errors appear surprisingly low.
            2. plotting it (setting plotCalculation=True)
            3. inspecting the fit's stats using lmfit.fit_report(self.fitInfo).

        @param xAdatomCentre: x coordinate of the artificial atom's centre.
        @type xAdatomCentre: float, optional
        @param yAdatomCentre: y coordinate of the artificial atom's centre.
        @type yAdatomCentre: float, optional
        If xAdatomCentre and yAdatomCentre are specified, each spectrum's distance
        from the artificial atom's centre will be calculated and added as self.r (float).
        @param plotCalculation: If True, plot of the spectrum, its found fit and its corresponding 2
        sigma conf band; the fit's minimum and its correspoding error bar derived by propagating
        the error on the fitting parameters. The default is False.
        @type plotCalculation: bool
        @param axFit:
        @type axFit:
        @param axResiduals:
        @type axResiduals:
        @param yChannel: Use if you want to analyse only one of the repeat sweeps e.g.
        'OC M1 Freq. Shift [00002] (Hz)'. Otherwise, the default is 'OC M1 Freq. Shift [AVG] (Hz)' or, in its absence,
        'OC M1 Freq. Shift (Hz)'.
        @type yChannel: str
        @return: if plotCalculation == True, the matplotlib fig and axs objects will be returned
        @rtype: a matplotlib Figure and two Axes objects
        The useful info from the analysis is added as Spectrum attributes: self.vContact (float), self.fit (arr),
        self.residuals (arr), self.dfAtVContact (float), self.vContactErr (float), self.dfAtVContactErr (float),
        self.fitInfo (Lmfit ModelResult instance), self.meanAbsRes (float), self.fitA (float)
        """
        # read the x and y channels
        self.x = self.ReadChannel('Bias calc (V)')

        if yChannel is not None:
            self.y = self.ReadChannel(yChannel)
        else:
            try: self.y = self.ReadChannel('OC M1 Freq. Shift [AVG] (Hz)')
            except: self.y = self.ReadChannel('OC M1 Freq. Shift [Hz]')

        # analyse
        kpfmAnalysis = KPFMSpectrumAnalysis(bias=self.x, df=self.y)

        if xAdatomCentre != None and yAdatomCentre != None:
            self.r = kpfmAnalysis.CalcR(self.x_pos, self.y_pos, xAdatomCentre, yAdatomCentre)

        # store useful variables
        self.vContact = kpfmAnalysis.CalcVContact()
        self.fit = kpfmAnalysis.fit
        self.residuals = kpfmAnalysis.fitInfo.residual
        self.dfAtVContact = kpfmAnalysis.dfAtVContact
        self.vContactErr = kpfmAnalysis.vContactErr
        self.dfAtVContactErr = kpfmAnalysis.dfAtVContactErr
        self.fitInfo = kpfmAnalysis.fitInfo
        self.meanAbsRes = kpfmAnalysis.meanAbsRes
        self.fitA = kpfmAnalysis.fitA


        if plotCalculation == True:
            fig, axFit, axResiduals = kpfmAnalysis.PlotVContactCalculation(axFit, axResiduals)
            return fig, axFit, axResiduals



    def ForceAnalysis(self, threshold=7.5e-12, plotCalculation=False):
        """
        work in progress
        @param threshold:
        @type threshold:
        @param plotCalculation:
        @type plotCalculation:
        """
        # read the needed channels
        self.x = self.ReadChannel('Index') # we might convert this to position or time
        self.y = self.ReadChannel('OC M1 Freq. Shift (Hz)')
        self.I = self.ReadChannel('Current (A)')

        # analyse
        analysis = ForceSpectrumAnalysis(self.y, self.I, threshold)
        analysis.FindEvent()

        # store useful variables
        self.eventMask = analysis.eventMask
        self.absDI = analysis.absDI
        self.threshold = analysis.threshold

        if plotCalculation == True:
            analysis.PlotEventCalculation()

def calc_force_trapz (z, df, A, k, f_0, abs_YN = True):
    Omega = df/f_0
    dOmega_dz = np.diff(Omega)/np.diff(z)

    z = z[:-1]
    # Delta_f = Delta_f[:-1]
    Omega = Omega[:-1]
    force = np.zeros(len(z) - 2)

    for j in range(len(z) - 2):
        # start at j+1 due to pole at t=z
        t = z[j+1:]
        
        # adjust length of Omega and dOmega_dz to length of t
        Omega_tmp = Omega[j+1:]
        dOmega_dz_tmp = dOmega_dz[j+1:]

        ### Abs to stop negative values, added 11:05, 29/10/24
        # abs_YN = True
        if abs_YN:

            integral = np.trapezoid((1 + np.sqrt(A) / (8 * np.sqrt(np.pi * abs(t - z[j])))) * Omega_tmp - 
                            A**(3/2) / np.sqrt(2 * abs(t - z[j])) * dOmega_dz_tmp, t)
            
            # correction terms for t=z from [2]
            corr1 = Omega[j] * (z[j+1] - z[j])
            corr2 = 2 * (np.sqrt(A) / (8 * np.sqrt(np.pi))) * Omega[j] * np.sqrt(abs(z[j+1] - z[j]))
            corr3 = (-2) * (A**(3/2) / np.sqrt(2)) * dOmega_dz[j] * np.sqrt(abs(z[j+1] - z[j]))

        else:
            
            inner = (1 + np.sqrt(A) / (8 * np.sqrt(np.pi * (t - z[j])))) * Omega_tmp - A**(3/2) / np.sqrt(2 * (t - z[j])) * dOmega_dz_tmp
        
            integral = np.trapz(inner, t)
            
            # correction terms for t=z from [2]
            corr1 = Omega[j] * (z[j+1] - z[j])
            corr2 = 2 * (np.sqrt(A) / (8 * np.sqrt(np.pi))) * Omega[j] * np.sqrt((z[j+1] - z[j]))
            corr3 = (-2) * (A**(3/2) / np.sqrt(2)) * dOmega_dz[j] * np.sqrt((z[j+1] - z[j]))
        force[j] = 2 * k * (corr1 + corr2 + corr3 + integral)
    return force


def file_to_force(file_path, output_path = None):

    if not os.path.exists(file_path):
        print("File path does not exist")
        raise FileNotFoundError  
    path = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    # print(path, file_name)

    spectrum = Spectrum(path, file_name)

    z = np.array(spectrum.ReadChannel('Z rel (m)'))
    index = spectrum.ReadChannel('Index')
    df = np.array(spectrum.ReadChannel('OC M1 Freq. Shift (Hz)'))
    
    amplitude = 100e-12 #1 angstrom Amplitude of the oscillation (m)
    k_spring = 1900 #Spring constant of cantilever (N/m)
    frequency_res = 20000 #Resonant frequency far from surface (Hz)
    

    Py_force = calc_force_trapz(z, df, amplitude, k_spring, frequency_res)
    
    with oct2py.Oct2Py() as oc:
        oc.addpath(r"C:\Users\ppxfc1\OneDrive - The University of Nottingham\Desktop\PhD\Code\PhD-Codes\seder-jarvis\Matlab_codes")
        M_f= oc.saderF(z, df,  frequency_res, k_spring, amplitude)

        # print("Matlab force: ", M_f)
        
        
    z_force = z[:len(Py_force)] 
    # print(len(z), len(Py_force), len(M_f))
    # print("Python force, matlab force, difference, Z")
    # for i in range(len(Py_force)):
    #     print(Py_force[i], np.real(M_f[i]), (Py_force[i] - M_f[i]),z_force[i])
    print("Python force, matlab force")
    print(Py_force[-1], M_f[-1])
    M_f = M_f.ravel()



    difference = Py_force - np.real(M_f)
    # print(difference)
    # print("Z, Py_force, M_f, difference")
    # for i, value in enumerate(difference):
    #     if i % 10 == 0:
    #         print(z_force[i], Py_force[i], M_f[i], value)
    fig, [dfAx, forceAx, differenceAx] = plt.subplots(3,1, sharex = True)
    
    plt.title(file_name)
    dfAx.plot(z, df, label = 'Frequency shift')
    dfAx.set_ylabel('Frequency shift (Hz)')
    dfAx.legend()

    forceAx.plot(z_force, np.real(Py_force)*1e9, label = 'Python force')
    forceAx.plot(z_force, np.real(M_f)*1e9, label = 'Matlab force')
    forceAx.set_ylabel('Force (nN)')
    forceAx.legend()

    differenceAx.plot(z_force, difference*1e9, label = 'Difference')
    differenceAx.set_ylabel('Difference (nN)')
    differenceAx.legend()

    plt.xlabel('Z (m)')
    plt.show()
    plt.close()
    if output_path == None:
        output_path = os.path.join(path, file_name[:-4] + '_force.txt')
    if os.path.exists(output_path):
        with open(output_path, 'w') as f:
            f.write('Z (m) , Force (N) , ' + file_name + '\n' )
            for j in range(len(z_force)):
                f.write(str(z_force[j]) + " , "+ str(Py_force[j]) +'\n')
    else:
        with open(output_path, 'x') as f:
            f.write('Z (m) , Force (N) , ' + file_name + '\n' )
            for j in range(len(z_force)):
                f.write(str(z_force[j]) + " , "+ str(Py_force[j]) +'\n')
    print("Force calculated and saved to ", output_path)
    return 1



if __name__ == "__main__":
    data_dir =r"C:\Users\ppxfc1\OneDrive - The University of Nottingham\Desktop\PhD\Code\coursework"
    file_names = [f.name for f in Path(data_dir).iterdir() if f.is_file()]
    file_names.sort()
    for file in file_names:
        if not file.endswith('.dat'):
            continue
        if not file[0] == "Z":
            continue
        file_path = data_dir + '\\' + file
        output = file_path[:-4] + '_force.txt'
        file_to_force(file_path, output)

    
   