# read an optotrack Volatage file.
import struct
import numpy as np
import matplotlib.pyplot as plt

class OptoHeader:
    FType = 1
    NItems = 0
    FItems = 0
    NFrames = 0
    Freq = 0
    Junk = 0
    DataMatrix = 0


def calcRMS(signal):

    #Returns the RMS for a 2d signal (returns the RMS of the columns
    #Means are subtracted before returning the RMS

    #Subtract the mean
    signal = signal - signal.mean(axis = 0)
    signal = np.power(signal, 2)
    signal = signal.mean(axis = 0)
    signal = np.power(signal, 0.5)
    return signal


def calcMPF(inp_signal, SR):

    #Calculates the PSDF for a 1D signal and returns the mean power freq.
    #Inputs are the 1D signal, and the sample rate of the signal.

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import signal

    #subtract Mean
    inp_signal = inp_signal - inp_signal.mean(axis = 0)
    #plt.plot(inp_signal)
    #plt.show()

    N = len(inp_signal)

    #sample spacing
    T = 1 / SR
    x = np.linspace(0.0,N*T, N)

    f, yf = signal.welch(inp_signal, SR, axis=0, scaling='spectrum')

    #xf = np.linspace(0,1/(2.0*T), N//2)


    #Uncomment below to get a plot of the power spectrum
    #plt.plot(f, yf)
    #plt.grid()
    #plt.show()

    top_sum = 0
    bot_sum = 0
    for idx in range(0,len(f)):
        top_sum = top_sum + f[idx]*yf[idx]
        bot_sum = bot_sum + yf[idx]

    MPF = top_sum/bot_sum

    return MPF

def ButterWorthFilt2(inp_signal, fs, fc):
    #returns a dual pass second order (4th order) butterworth filter for a 1D signal.
    #inputs are the signal, sampling frequency and the cutt-off frequency

    from scipy import signal
    import numpy as np
    # import matplotlib.pyplot as plt ---> uncomment if a plot is required.

    inp_signal = np.asarray(inp_signal)
    #plt.plot(inp_signal, label = 'pre') # plot the input signal

    #specify the filter properties
    nyq = 0.5 * fs
    N = 2 # Filter order
    w = (fc/nyq) #normalized cuttoff
    b, a = signal.butter(N, w) # order, w,

    #filter the signal.
    y = signal.filtfilt(b, a, inp_signal, axis = 0)
    #plt.plot(y, label = 'post') # plot the resulting figure
    #plt.show()
    return y

def ReadFPFile(Tracking_File):

    #open the binary Optotrack file
    # Tracking_File = 'v1#022.dat'
    file = open(Tracking_File, 'rb')

    FileInfo = OptoHeader()
    FileInfo.FType = file.read(1)
    FileInfo.NItems = int.from_bytes(file.read(2), byteorder='little')
    FileInfo.FItems = int.from_bytes(file.read(2), byteorder = 'little')
    FileInfo.NFrames = int.from_bytes(file.read(4), byteorder = 'little')
    FileInfo.Freq = struct.unpack('f', file.read(4))
    FileInfo.Junk = file.read(243)
    #FileInfo.DataMatrix = file.read()

    #Check Data
    print('File Name: ',Tracking_File)
    print("FType: ", FileInfo.FType)
    print("NItems: ", FileInfo.NItems)
    print("FItems: ",FileInfo.FItems)
    print("NFrames: ",FileInfo.NFrames)
    print("Freq: ", FileInfo.Freq)

    Sorted_Matrix = np.zeros((FileInfo.NFrames,FileInfo.NItems * FileInfo.FItems))
    for row in range(FileInfo.NFrames):
        for col in range(FileInfo.NItems*FileInfo.FItems):
            Temp_Value = struct.unpack('f', file.read(4))
            Sorted_Matrix[row,col] = Temp_Value[0]

            #Close the tracking file
    file.close()

    return(FileInfo, Sorted_Matrix)
#np.savetxt("test.csv", Sorted_Matrix, delimiter = ',')


def ReadKinoveaBarPosition():
    import pandas as pd
    import os 
    import tkinter as tk
    from tkinter import filedialog
    import numpy as np
    
    root = tk.Tk()
    root.withdraw()
    FilePath = filedialog.askopenfilename()
    Path, File = os.path.split(FilePath)
    #print(Path, File)
    os.chdir(Path)
    xl = pd.ExcelFile(File)
    AllData = xl.parse(0)
    Sheet_Info = AllData.shape
    BarData = AllData.iloc[9:Sheet_Info[0],0:2] #Pull out the bar data
    BarPosition = BarData.to_numpy()

    return(BarPosition)


#Function to Correct Treadmill data from Xsens
def Treadmill_Correction(DataIn):
    from scipy import signal
    Output = []
    total_frames = len(DataIn)
    for ii in range(total_frames):

        #Determine the angular offset of the x-axis (angle of progression)
        Angle = np.arctan2(DataIn[ii][1],DataIn[ii][0]) #1 = y , 0 = X
        Rot_Matrix = [[np.cos(Angle), -np.sin(Angle), 0],
        [np.sin(Angle), np.cos(Angle), 0],
        [0, 0, 1]]
        Output.append(np.matmul(DataIn[ii],Rot_Matrix))

    #Pull out the XColumn
    X_Col = [row[0] for row in Output]
    X_Col = signal.detrend(X_Col)
   #print(X_Col)

    #replace Xcolumn with detrended version
    for row, item in enumerate(Output):
        Output[row][0] = X_Col[row]
    return(Output)



def rubber_band(SignalIn, NPoints):
    from scipy.signal import resample
    f = resample(SignalIn, NPoints)
    return(f)

def QuickPlot(DataIn):
    import matplotlib.pyplot as plt
    plt.plot(range(len(DataIn)), DataIn)
    plt.show()

def LoadMVNX():
    import os
    import tkinter as tk
    import mvnx
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    FilePath = filedialog.askopenfilename()
    Path, File = os.path.split(FilePath)

    #Read the Position Data
    os.chdir(Path)
    data = mvnx.MVNX(File)
    position_info = data.get_info('position')
    #mvnx.plotData(position_info, 'right_hand')
    return(position_info)

def GetPositionData(data, tag):
    total_frames = len(data)
    DataOut = []
    for i in range(total_frames): DataOut.append(data[i][tag])

    return(DataOut)

def LoadMVNX_Angles():
    import os
    import tkinter as tk
    import mvnx
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    FilePath = filedialog.askopenfilename()
    Path, File = os.path.split(FilePath)

    #Read the Position Data
    os.chdir(Path)
    data = mvnx.MVNX(File)
    angle_info = data.get_info('jointAngleXZY')
    #mvnx.plotData(position_info, 'right_hand')
    return(angle_info, File)

def GetAngleData(data, tag):
    total_frames = len(data)
    DataOut = []
    for i in range(total_frames): DataOut.append(data[i][tag])

    return(DataOut)

def SHOT_PLOT(mean_values, std_dev_values, REL_mean_downswing_position, REL_mean_contact_position):
    fig, ax = plt.subplots(facecolor='darkslategrey')

    # Convert lists to numpy arrays for element-wise operations
    mean_values = np.array(mean_values)
    std_dev_values = np.array(std_dev_values)
    
    # Plot the mean line
    ax.plot(mean_values, label='Mean', color='hotpink', linewidth=2)

    # Fill the area between (mean - std) and (mean + std) with a colored background
    ax.fill_between(range(len(mean_values)), mean_values - std_dev_values, mean_values + std_dev_values, color='aqua', alpha=0.5, label='Standard Deviation')

    # Fill the thick band at the bottom of the graph for Zone 1 (0-REL_mean_downswing_position)
    ax.fill_between(range(0, int(REL_mean_downswing_position)), min(mean_values - std_dev_values), y2=min(mean_values - std_dev_values), linewidth=1, color='red', alpha=1, label='BackSwing')

    # Fill the thick band at the bottom of the graph for Zone 2 (REL_mean_downswing_position-REL_mean_contact_position)
    ax.fill_between(range(int(REL_mean_downswing_position), int(REL_mean_contact_position)), min(mean_values - std_dev_values), y2=min(mean_values - std_dev_values), linewidth=1, color='yellow', alpha=1, label='DownSwing')

    # Fill the thick band at the bottom of the graph for Zone 3 (REL_mean_contact_position-end)
    ax.fill_between(range(int(REL_mean_contact_position), len(mean_values)), min(mean_values - std_dev_values), y2=min(mean_values - std_dev_values), linewidth=1, color='green', alpha=1, label='Follow-Thru')

    # Add a vertical line at x = REL_mean_contact_position
    ax.axvline(x=int(REL_mean_contact_position), color='yellow', linestyle='-.', label='Point of Puck Contact')
    # Add a vertical line at x = REL_mean_downswing_position
    ax.axvline(x=int(REL_mean_downswing_position), color='red', linestyle='-.', label='Transition from Back Swing to Down Swing')

    # Customize the plot
    ax.set_xlabel('% Shooting Cycle', fontsize=12)
    ax.set_ylabel('Some Joint (Degrees)', fontsize=12)

    # Commented out Y-axis limiter
    # ax.set_ylim(-1.5, 1.5)  # Adjust these values as needed

    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=7, facecolor='darkslategrey', edgecolor='none', borderaxespad=0, fontsize=8)
    for text in legend.get_texts():
        text.set_color("white")  # Set color for all legend texts

    plt.subplots_adjust(bottom=0.2)

    ax.set_facecolor('darkslategrey')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')

    # Show the plot
    plt.show()

# Example usage:
# SHOT_PLOT(mean_values, std_dev_values, REL_mean_downswing_position, REL_mean_contact_position)