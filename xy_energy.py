# calculate DFS coefficients of each x vs. t and y vs. t signal
# calculate percent energy in peak, std. deviation of energy distribution,
# and percent energy within 1 std. deviation of w = 0

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# save interesting quantities
# fraction of energy contained in largest DFS coefficient
Epeak_x_copy,Epeak_y_copy = [],[]
Epeak_x_command,Epeak_y_command = [],[]

# fraction of energy contained in DFS coefficients that
# within 1 std. deviation of w = 0
Estd_x_copy,Estd_y_copy = [],[]
Estd_x_command,Estd_y_command = [],[]
Ecentral_x_copy,Ecentral_y_copy = [],[]
Ecentral_x_command,Ecentral_y_command = [],[]

for fname in dirs:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # copy or command clock?
    clock_type = ''
    x_copy,y_copy,t_copy = [],[],[]
    x_command,y_command,t_command = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found copy clock?
        if found_clock==False:
            if 'COPY' in line:
                clock_type = 'COPY'
                found_clock = True
            elif 'COMMAND' in line:
                clock_type = 'COMMAND'
                found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'symbol label' in line and clock_type=='COPY' and len(x_copy)>0:
                found_clock = False
                record = False
                continue
            elif 'symbol label' in line and clock_type=='COMMAND' and len(x_command)>0:
                found_clock = False
                record = False
                continue
            elif 'point' not in line:
                continue
        # other?
        else:
            continue
        line = line.split('"')
        
        # position and time
        xcoord = double(line[3])
        ycoord = double(line[1])
        timestamp = double(line[7])
        if clock_type=='COPY':
            x_copy.append(xcoord)
            y_copy.append(ycoord)
            t_copy.append(timestamp)
        elif clock_type=='COMMAND':
            x_command.append(xcoord)
            y_command.append(ycoord)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
    
    f.close()

    x_copy,y_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(t_copy)
    x_command,y_command,t_command = np.array(x_command),np.array(y_command),np.array(t_command)

    # subtract mean values (zm = zero mean)
    x_copy_zm,y_copy_zm = x_copy-mean(x_copy),y_copy-mean(y_copy)
    x_command_zm,y_command_zm = x_command-mean(x_command),y_command-mean(y_command)
    
    # dft's
    dft_size_copy,dft_size_command = len(x_copy),len(x_command)
    k_copy,k_command = np.arange(dft_size_copy),np.arange(dft_size_command)
    dftx_copy,dfty_copy = np.fft.fft(x_copy_zm,n=dft_size_copy),np.fft.fft(y_copy_zm,n=dft_size_copy)
    dftx_command,dfty_command = np.fft.fft(x_command_zm,n=dft_size_command),np.fft.fft(y_command_zm,n=dft_size_command)
    
    # k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than or equal to pi
    k_near_pi_copy,k_near_pi_command = 0,0
    if dft_size_copy%2==0:
        k_near_pi_copy = dft_size_copy/2+1
    else:
        k_near_pi_copy = math.ceil(dft_size_copy/2)
    if dft_size_command%2==0:
        k_near_pi_command = dft_size_command/2+1
    else:
        k_near_pi_command = math.ceil(dft_size_command/2)
    
    # center the coefficients around w = 0
    k_centered_copy = np.linspace(-dft_size_copy/2,dft_size_copy/2,dft_size_copy)
    dftx_centered_copy = np.concatenate((dftx_copy[k_near_pi_copy:],dftx_copy[:k_near_pi_copy]))
    dfty_centered_copy = np.concatenate((dfty_copy[k_near_pi_copy:],dfty_copy[:k_near_pi_copy]))
    
    k_centered_command = np.linspace(-dft_size_command/2,dft_size_copy/2,dft_size_command)
    dftx_centered_command = np.concatenate((dftx_command[k_near_pi_command:],dftx_command[:k_near_pi_command]))
    dfty_centered_command = np.concatenate((dfty_command[k_near_pi_command:],dfty_command[:k_near_pi_command]))

    # percent energy in peak
    # copy clocks
    Ex_copy,Ey_copy = np.abs(dftx_centered_copy)**2,np.abs(dfty_centered_copy)**2
    Ex_total_copy,Ey_total_copy = sum(Ex_copy),sum(Ey_copy)
    Ex_peak_copy,Ey_peak_copy = 2*max(Ex_copy)/Ex_total_copy,2*max(Ey_copy)/Ey_total_copy
    Epeak_x_copy.append((ftype,Ex_peak_copy))
    Epeak_y_copy.append((ftype,Ey_peak_copy))

    # command clocks
    Ex_command,Ey_command = np.abs(dftx_centered_command)**2,np.abs(dfty_centered_command)**2
    Ex_total_command,Ey_total_command = sum(Ex_command),sum(Ey_command)
    Ex_peak_command,Ey_peak_command = 2*max(Ex_command)/Ex_total_command,2*max(Ey_command)/Ey_total_command
    Epeak_x_command.append((ftype,Ex_peak_command))
    Epeak_y_command.append((ftype,Ey_peak_command))

    # percent energy within 1 std. deviation of the center of the energy distribution
    # copy clocks
    # percent energy within 1 std. deviation of the center of the energy distribution
    Ex_var_copy,Ey_var_copy = 0,0
    for w in range(dft_size_copy):
        Ex_var_copy += k_centered_copy[w]**2*Ex_copy[w]/Ex_total_copy
        Ey_var_copy += k_centered_copy[w]**2*Ey_copy[w]/Ey_total_copy
    Ex_std_copy,Ey_std_copy = math.sqrt(Ex_var_copy),math.sqrt(Ey_var_copy)
    Estd_x_copy.append((ftype,Ex_std_copy))
    Estd_y_copy.append((ftype,Ey_std_copy))
    
    Ex_central_copy,Ey_central_copy = 0,0
    for d in range(dft_size_copy):
        if abs((dft_size_copy-1)/2-d)<=Ex_std_copy:
            Ex_central_copy += Ex_copy[d]/Ex_total_copy
        if abs((dft_size_copy-1)/2-d)<=Ey_std_copy:
            Ey_central_copy += Ey_copy[d]/Ey_total_copy
    Ecentral_x_copy.append((ftype,Ex_central_copy))
    Ecentral_y_copy.append((ftype,Ey_central_copy))

    # command clocks
    Ex_var_command,Ey_var_command = 0,0
    for w in range(dft_size_command):
        Ex_var_command += k_centered_command[w]**2*Ex_command[w]/Ex_total_command
        Ey_var_command += k_centered_command[w]**2*Ey_command[w]/Ey_total_command
    Ex_std_command,Ey_std_command = math.sqrt(Ex_var_command),math.sqrt(Ey_var_command)
    Estd_x_command.append((ftype,Ex_std_command))
    Estd_y_command.append((ftype,Ey_std_command))
    
    Ex_central_command,Ey_central_command = 0,0
    for d in range(dft_size_command):
        if abs((dft_size_command-1)/2-d)<=Ex_std_command:
            Ex_central_command += Ex_command[d]/Ex_total_command
        if abs((dft_size_command-1)/2-d)<=Ey_std_command:
            Ey_central_command += Ey_command[d]/Ey_total_command
    Ecentral_x_command.append((ftype,Ex_central_command))
    Ecentral_y_command.append((ftype,Ey_central_command))

# compare energy properties for the drawings of healthy vs. impaired patients
binedges = ct.get_bins(Epeak_x_copy,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_x_copy if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_x_copy if elt[0]=='impaired'],
             binedges,'Fraction of Energy in Largest DFS Coefficient','Epeak_x_copy',path)
binedges = ct.get_bins(Epeak_y_copy,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_y_copy if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y_copy if elt[0]=='impaired'],
             binedges,'Fraction of Energy in Largest DFS Coefficient','Epeak_y_copy',path)
binedges = ct.get_bins(Epeak_x_command,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_x_command if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_x_command if elt[0]=='impaired'],
             binedges,'Fraction of Energy in Largest DFS Coefficient','Epeak_x_command',path)
binedges = ct.get_bins(Epeak_y_command,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_y_command if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y_command if elt[0]=='impaired'],
             binedges,'Fraction of Energy in Largest DFS Coefficient','Epeak_y_command',path)

binedges = ct.get_bins(Estd_x_copy,nbins=10)
ct.make_hist([elt[1] for elt in Estd_x_copy if elt[0]=='healthy'],
             [elt[1] for elt in Estd_x_copy if elt[0]=='impaired'],
             binedges,'Std. Deviation of Energy Distribution','Estd_x_copy',path)
binedges = ct.get_bins(Estd_y_copy,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_y_copy if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y_copy if elt[0]=='impaired'],
             binedges,'Std. Deviation of Energy Distribution','Estd_y_copy',path)
binedges = ct.get_bins(Estd_x_command,nbins=10)
ct.make_hist([elt[1] for elt in Estd_x_command if elt[0]=='healthy'],
             [elt[1] for elt in Estd_x_command if elt[0]=='impaired'],
             binedges,'Std. Deviation of Energy Distribution','Estd_x_command',path)
binedges = ct.get_bins(Estd_y_copy,nbins=10)
ct.make_hist([elt[1] for elt in Epeak_y_command if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y_command if elt[0]=='impaired'],
             binedges,'Std. Deviation of Energy Distribution','Estd_y_command',path)

binedges = ct.get_bins(Ecentral_x_copy,nbins=10)
ct.make_hist([elt[1] for elt in Ecentral_x_copy if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_x_copy if elt[0]=='impaired'],
             binedges,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_x_copy',path)
binedges = ct.get_bins(Ecentral_y_copy,nbins=10)
ct.make_hist([elt[1] for elt in Ecentral_y_copy if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_y_copy if elt[0]=='impaired'],
             binedges,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_y_copy',path)
binedges = ct.get_bins(Ecentral_x_command,nbins=10)
ct.make_hist([elt[1] for elt in Ecentral_x_command if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_x_command if elt[0]=='impaired'],
             binedges,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_x_command',path)
binedges = ct.get_bins(Ecentral_y_command,nbins=10)
ct.make_hist([elt[1] for elt in Ecentral_y_command if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_y_command if elt[0]=='impaired'],
             binedges,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_y_command',path)
