#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""
Data_Sort_and_Plot_Script.py
Script created by Jeffrey N. Murphy (2015); Email: jnmurphy@ualberta.ca
Provided without any warranty.
The directory locations will need to be modified, depending on the location of the input file.
The input for this script is the file produced by Data_Amalgamation_Script.py
"""



Save_Directory = "C:\\Users\\Jeffrey\\Plot Files" #Where plot files iwll be saved

CSV_file = "C:\\Users\\Jeffrey\\Summary Files\\Summary_2.csv"
#CSV_file = "/home/jeffrey/Summary Files/Summary_2.csv" #tux
CSV_directory = CSV_file[0:CSV_file.find("Summary_2.csv")]

output_filename = "Summary_5.csv"





import os
windows = "\\"
linux = "//"
filesep = os.path.sep

import csv
import numpy as np

dataset = []
dataset_index = -1

labels_old = [
        "Image_Title.String","Version","Defect_Density_um",
        "Total_Area_nm","correlation_length","opa_Hermans",
        "LER_width_avg","LER_sigma_avg","wfft_Period_px",
        "wfft_Period_nm","nm_per_pixel"
    ]

labels = ["Image_Number", "Image_Title.String",   "Output_Folder",    "Date", "Time", "Version",  "nm_per_pixel", "Width_initial",    "Height_initial",   "Crop_1",   "wfft_Period_px",   "wfft_Period_nm",   "Smoothing",    "Smoothing_radius", "Crop_2",   "Width_final",  "Height_final", "Threshold.String", "Threshold",    "Threshold.Auto-Local.String",  "Up_Thresh",    "Low_Thresh",   "PA.pos_nPixels.1", "PA.pos_mean.1",    "PA.pos_min.1", "PA.pos_max.1", "PA.pos_std.1", "PA.neg_nPixels.1", "PA.neg_mean.1",    "PA.neg_min.1", "PA.neg_max.1", "PA.neg_std.1", "PA.nPositive.1",   "PA.nNegative.1",   "PA.nTotal.1",  "PA.pos_nPixels.2", "PA.pos_mean.2",    "PA.pos_min.2", "PA.pos_max.2", "PA.pos_std.2", "PA.neg_nPixels.2", "PA.neg_mean.2",    "PA.neg_min.2", "PA.neg_max.2", "PA.neg_std.2", "PA.nPositive.2",   "PA.nNegative.2",   "PA.nTotal.2",  "PA.SET.large_drops",   "PA.SET.min_area",  "PA.SET.wlsq_iterations",   "PA.WLSQ.positive.0",   "PA.WLSQ.positive.1",   "PA.WLSQ.positive.2",   "PA.WLSQ.positive.3",   "PA.WLSQ.positive.4",   "PA.WLSQ.positive.5",   "PA.WLSQ.positive.6",   "PA.WLSQ.negative.0",   "PA.WLSQ.negative.1",   "PA.WLSQ.negative.2",   "PA.WLSQ.negative.3",   "PA.WLSQ.negative.4",   "PA.WLSQ.negative.5",   "PA.WLSQ.negative.6",   "PA.Width.positive",    "PA.Width.negative",    "PA.Width.proportion",  "PA.Blobs.p.i.count",   "PA.Blobs.p.i.area",    "PA.Blobs.p.f.area",    "PA.Blobs.p.f.count",   "PA.Blobs.n.i.count",   "PA.Blobs.n.i.area",    "PA.Blobs.n.f.area",    "PA.Blobs.n.f.count",   "dot_max_area_pos", "pos_edge_dot_count",   "pos_dot_count",    "dot_max_area_neg", "neg_edge_dot_count",   "neg_dot_count",    "line_min_area_pos",    "pos_edge_line_count",  "pos_line_count",   "line_min_area_neg",    "neg_edge_line_count",  "neg_line_count",   "pos_min_area", "neg_min_area", "pos_mDist",    "neg_mDist",    "pos_t_defects.L",  "pos_j_defects.L",  "neg_t_defects.L",  "neg_j_defects.L",  "Skel.Coverage.Metric.pos", "Skel.Coverage.Metric.neg", "Correlation.Phase",    "opa_factor",   "opa_set_ds",   "opa_set_ds_points",    "opa_Hermans",  "Array_Length_initial", "Array_Length_downsampled", "correlation_length",   "correlation_length_linear",    "opa_bar_R_squared",    "loop_count",   "P.LER_sigma_avg",  "N.LER_sigma_avg",  "LER_length_total_px",  "LER_N_width_avg",  "LER_P_width_avg",  "LER_width_avg",    "P.E.sigma.avg",    "P.E.sigma.sigma",  "P.E.sigma.min",    "P.E.sigma.max",    "P.E.sigma.count",  "P.E.avg.avg",  "P.E.avg.sigma",    "P.E.avg.min",  "P.E.avg.max",  "P.E.avg.count",    "P.E.count.avg",    "P.E.count.sigma",  "P.E.count.min",    "P.E.count.max",    "P.E.count.count",  "P.E.sum.avg",  "P.E.sum.sigma",    "P.E.sum.min",  "P.E.sum.max",  "P.E.sum.count",    "P.E.min.avg",  "P.E.min.sigma",    "P.E.min.min",  "P.E.min.max",  "P.E.min.count",    "P.E.max.avg",  "P.E.max.sigma",    "P.E.max.min",  "P.E.max.max",  "P.E.max.count",    "P.W.sigma.avg",    "P.W.sigma.sigma",  "P.W.sigma.min",    "P.W.sigma.max",    "P.W.sigma.count",  "P.W.avg.avg",  "P.W.avg.sigma",    "P.W.avg.min",  "P.W.avg.max",  "P.W.avg.count",    "P.W.count.avg",    "P.W.count.sigma",  "P.W.count.min",    "P.W.count.max",    "P.W.count.count",  "P.W.sum.avg",  "P.W.sum.sigma",    "P.W.sum.min",  "P.W.sum.max",  "P.W.sum.count",    "P.W.min.avg",  "P.W.min.sigma",    "P.W.min.min",  "P.W.min.max",  "P.W.min.count",    "P.W.max.avg",  "P.W.max.sigma",    "P.W.max.min",  "P.W.max.max",  "P.W.max.count",    "N.E.sigma.avg",    "N.E.sigma.sigma",  "N.E.sigma.min",    "N.E.sigma.max",    "N.E.sigma.count",  "N.E.avg.avg",  "N.E.avg.sigma",    "N.E.avg.min",  "N.E.avg.max",  "N.E.avg.count",    "N.E.count.avg",    "N.E.count.sigma",  "N.E.count.min",    "N.E.count.max",    "N.E.count.count",  "N.E.sum.avg",  "N.E.sum.sigma",    "N.E.sum.min",  "N.E.sum.max",  "N.E.sum.count",    "N.E.min.avg",  "N.E.min.sigma",    "N.E.min.min",  "N.E.min.max",  "N.E.min.count",    "N.E.max.avg",  "N.E.max.sigma",    "N.E.max.min",  "N.E.max.max",  "N.E.max.count",    "N.W.sigma.avg",    "N.W.sigma.sigma",  "N.W.sigma.min",    "N.W.sigma.max",    "N.W.sigma.count",  "N.W.avg.avg",  "N.W.avg.sigma",    "N.W.avg.min",  "N.W.avg.max",  "N.W.avg.count",    "N.W.count.avg",    "N.W.count.sigma",  "N.W.count.min",    "N.W.count.max",    "N.W.count.count",  "N.W.sum.avg",  "N.W.sum.sigma",    "N.W.sum.min",  "N.W.sum.max",  "N.W.sum.count",    "N.W.min.avg",  "N.W.min.sigma",    "N.W.min.min",  "N.W.min.max",  "N.W.min.count",    "N.W.max.avg",  "N.W.max.sigma",    "N.W.max.min",  "N.W.max.max",  "N.W.max.count",    "V.E.sigma.avg",    "V.E.sigma.sigma",  "V.E.sigma.min",    "V.E.sigma.max",    "V.E.sigma.count",  "V.E.avg.avg",  "V.E.avg.sigma",    "V.E.avg.min",  "V.E.avg.max",  "V.E.avg.count",    "V.E.count.avg",    "V.E.count.sigma",  "V.E.count.min",    "V.E.count.max",    "V.E.count.count",  "V.E.sum.avg",  "V.E.sum.sigma",    "V.E.sum.min",  "V.E.sum.max",  "V.E.sum.count",    "V.E.min.avg",  "V.E.min.sigma",    "V.E.min.min",  "V.E.min.max",  "V.E.min.count",    "V.E.max.avg",  "V.E.max.sigma",    "V.E.max.min",  "V.E.max.max",  "V.E.max.count",    "V.W.sigma.avg",    "V.W.sigma.sigma",  "V.W.sigma.min",    "V.W.sigma.max",    "V.W.sigma.count",  "V.W.avg.avg",  "V.W.avg.sigma",    "V.W.avg.min",  "V.W.avg.max",  "V.W.avg.count",    "V.W.count.avg",    "V.W.count.sigma",  "V.W.count.min",    "V.W.count.max",    "V.W.count.count",  "V.W.sum.avg",  "V.W.sum.sigma",    "V.W.sum.min",  "V.W.sum.max",  "V.W.sum.count",    "V.W.min.avg",  "V.W.min.sigma",    "V.W.min.min",  "V.W.min.max",  "V.W.min.count",    "V.W.max.avg",  "V.W.max.sigma",    "V.W.max.min",  "V.W.max.max",  "V.W.max.count",    "Skel.Dist.avg",    "Skel.Dist.sigma",  "Skel.Dist.min",    "Skel.Dist.max",    "Skel.Dist.count",  "DRAW_edge_limit",  "DRAW_jn_radius",   "edge_limit",   "PTE",  "NTE",  "PT",   "NT",   "PDE",  "NDE",  "PD",   "ND",   "PJ3",  "NJ3",  "PJ4",  "NJ4",  "PJx",  "NJx",  "Ptot", "Ntot", "Total_Defects",    "Total_Area_px",    "Total_Area_nm",    "Defect_Density_nm",    "Defect_Density_um"]


def clean_number(s):
    if type(s) != type(''):
        return True
    else:
        if s[0] == '-' and len(s)>1:
            s = s[1:]
        s = s.replace('"','')
        return s.replace('.','',1).isdigit()

def clean_my_data(x_series,y_series):
    x_new , y_new = [] , []
    for i in range(0,len(x_series)):
        if type(x_series[i]) != type('') and type(y_series[i]) != type(''):
            x_new.append(x_series[i])
            y_new.append(y_series[i])
    return x_new , y_new

def clean_my_data_3(x_series,y_series,z_series):
    x_new , y_new , z_new = [] , [] , []
    for i in range(0,len(x_series)):
        if type(x_series[i]) != type('') and type(y_series[i]) != type('') and type(z_series[i]) != type(''):
            x_new.append(x_series[i])
            y_new.append(y_series[i])
            z_new.append(z_series[i])
    return x_new , y_new , z_new

def clean_my_data_4(x_series,y_series,z_series, a_series):
    x_new , y_new , z_new , a_new = [] , [] , [] , []
    for i in range(0,len(x_series)):
        if type(x_series[i]) != type('') and type(y_series[i]) != type('') and type(z_series[i]) != type('') and type(a_series[i]) != type(''):
            x_new.append(x_series[i])
            y_new.append(y_series[i])
            z_new.append(z_series[i])
            a_new.append(a_series[i])
    return x_new , y_new , z_new , a_new

def clean_my_data_5(x_series,y_series,z_series, a_series, b_series):
    x_new , y_new , z_new , a_new , b_new = [] , [] , [] , [] , []
    for i in range(0,len(x_series)):
        if type(x_series[i]) != type('') and type(y_series[i]) != type('') and type(z_series[i]) != type('') and type(a_series[i]) != type('') and type(b_series[i]) != type(''):
            x_new.append(x_series[i])
            y_new.append(y_series[i])
            z_new.append(z_series[i])
            a_new.append(a_series[i])
            b_new.append(b_series[i])
    return x_new , y_new , z_new , a_new , b_new

rfile = open(CSV_file,"rb")
reader = csv.reader(rfile)

temporary_data = []
for row in reader:
    temp = row[0].split("\t")
    for i in range(0,len(temp)):
        temp[i] = temp[i].replace('"""', '"')
        temp[i] = temp[i].replace('""', '"')
        if clean_number(temp[i]):
            temp[i] = float(temp[i].replace('"',''))
    #print(clean_number(temp[1]))
    #print(temp[1])
    #print("\n")
    temporary_data.append(temp)


rfile.close()
ofile  = open(CSV_directory + output_filename, "wb")
writer = csv.writer(ofile, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL) #NONNUMERIC) #quoting=csv.QUOTE_NONE)


for i in range(0,len(temporary_data)):
    data = []
    for j in range(0,len(temporary_data[i])):
        data.append(temporary_data[i][j])
    #print("\n")
    writer.writerow(temporary_data[i])
    #writer.writerow(data)
    
ofile.close() 

print("DATA EXPORTED")

"""
Assembling Data 
"""

''' Get BCP List '''
bcp_set = []
bcp_legend = []
bcp_row = []
bcp_datasets = []
data_labels = labels #["nm_per_pixel", "wfft_Period_px", "wfft_Period_nm", "V.W.sigma.avg", "V.W.avg.avg","V.E.sigma.avg","V.E.avg.avg","V.W.avg.sigma"]
data_label_index = []


for i in range(1,len(temporary_data)):
    # 1: skip the first line with column headings
    bcp = temporary_data[i][0].split(' ')
    if bcp[0] not in bcp_set:
        bcp_set.append(bcp[0])
        bcp_datasets.append([])
        bcp_row.append([])
        #for k in range(0,len(data_labels)):
        #    bcp_datasets[len(bcp_datasets)-1].append([])

bcp_set.sort()

for i in range(0,len(bcp_set)):
    bcp_legend.append(bcp_set[i].replace("-",".").replace("_","-"))



print(bcp_set) 
print(bcp_datasets)

''' Find all of the relevant rows for each BCP set
    No need to re-confirm that it matches the BCP  '''
for i in range(0,len(bcp_set)):
    values = []
    for row in range(1,len(temporary_data)):
        # 1: skip the first line with column headings
        bcp = temporary_data[row][0].split(' ')
        if bcp[0] == bcp_set[i]:
            values.append(row)
    bcp_row[i] = values
    print(bcp_row[i])

'''Create data_label_index to avoid need to match each column'''
for k in range(0,len(data_labels)):
    #print(data_labels[k])
    for i in range(0,len(temporary_data[0])):
        #print(temporary_data[0][i])
        if data_labels[k] == temporary_data[0][i].replace('"',''):
            data_label_index.append(i)

print(data_label_index)

## This would be going across the columns... why not go down the columns?
# for i in range(1,len(temporary_data)):
#     # 1: skip the first line with column headings
#     bcp = temporary_data[i][0].split(' ')
#     for j in range(0,len(bcp_set)):
#         if bcp[0] == bcp_set[j]:
#             J = j

"""Collect Values and Attach as a list"""
for i in range(0,len(bcp_set)):
    #print(bcp_set[i])
    for m in range(0,len(data_labels)):
        values = []
        for k in range(0,len(bcp_row[i])):
            values.append(temporary_data[bcp_row[i][k]][data_label_index[m]])
        bcp_datasets[i].append(values)
        #print(data_labels[m])
        #print(values)






"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PLOTS BELOW HERE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
"""






"""
Plotting Data 
"""

import matplotlib.pyplot as plt
# http://matplotlib.org/api/pyplot_api.html

marker_style = ["o","v","s","p","D"] #A list of several possible data markers
marker_color = ["SteelBlue","SeaGreen","GoldenRod","OrangeRed","FireBrick"]
marker_size = [9,10,9,11,9]

#axis_font = {'fontname':'Droid Sans', 'size':'18', 'color':'black', 'weight':'normal'} 

font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 27.5,
        'sans-serif' : 'Arial'}
plt.rc('font', **font)
axes_style = {}
#http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
#http://matplotlib.org/users/customizing.html

# set tick width
plt.rcParams['axes.linewidth'] = 4
#X-axis
plt.rcParams['xtick.major.size'] = 12
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.minor.size'] = 6
plt.rcParams['xtick.minor.width'] = 2
#Y-axis
plt.rcParams['ytick.major.size'] = 12
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['ytick.minor.size'] = 6
plt.rcParams['ytick.minor.width'] = 2
#http://stackoverflow.com/questions/14705904/matplotlib-ticks-thickness


plotlegend = 0 #Set to 1 in order for legends to be turned on.



"""
PLOT 1:

X: wfft_Period_nm
Y: V.E.sigma.avg / V.W.avg.avg

"""

plot_file_1 = "1_LineEdge_StDevNorm_vs_Period"

x_axis_label = "Period (nm)"
y_axis_label = "Line Edge: St.Dev/Mean.Width (nm/nm)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "wfft_Period_nm"
b_label = "V.E.sigma.avg"
c_label = "V.W.avg.avg"
d_label = "nm_per_pixel"
a_data = []
b_data = []
c_data = []
d_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    x_points , b_points , c_points = clean_my_data_3(a_data[k],b_data[k],c_data[k])
    
    y_points = []
    for i in range(0,len(x_points)):
        x_points[i] = float(x_points[i])
        y_points.append( float(b_points[i]) / float(c_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)

    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.ylim([0,0.25])

plt.savefig(CSV_directory + plot_file_1 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_1 + '.svg')

#plt.show()
plt.cla()
plt.clf()

"""BOX PLOT"""

#convert x values to displacements
bp_scale_x = 5
for i in range(0,len(x_boxplot)):
    x_avg = sum(x_boxplot[i])/len(x_boxplot[i])
    print(x_avg)
    for j in range(0,len(x_boxplot[i])):
        x_boxplot[i][j] = (x_boxplot[i][j] - x_avg)/bp_scale_x + i + 1

#print(y_boxplot)
#print(x_boxplot)

plt.boxplot(y_boxplot)
plt.ylim([0,0.25])

for k in range(0,len(y_boxplot)):
    plt.plot(x_boxplot[k],y_boxplot[k], marker_style[k], linestyle="none", markersize=marker_size[k]/2, label=bcp_legend[k], color=marker_color[k], alpha=0.75)

"""Boxplot Labels, etc..."""
if plotlegend: plt.legend(loc='best', numpoints=1)
x_axis_label_bp = "Grouped by Polymer"
plt.xlabel(x_axis_label_bp, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

plt.savefig(CSV_directory + "BP" + plot_file_1 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + "BP" + plot_file_1 + '.svg')
#plt.show()
plt.cla()
plt.clf()






"""
PLOT 2:

X: wfft_Period_nm
Y: V.W.sigma.avg / V.W.avg.avg

"""

plot_file_2 = "2_LineWidth_StDevNorm_vs_Period"

x_axis_label = "Period (nm)"
y_axis_label = "Line Width: St.Dev/Mean (nm/nm)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "wfft_Period_nm"
b_label = "V.W.sigma.avg"
c_label = "V.W.avg.avg"
d_label = "nm_per_pixel"
a_data = []
b_data = []
c_data = []
d_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)

y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    x_points , b_points , c_points = clean_my_data_3(a_data[k],b_data[k],c_data[k])
    
    y_points = []
    for i in range(0,len(x_points)):
        x_points[i] = float(x_points[i])
        y_points.append( float(b_points[i]) / float(c_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.ylim([0,0.25])

plt.savefig(CSV_directory + plot_file_2 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_2 + '.svg')

#plt.show()
plt.cla()
plt.clf()

"""BOX PLOT"""

#convert x values to displacements
bp_scale_x = 5
for i in range(0,len(x_boxplot)):
    x_avg = sum(x_boxplot[i])/len(x_boxplot[i])
    print(x_avg)
    for j in range(0,len(x_boxplot[i])):
        x_boxplot[i][j] = (x_boxplot[i][j] - x_avg)/bp_scale_x + i + 1

#print(y_boxplot)
#print(x_boxplot)

plt.boxplot(y_boxplot)
plt.ylim([0,0.25])

for k in range(0,len(y_boxplot)):
    plt.plot(x_boxplot[k],y_boxplot[k], marker_style[k], linestyle="none", markersize=marker_size[k]/2, label=bcp_legend[k], color=marker_color[k], alpha=0.75)

"""Boxplot Labels, etc..."""
if plotlegend: plt.legend(loc='best', numpoints=1)
x_axis_label_bp = "Grouped by Polymer"
plt.xlabel(x_axis_label_bp, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

plt.savefig(CSV_directory + "BP" + plot_file_2 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + "BP" + plot_file_2 + '.svg')
#plt.show()
plt.cla()
plt.clf()







"""
PLOT 3:

X: wfft_Period_nm
Y: V.W.avg.avg +/- V.W.avg.sigma

"""

plot_file_3 = "3_LineWidth_vs_Period"

x_axis_label = "Period (nm)"
y_axis_label = "Line Width (nm)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "wfft_Period_nm"
b_label = "V.W.avg.avg"
c_label = "V.W.avg.sigma"
d_label = "nm_per_pixel"
a_data = []
b_data = []
c_data = []
d_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    x_points , b_points , c_points, d_points = clean_my_data_4(a_data[k],b_data[k],c_data[k],d_data[k])
    
    y_points = []
    yerr_points = []
    for i in range(0,len(x_points)):
        x_points[i] = float(x_points[i])
        y_points.append( float(b_points[i]) * float(d_points[i]) )
        yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_points,y_points,yerr=yerr_points, linestyle="none", color=marker_color[k])
#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.ylim([0,20])

plt.savefig(CSV_directory + plot_file_3 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_3 + '.svg')

#plt.show()
plt.cla()
plt.clf()


"""BOX PLOT"""

#convert x values to displacements
bp_scale_x = 5
for i in range(0,len(x_boxplot)):
    x_avg = sum(x_boxplot[i])/len(x_boxplot[i])
    print(x_avg)
    for j in range(0,len(x_boxplot[i])):
        x_boxplot[i][j] = (x_boxplot[i][j] - x_avg)/bp_scale_x + i + 1

#print(y_boxplot)
#print(x_boxplot)

plt.boxplot(y_boxplot)
plt.ylim([0,20])

for k in range(0,len(y_boxplot)):
    plt.plot(x_boxplot[k],y_boxplot[k], marker_style[k], linestyle="none", markersize=marker_size[k]/2, label=bcp_legend[k], color=marker_color[k], alpha=0.75)

"""Boxplot Labels, etc..."""
if plotlegend: plt.legend(loc='best', numpoints=1)
x_axis_label_bp = "Grouped by Polymer"
plt.xlabel(x_axis_label_bp, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

plt.savefig(CSV_directory + "BP" + plot_file_3 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + "BP" + plot_file_3 + '.svg')
#plt.show()
plt.cla()
plt.clf()









"""
PLOT 4:

X: nm_per_pixel
Y: opa_Hermans

"""

plot_file_4 = "4_Hermans_vs_Resolution"

x_axis_label = "Resolution (nm/pixel)"
y_axis_label = "Herman's 2D Order Parameter (unitless)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "opa_Hermans"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "wfft_Period_px"
a_data = []
b_data = []
c_data = []
d_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    y_points , b_points, c_points = clean_my_data_3(a_data[k],b_data[k], c_data[k])
    
    x_points = []
    # yerr_points = []
    for i in range(0,len(y_points)):
        #x_points[i] = float(x_points[i])
        x_points.append( float(b_points[i]) / float(c_points[i]) )
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.ylim([0,1])

plt.savefig(CSV_directory + plot_file_4 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_4 + '.svg')

#plt.show()
plt.cla()
plt.clf()










"""
PLOT 5:

X: total wire length measured
Y: opa_Hermans

"""

plot_file_5 = "5_Hermans_vs_Total_Length_nm"

x_axis_label = "Total Length of Lines (um)"
y_axis_label = "Herman's 2D Order Parameter (unitless)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "opa_Hermans"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    y_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(y_points)):
        #x_points[i] = float(x_points[i])
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / float(c_points[i]) /1000 )
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.ylim([0,1])
plt.xscale('log')
plt.xlim([1,10000])

plt.savefig(CSV_directory + plot_file_5 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_5 + '.svg')

#plt.show()
plt.cla()
plt.clf()




"""
PLOT 6:

X: total wire length measured
Y: correlation_length

"""

plot_file_6 = "6_Correlation_nm_vs_Total_Length_nm"

x_axis_label = "Total Length of Lines (um)"
y_axis_label = "Correlation Length (nm)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "correlation_length"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    a_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    y_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(a_points)):
        #x_points[i] = float(x_points[i])
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / float(c_points[i]) /1000 )
        y_points.append(float(a_points[i])*float(b_points[i]))
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.yscale('log')
plt.ylim([10,2000])
plt.xscale('log')
plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_6 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_6 + '.svg')

#plt.show()
plt.cla()
plt.clf()







"""
PLOT 7:

X: Image Area (um^2)
Y: correlation_length

"""

plot_file_7 = "7_Correlation_nm_vs_ImageArea_um2"

x_axis_label = "Image Area (um )"  #insert superscript-2 later
y_axis_label = "Correlation Length (nm)"
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "correlation_length"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    a_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    y_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(a_points)):
        #x_points[i] = float(x_points[i])
        y_points.append(float(a_points[i])*float(b_points[i]))
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 )
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.yscale('log')
plt.ylim([10,2000])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_7 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_7 + '.svg')

#plt.show()
plt.cla()
plt.clf()











"""
PLOT 8:

X: Image Area (um^2)
Y: Defect Density

"""

plot_file_8 = "8_DefectDensity_vs_ImageArea_um2"

x_axis_label = "Image Area (um )"  #insert superscript 2 later
y_axis_label = "Defect Pair Density (um  )" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "Defect_Density_um"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
for k in range(0,len(bcp_set)):
    y_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(y_points)):
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        print(stagger)
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)

if plotlegend: plt.legend(loc='best', numpoints=1)
plt.yscale('log')
plt.ylim([10,3000])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_8 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_8 + '.svg')

#plt.show()
plt.cla()
plt.clf()









"""
PLOT 9: (Averaged)

X: Image Area (um^2)
Y: Defect Density

"""

plot_file_9 = "9_DefectDensity_vs_ImageArea_um2"

x_axis_label = "Image Area (um )"  #insert superscript 2 later
y_axis_label = "Defect Pair Density (um  )" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "Defect_Density_um"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
yerr_points = []
for k in range(0,len(bcp_set)):
    y_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(y_points)):
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        #print(stagger)
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    # GROUP & AVERAGE
    x_new = []
    for i in range(0,len(x_points)):
        count = 0
        for j in range(0,len(x_points)):
            if x_points[i] == x_points[j]:
                count += 1
        if x_points[i] not in x_new and count > 1:
            x_new.append(x_points[i])
    y_new = []
    yerr_new_neg = [] #for log plots
    yerr_new = []
    for j in range(0,len(x_new)):
        y_temp = []
        for i in range(0,len(x_points)):
            if x_points[i] == x_new[j]:
                y_temp.append(y_points[i])
        y_new.append(np.mean(y_temp))
        yerr_new.append(np.std(y_temp))
        if np.mean(y_temp) - np.std(y_temp) <= 0:
            y_neg_std = np.mean(y_temp) - 0.01
        else:
            y_neg_std = np.std(y_temp)
        yerr_new_neg.append(y_neg_std)
    x_points = x_new
    y_points = y_new
    yerr_points.append(yerr_new)

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_points,y_points,yerr=[yerr_new_neg, yerr_new], linestyle="none", color=marker_color[k])
    print(bcp_legend[k])
    print(x_points)
    print(yerr_points[k])

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)


if plotlegend: plt.legend(loc='best', numpoints=1)
plt.yscale('log')
plt.ylim([10,3000])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_9 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_9 + '.svg')

#plt.show()
plt.cla()
plt.clf()















"""
PLOT 10: (Averaged)

X: Image Area (um^2)
Y: Defect Density

"""

plot_file_10 = "10_DefectDensity_vs_ImageArea_um2"

x_axis_label = "Image Area (um )"  #insert superscript 2 later
y_axis_label = "Defect Pair Density (um  )" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "Defect_Density_um"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
yerr_points = []
for k in range(0,len(bcp_set)):
    y_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(y_points)):
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        #print(stagger)
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    # GROUP & AVERAGE
    x_new = []
    for i in range(0,len(x_points)):
        count = 0
        for j in range(0,len(x_points)):
            if x_points[i] == x_points[j]:
                count += 1
        if x_points[i] not in x_new and count > 1:
            x_new.append(x_points[i])
    x_new.sort()
    y_new = []
    yerr_new_neg = [] #for log plots
    yerr_new = []
    for j in range(0,len(x_new)):
        y_temp = []
        for i in range(0,len(x_points)):
            if x_points[i] == x_new[j]:
                y_temp.append(y_points[i])
        y_new.append(np.mean(y_temp))
        yerr_new.append(np.std(y_temp))
        if np.mean(y_temp) - np.std(y_temp) <= 0:
            y_neg_std = np.mean(y_temp) - 0.01
        else:
            y_neg_std = np.std(y_temp)
        yerr_new_neg.append(y_neg_std)
    #x_points = x_new
    #y_points = y_new
    yerr_points.append(yerr_new)

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.2)
    plt.plot(x_new, y_new, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_new,y_new,yerr=[yerr_new_neg, yerr_new], linestyle="--", color=marker_color[k])
    print(bcp_legend[k])
    print(x_points)
    print(yerr_points[k])

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)


if plotlegend: plt.legend(loc='top right', numpoints=1, bbox_to_anchor=(1, 1))
plt.yscale('log')
plt.ylim([10,3000])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_10 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_10 + '.svg')

#plt.show()
plt.cla()
plt.clf()





"""
PLOT 11: (Averaged)

X: Image Area (um^2)
Y: Correlation Length

Description: Averaged with y-error bars, overlaid on semi-transparent raw datapoints.

"""

plot_file_11 = "11_CorrelationLength_vs_ImageArea_um2"

x_axis_label = "Image Area (um^2)"  #insert superscript 2 later
y_axis_label = "Correlation Length (nm)" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "correlation_length"
b_label = "nm_per_pixel"
c_label = "wfft_Period_nm"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
yerr_points = []
for k in range(0,len(bcp_set)):
    a_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    y_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(a_points)):
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        #print(stagger)
        y_points.append(float(a_points[i])*float(b_points[i]))
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    # GROUP & AVERAGE
    x_new = []
    for i in range(0,len(x_points)):
        count = 0
        for j in range(0,len(x_points)):
            if x_points[i] == x_points[j]:
                count += 1
        if x_points[i] not in x_new and count > 1:
            x_new.append(x_points[i])
    x_new.sort()
    y_new = []
    yerr_new_neg = [] #for log plots
    yerr_new = []
    for j in range(0,len(x_new)):
        y_temp = []
        for i in range(0,len(x_points)):
            if x_points[i] == x_new[j]:
                y_temp.append(y_points[i])
        y_new.append(np.mean(y_temp))
        yerr_new.append(np.std(y_temp))
        if np.mean(y_temp) - np.std(y_temp) <= 0:
            y_neg_std = np.mean(y_temp) - 0.01
        else:
            y_neg_std = np.std(y_temp)
        yerr_new_neg.append(y_neg_std)
    #x_points = x_new
    #y_points = y_new
    yerr_points.append(yerr_new)

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.2)
    plt.plot(x_new, y_new, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_new,y_new,yerr=[yerr_new_neg, yerr_new], linestyle="--", color=marker_color[k])
    print(bcp_legend[k])
    print(x_points)
    print(yerr_points[k])

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)


if plotlegend: plt.legend(loc='top right', numpoints=1, bbox_to_anchor=(1, 1))
plt.yscale('log')
plt.ylim([10,3000])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_11 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_11 + '.svg')

#plt.show()
plt.cla()
plt.clf()













"""
PLOT 12: (Averaged)

X: Image Area (um^2)
Y: LER

Description: Averaged with y-error bars, overlaid on semi-transparent raw datapoints.

"""

plot_file_12 = "12_LineEdge_StDevNorm_vs_ImageArea_um2"

x_axis_label = "Image Area (um^2)"  #insert superscript 2 later
y_axis_label = "Line Edge: St.Dev/Mean.Width (nm/nm)" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "V.E.sigma.avg"
b_label = "nm_per_pixel"
c_label = "V.W.avg.avg"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
yerr_points = []
for k in range(0,len(bcp_set)):
    a_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    y_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(a_points)):
        y_points.append( float(a_points[i]) / float(c_points[i]) )
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        #print(stagger)
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    # GROUP & AVERAGE
    x_new = []
    for i in range(0,len(x_points)):
        count = 0
        for j in range(0,len(x_points)):
            if x_points[i] == x_points[j]:
                count += 1
        if x_points[i] not in x_new and count > 1:
            x_new.append(x_points[i])
    x_new.sort()
    y_new = []
    yerr_new_neg = [] #for log plots
    yerr_new = []
    for j in range(0,len(x_new)):
        y_temp = []
        for i in range(0,len(x_points)):
            if x_points[i] == x_new[j]:
                y_temp.append(y_points[i])
        y_new.append(np.mean(y_temp))
        yerr_new.append(np.std(y_temp))
        if np.mean(y_temp) - np.std(y_temp) <= 0:
            y_neg_std = np.mean(y_temp) - 0.01
        else:
            y_neg_std = np.std(y_temp)
        yerr_new_neg.append(y_neg_std)
    #x_points = x_new
    #y_points = y_new
    yerr_points.append(yerr_new)

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.2)
    plt.plot(x_new, y_new, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_new,y_new,yerr=[yerr_new_neg, yerr_new], linestyle="--", color=marker_color[k])
    print(bcp_legend[k])
    print(x_points)
    print(yerr_points[k])

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)


if plotlegend: plt.legend(loc='top right', numpoints=1, bbox_to_anchor=(1, 1))
#plt.yscale('log')
plt.ylim([0.0,0.25])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_12 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_12 + '.svg')

#plt.show()
plt.cla()
plt.clf()








"""
PLOT 13: (Averaged)

X: Image Area (um^2)
Y: LWR/W

Description: Averaged with y-error bars, overlaid on semi-transparent raw datapoints.

"""

plot_file_13 = "13_LineWidth_StDevNorm_vs_ImageArea_um2"

x_axis_label = "Image Area (um^2)"  #insert superscript 2 later
y_axis_label = "Line Edge: St.Dev/Mean.Width (nm/nm)" #insert superscript -2 later
x_points = []
y_points = []

''' Need to ensure that these labels are included in data_labels
    Could just set data_labels = labels, but that slows down re-compiling the data'''
a_label = "V.W.sigma.avg"
b_label = "nm_per_pixel"
c_label = "V.W.avg.avg"
d_label = "Width_initial"
e_label = "Height_initial"
a_data = []
b_data = []
c_data = []
d_data = []
e_data = []

for i in range(0,len(data_labels)):
    if a_label == data_labels[i]:
        i_a = data_label_index[i]
    if b_label == data_labels[i]:
        i_b = data_label_index[i]
    if c_label == data_labels[i]:
        i_c = data_label_index[i]
    if d_label == data_labels[i]:
        i_d = data_label_index[i]
    if e_label == data_labels[i]:
        i_e = data_label_index[i]

print(data_labels)
print("i_a: " + str(i_a))
print("i_b: " + str(i_b))
print("i_c: " + str(i_c))
print("i_d: " + str(i_d))
print("i_e: " + str(i_e))
## put some error in here if it's still -1,-1

#marker_style = ["s","v","o","p","D"] #A list of several possible data markers
#marker_color = ["red","green","blue","cyan","magenta"]
#marker_size = [9,10,9,11,9]
#colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

## Simple case; no clean-up / filtering required
for k in range(0,len(bcp_set)):
    a_values = []
    b_values = []
    c_values = []
    d_values = []
    e_values = []
    for i in range(0,len(bcp_row[k])):
        #print(bcp_row[k])
        #print(ix)
        #print(iy)
        #a_values.append(float(temporary_data[bcp_row[k][i]][i_a]))
        #b_values.append(float(temporary_data[bcp_row[k][i]][i_b]))
        #c_values.append(float(temporary_data[bcp_row[k][i]][i_c]))
        #d_values.append(float(temporary_data[bcp_row[k][i]][i_d]))
        a_values.append(temporary_data[bcp_row[k][i]][i_a])
        b_values.append(temporary_data[bcp_row[k][i]][i_b])
        c_values.append(temporary_data[bcp_row[k][i]][i_c])
        d_values.append(temporary_data[bcp_row[k][i]][i_d])
        e_values.append(temporary_data[bcp_row[k][i]][i_e])
    a_data.append(a_values)
    b_data.append(b_values)
    c_data.append(c_values)
    d_data.append(d_values)
    e_data.append(e_values)


y_boxplot = []
x_boxplot = []
yerr_points = []
for k in range(0,len(bcp_set)):
    a_points , b_points, c_points, d_points, e_points = clean_my_data_5(a_data[k], b_data[k], c_data[k], d_data[k], e_data[k])
    
    x_points = []
    y_points = []
    # a_label = "opa_Hermans"
    # b_label = "nm_per_pixel"
    # c_label = "wfft_Period_nm"
    # d_label = "Width_initial"
    # e_label = "Height_initial"
    for i in range(0,len(a_points)):
        y_points.append( float(a_points[i]) / float(c_points[i]) )
        #x_points[i] = float(x_points[i])
        stagger = (100.0 + 20.0*( float(k+1) - float(len(bcp_set))/2.0 ) / ( float(len(bcp_set))/2.0 ))/100.0
        #print(stagger)
        x_points.append( float(d_points[i]) * float(e_points[i]) * float(b_points[i])* float(b_points[i]) / 1000000 * stagger)
        #yerr_points.append( float(c_points[i]) * float(d_points[i]) )

    # GROUP & AVERAGE
    x_new = []
    for i in range(0,len(x_points)):
        count = 0
        for j in range(0,len(x_points)):
            if x_points[i] == x_points[j]:
                count += 1
        if x_points[i] not in x_new and count > 1:
            x_new.append(x_points[i])
    x_new.sort()
    y_new = []
    yerr_new_neg = [] #for log plots
    yerr_new = []
    for j in range(0,len(x_new)):
        y_temp = []
        for i in range(0,len(x_points)):
            if x_points[i] == x_new[j]:
                y_temp.append(y_points[i])
        y_new.append(np.mean(y_temp))
        yerr_new.append(np.std(y_temp))
        if np.mean(y_temp) - np.std(y_temp) <= 0:
            y_neg_std = np.mean(y_temp) - 0.01
        else:
            y_neg_std = np.std(y_temp)
        yerr_new_neg.append(y_neg_std)
    #x_points = x_new
    #y_points = y_new
    yerr_points.append(yerr_new)

    y_boxplot.append(y_points)
    x_boxplot.append(x_points)
    #print(x_points)
    #print(y_points)
    plt.plot(x_points, y_points, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.2)
    plt.plot(x_new, y_new, marker_style[k], linestyle="none", markersize=marker_size[k], label=bcp_legend[k], color=marker_color[k], alpha=0.75)
    plt.errorbar(x_new,y_new,yerr=[yerr_new_neg, yerr_new], linestyle="--", color=marker_color[k])
    print(bcp_legend[k])
    print(x_points)
    print(yerr_points[k])

#axis labels
axis_font = {'fontname':'Arial', 'size':'28', 'color':'black', 'weight':'bold'} 
plt.xlabel(x_axis_label, **axis_font)
#####plt.ylabel(y_axis_label, **axis_font)


if plotlegend: plt.legend(loc='top right', numpoints=1, bbox_to_anchor=(1, 1))
#plt.yscale('log')
plt.ylim([0.0,0.25])
plt.xscale('log')
#plt.xlim([1,1000])

plt.savefig(CSV_directory + plot_file_13 + '.png', bbox_inches='tight')
plt.savefig(CSV_directory + plot_file_13 + '.svg')

#plt.show()
plt.cla()
plt.clf()
