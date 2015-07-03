#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""
Data_Amalgamation_Script.py
Script created by Jeffrey N. Murphy (2015); Email: jnmurphy@ualberta.ca
Provided without any warranty.
The directory locations will need to be modified, depending on the location of the input folders.
"""


#Directory = '/home/jeffrey/ownCloud/Algorithm Output/Paper Image Set/all2-output/Output_20140514_111653'
#SaveFile = Directory[Directory.rindex("/"):]

#Save_Directory = "D:\SEM Data\20150131\Output_20150131_103030"
Save_Directory = "/home/jeffrey/Data/Image Analysis"

windows = "\\"
linux = "//"
from os.path import sep
winlinux = os.path.sep

List_of_Directories = [Save_Directory]

List_of_Directories = []

from os import listdir
from os.path import isfile
d = listdir(Save_Directory)
for i in range(0,len(d)):
    if not isfile(d[i]):
        List_of_Directories.append(d[i])

print(List_of_Directories)

import os
import csv

dataset = []
dataset_index = -1

labels_old = [
        "Image_Title.String","Version","Defect_Density_um",
        "Total_Area_nm","correlation_length","opa_Hermans",
        "LER_width_avg","LER_sigma_avg","wfft_Period_px",
        "wfft_Period_nm","nm_per_pixel"


    ]

labels = ["Image_Number", "Image_Title.String",   "Output_Folder",    "Date", "Time", "Version",  "nm_per_pixel", "Width_initial",    "Height_initial",   "Crop_1",   "wfft_Period_px",   "wfft_Period_nm",   "Smoothing",    "Smoothing_radius", "Crop_2",   "Width_final",  "Height_final", "Threshold.String", "Threshold",    "Threshold.Auto-Local.String",  "Up_Thresh",    "Low_Thresh",   "PA.pos_nPixels.1", "PA.pos_mean.1",    "PA.pos_min.1", "PA.pos_max.1", "PA.pos_std.1", "PA.neg_nPixels.1", "PA.neg_mean.1",    "PA.neg_min.1", "PA.neg_max.1", "PA.neg_std.1", "PA.nPositive.1",   "PA.nNegative.1",   "PA.nTotal.1",  "PA.pos_nPixels.2", "PA.pos_mean.2",    "PA.pos_min.2", "PA.pos_max.2", "PA.pos_std.2", "PA.neg_nPixels.2", "PA.neg_mean.2",    "PA.neg_min.2", "PA.neg_max.2", "PA.neg_std.2", "PA.nPositive.2",   "PA.nNegative.2",   "PA.nTotal.2",  "PA.SET.large_drops",   "PA.SET.min_area",  "PA.SET.wlsq_iterations",   "PA.WLSQ.positive.0",   "PA.WLSQ.positive.1",   "PA.WLSQ.positive.2",   "PA.WLSQ.positive.3",   "PA.WLSQ.positive.4",   "PA.WLSQ.positive.5",   "PA.WLSQ.positive.6",   "PA.WLSQ.negative.0",   "PA.WLSQ.negative.1",   "PA.WLSQ.negative.2",   "PA.WLSQ.negative.3",   "PA.WLSQ.negative.4",   "PA.WLSQ.negative.5",   "PA.WLSQ.negative.6",   "PA.Width.positive",    "PA.Width.negative",    "PA.Width.proportion",  "PA.Blobs.p.i.count",   "PA.Blobs.p.i.area",    "PA.Blobs.p.f.area",    "PA.Blobs.p.f.count",   "PA.Blobs.n.i.count",   "PA.Blobs.n.i.area",    "PA.Blobs.n.f.area",    "PA.Blobs.n.f.count",   "dot_max_area_pos", "pos_edge_dot_count",   "pos_dot_count",    "dot_max_area_neg", "neg_edge_dot_count",   "neg_dot_count",    "line_min_area_pos",    "pos_edge_line_count",  "pos_line_count",   "line_min_area_neg",    "neg_edge_line_count",  "neg_line_count",   "pos_min_area", "neg_min_area", "pos_mDist",    "neg_mDist",    "pos_t_defects.L",  "pos_j_defects.L",  "neg_t_defects.L",  "neg_j_defects.L",  "Skel.Coverage.Metric.pos", "Skel.Coverage.Metric.neg", "Correlation.Phase",    "opa_factor",   "opa_set_ds",   "opa_set_ds_points",    "opa_Hermans",  "Array_Length_initial", "Array_Length_downsampled", "correlation_length",   "correlation_length_linear",    "opa_bar_R_squared",    "loop_count",   "P.LER_sigma_avg",  "N.LER_sigma_avg",  "LER_length_total_px",  "LER_N_width_avg",  "LER_P_width_avg",  "LER_width_avg",    "P.E.sigma.avg",    "P.E.sigma.sigma",  "P.E.sigma.min",    "P.E.sigma.max",    "P.E.sigma.count",  "P.E.avg.avg",  "P.E.avg.sigma",    "P.E.avg.min",  "P.E.avg.max",  "P.E.avg.count",    "P.E.count.avg",    "P.E.count.sigma",  "P.E.count.min",    "P.E.count.max",    "P.E.count.count",  "P.E.sum.avg",  "P.E.sum.sigma",    "P.E.sum.min",  "P.E.sum.max",  "P.E.sum.count",    "P.E.min.avg",  "P.E.min.sigma",    "P.E.min.min",  "P.E.min.max",  "P.E.min.count",    "P.E.max.avg",  "P.E.max.sigma",    "P.E.max.min",  "P.E.max.max",  "P.E.max.count",    "P.W.sigma.avg",    "P.W.sigma.sigma",  "P.W.sigma.min",    "P.W.sigma.max",    "P.W.sigma.count",  "P.W.avg.avg",  "P.W.avg.sigma",    "P.W.avg.min",  "P.W.avg.max",  "P.W.avg.count",    "P.W.count.avg",    "P.W.count.sigma",  "P.W.count.min",    "P.W.count.max",    "P.W.count.count",  "P.W.sum.avg",  "P.W.sum.sigma",    "P.W.sum.min",  "P.W.sum.max",  "P.W.sum.count",    "P.W.min.avg",  "P.W.min.sigma",    "P.W.min.min",  "P.W.min.max",  "P.W.min.count",    "P.W.max.avg",  "P.W.max.sigma",    "P.W.max.min",  "P.W.max.max",  "P.W.max.count",    "N.E.sigma.avg",    "N.E.sigma.sigma",  "N.E.sigma.min",    "N.E.sigma.max",    "N.E.sigma.count",  "N.E.avg.avg",  "N.E.avg.sigma",    "N.E.avg.min",  "N.E.avg.max",  "N.E.avg.count",    "N.E.count.avg",    "N.E.count.sigma",  "N.E.count.min",    "N.E.count.max",    "N.E.count.count",  "N.E.sum.avg",  "N.E.sum.sigma",    "N.E.sum.min",  "N.E.sum.max",  "N.E.sum.count",    "N.E.min.avg",  "N.E.min.sigma",    "N.E.min.min",  "N.E.min.max",  "N.E.min.count",    "N.E.max.avg",  "N.E.max.sigma",    "N.E.max.min",  "N.E.max.max",  "N.E.max.count",    "N.W.sigma.avg",    "N.W.sigma.sigma",  "N.W.sigma.min",    "N.W.sigma.max",    "N.W.sigma.count",  "N.W.avg.avg",  "N.W.avg.sigma",    "N.W.avg.min",  "N.W.avg.max",  "N.W.avg.count",    "N.W.count.avg",    "N.W.count.sigma",  "N.W.count.min",    "N.W.count.max",    "N.W.count.count",  "N.W.sum.avg",  "N.W.sum.sigma",    "N.W.sum.min",  "N.W.sum.max",  "N.W.sum.count",    "N.W.min.avg",  "N.W.min.sigma",    "N.W.min.min",  "N.W.min.max",  "N.W.min.count",    "N.W.max.avg",  "N.W.max.sigma",    "N.W.max.min",  "N.W.max.max",  "N.W.max.count",    "V.E.sigma.avg",    "V.E.sigma.sigma",  "V.E.sigma.min",    "V.E.sigma.max",    "V.E.sigma.count",  "V.E.avg.avg",  "V.E.avg.sigma",    "V.E.avg.min",  "V.E.avg.max",  "V.E.avg.count",    "V.E.count.avg",    "V.E.count.sigma",  "V.E.count.min",    "V.E.count.max",    "V.E.count.count",  "V.E.sum.avg",  "V.E.sum.sigma",    "V.E.sum.min",  "V.E.sum.max",  "V.E.sum.count",    "V.E.min.avg",  "V.E.min.sigma",    "V.E.min.min",  "V.E.min.max",  "V.E.min.count",    "V.E.max.avg",  "V.E.max.sigma",    "V.E.max.min",  "V.E.max.max",  "V.E.max.count",    "V.W.sigma.avg",    "V.W.sigma.sigma",  "V.W.sigma.min",    "V.W.sigma.max",    "V.W.sigma.count",  "V.W.avg.avg",  "V.W.avg.sigma",    "V.W.avg.min",  "V.W.avg.max",  "V.W.avg.count",    "V.W.count.avg",    "V.W.count.sigma",  "V.W.count.min",    "V.W.count.max",    "V.W.count.count",  "V.W.sum.avg",  "V.W.sum.sigma",    "V.W.sum.min",  "V.W.sum.max",  "V.W.sum.count",    "V.W.min.avg",  "V.W.min.sigma",    "V.W.min.min",  "V.W.min.max",  "V.W.min.count",    "V.W.max.avg",  "V.W.max.sigma",    "V.W.max.min",  "V.W.max.max",  "V.W.max.count",    "Skel.Dist.avg",    "Skel.Dist.sigma",  "Skel.Dist.min",    "Skel.Dist.max",    "Skel.Dist.count",  "DRAW_edge_limit",  "DRAW_jn_radius",   "edge_limit",   "PTE",  "NTE",  "PT",   "NT",   "PDE",  "NDE",  "PD",   "ND",   "PJ3",  "NJ3",  "PJ4",  "NJ4",  "PJx",  "NJx",  "Ptot", "Ntot", "Total_Defects",    "Total_Area_px",    "Total_Area_nm",    "Defect_Density_nm",    "Defect_Density_um"]




for Directory in List_of_Directories:
    print(Directory)
    if "/" in Directory:
        SaveFile = Directory[Directory.rindex(winlinux):]
    else:
        SaveFile = Directory
        Directory = Save_Directory + winlinux + Directory
    items = os.listdir(Directory)

    dirs = [d for d in os.listdir(Directory) if os.path.isdir(os.path.join(Directory, d))]
    # http://stackoverflow.com/questions/7781545/how-to-get-all-folder-only-in-a-given-path-in-python
    print(dirs)
    #dataset = []
    #dataset_index = -1
    

    list_of_folders_missing_output = []

    for d in dirs:
        print(Directory + winlinux + d)

        file_list = os.listdir(Directory + winlinux + d)
        if "outputTD_vertical.xls" in file_list:
            print("found")
            rfile = open(Directory + winlinux + d + winlinux + "outputTD_vertical.xls","rb")
            reader = csv.reader(rfile)

            temporary_data = []
            for row in reader:
                print(row)
                temporary_data.append(row)
            dataset.append({})
            dataset_index += 1
            dataset[dataset_index]["File"] = SaveFile
            for j in range(0,len(labels)):
                dataset[dataset_index][labels[j]] = -1
            for item in labels:
                unfound = True
                for stuff in temporary_data:
                    #print(stuff)
                    stuff_list = stuff[0].split('\t')
                    if item in stuff_list:
                        dataset[dataset_index][item] = stuff_list[2]
                        unfound = False
                if unfound:
                    unfoundx = 0
        else:
            list_of_folders_missing_output.append(Directory + winlinux + d)

    #print(dataset)


ofile  = open(Save_Directory + winlinux + "Summary_2.csv", "wb")
writer = csv.writer(ofile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

# Labels
spreadsheet_labels = ["File"] + labels
writer.writerow(spreadsheet_labels)

for i in range(0,len(dataset)):
    row = []
    row.append(dataset[i]["File"])
    for entry in labels:
        row.append(dataset[i][entry])
    #print(row)
    writer.writerow(row)
ofile.close() 




            
