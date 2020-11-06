#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fastqc Report Parser to output from FastQC file
Generate graphs for each FastQC modules.

"""
import os
import pandas as pd
import numpy as np
import errno #used to raise error in making dire
import matplotlib.pyplot as plt
import seaborn as sns
import argparse 

cwd = os.getcwd()

def directories(folders):
    
    """ Create directory for each fastqc module report. 
        Data is printed to the user command line.  
        Acceptes one Aguments :-
                                folder = List of Directories name for fastqc modules.  """
    
    for i in folders:
        try:
            os.makedirs(i)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
def filter_report_generator(file_path, output_dire, module):
    
    
    """ Read Fastqc file and generate filter report for all fastqc modules. 
          Returns plain txt file with in the specified directory.
          Acceptes one Aguments :-
                                  file_path = File name of fastqc
                                  output_dire = Output Folder | directory name
                                  module = Fastqc module name | filter"""
    
    try:
        # opening fastqc data
        fastqc_data = open(file_path, 'r')
    except OSError as exception:
        print('failed to open file', exception)
    # list to store modules
    data = []
    
    lines = fastqc_data.readlines()

    # reading module header from the fastqc file with for loop and if condition
    for line in lines:
        if line.startswith(">>")  and not line.startswith(">>END"):
            data.append(line)
    
    
    fname = "\t" + file_path
    
    df = [] 
    for d in data:
        d = d.strip()
        d += fname
        df.append(d)
    
    
    meta_data = "Module\t" + "Status\t" + "File name"
    df.insert(0, meta_data)   
    data_join = "\n".join(df)
    
    #creating directory
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    # creating file path for both image and text file output
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "Filter", 'Filter_summary.txt')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Filter_summary.txt')
            
    #write to file 
    fw = open(output_path_dire_txt, 'w' )
    fw.write(data_join)
    fw.close()
    fastqc_data.close()
    
    

def get_basic_statistics (file_path):
    
    """ Read Basic_statistics contents from a fastqc file. 
        Returns a list of data.
        Data is printed to the user command line.  
        Acceptes one Aguments :-
                                 file_path = File name of fastqc """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Basic Statistics
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop' 

    for line in fastqc_data:
        if line.startswith(">>Basic Statistics"):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per base sequence quality'):
            break

        elif read_line == 'start':
            data.append(line)
    fastqc_data.close()
    data = "\n".join(data)
    for line in data.splitlines():
        if not line.startswith('#'):
            print(line)
    return data

def Per_base_sequence_quality (file_path, output_dire, module):
    
    """ Read Per base sequence quality contents from a fastqc file and parses to output file. 
        Returns a list of data.
        Data is used for plotting graphs which are saved to output directory as png image file.  
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """
    cwd = os.getcwd()
    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Per base sequence quality
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        if line.startswith('>>Per base sequence quality'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per tile sequence quality'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.
        
        if not module == "ALL":
            try:
                os.makedirs(output_dire)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                        raise
        # creating file path for both image and text file output
        if module == "ALL":
            output_path_dire_txt = os.path.join(cwd, "PBSQ", 'Per_base_sequence_quality.txt')
            output_path_dire_img = os.path.join(cwd, "PBSQ", 'Per_base_sequence_quality.png')
        else:
            output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_base_sequence_quality.txt')
            output_path_dire_img = os.path.join(cwd, output_dire, 'Per_base_sequence_quality.png')
            
    data.insert(0,head)#insert header before writing 
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
    
   
    # Loading and manuplating data to dataframe to be used for ploting
    df = np.loadtxt(data, dtype='float', usecols = [0,1,2,3,4,5,6], skiprows=1)
    df1 = pd.DataFrame(df)#columns=Base', 'Mean', 'Median', 'Lower Quartile', 'Upper Quartile', '10th Percentile',  '90th Percentile'
    df_Me = df1.drop([0,2,3,4,5,6], axis=1)# Mean columns
    df_Me = df_Me.astype(np.float32)
    df_pd1 = df.astype(np.float32)
    df_pd2 = np.delete(df_pd1, [0,1], axis=1)
    df_pd1 = pd.DataFrame(df).astype(np.float32)
    df_pd1 = df_pd1.drop([0],axis=1)
    df_pd3 = np.transpose(df_pd1)
    
    #Plotting the data
    plt.figure(figsize=(25,10))
    plt.boxplot(df_pd3,whis= [0,20], sym='_')
    plt.plot(df_Me)    
    plt.yticks(np.arange(0, 45, step=4))
    plt.xlabel('Position in read(bp)')
    plt.ylabel('Quality score (Phred)')
    plt.ylim(0,45)
    plt.xlim([0,76])
    plt.axhspan(20, 28, facecolor='wheat', alpha=0.7)
    plt.axhspan(28, 45, facecolor='lightgreen', alpha=0.6)
    plt.axhspan(0, 20, facecolor='salmon', alpha=0.6)
    plt.title("Per Base Sequence Quality", fontsize=20)
    plt.savefig(output_path_dire_img, format='png', dpi=200)
    #plt.show()
    plt.close() 
    return data
    

def Per_tile_sequence_quality (file_path, output_dire, module):
    
    """ Read Per tile sequence quality contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot heatmap.
        Data and genereted graphs are saved to output directory as svg image file and plain text. 
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Per Tile Sequence Quality
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Per tile sequence quality'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per sequence quality scores'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "PTSQ", 'Per_tile_sequence_quality.txt')
        output_path_dire_img = os.path.join(cwd, "PTSQ", 'Per_tile_sequence_quality.svg')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_tile_sequence_quality.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Per_tile_sequence_quality.svg')
    
    data.insert(0,head)#insert header before writing 
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
   

    #load per tile data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1,2], skiprows=2)
    df1 = pd.DataFrame(df, columns= ['Tile', 'Base', 'Mean'])
    # group the df by Base position to get all the tile mean value for each base
    gb = df1.groupby("Base")
    df=[gb.get_group(x) for x in gb.groups]
    # iterate to get the mean for each base in tile 
    vb = []
    for i in df:
        d1 = i['Mean']
        vb.append(d1.reset_index(drop=True))
    # convert each mean for the number of tiles in each base position
    vf = pd.concat(vb, axis=1)
    col=np.arange(1,76)# Base position
    row=np.arange(1,113)# tiles
    vf1 = vf
    vf1.columns = col# rename columns by base position
    vf2 = vf1.rename(dict(enumerate(row)))
    vf2=vf2.sort_index(ascending=False)# orderig the rows for appopraite plotting
    
    #plot heatmap to see the tile quality across each base position
    sn = sns.heatmap(vf2, cmap="coolwarm_r", cbar=False, xticklabels= True, yticklabels = True)# plot heat map from seaborn
    sn.set_xticklabels(sn.get_xticklabels(), rotation=30)
    sn.set_title('Quality per tile', fontsize=15)
    sn.set_xlabel("Position in read(bp)", fontsize=20)
    sn.set_ylabel("Tiles", fontsize=20)
    sn.figure.set_figwidth(40)
    sn.figure.set_figheight(28)
    sn.figure.savefig(output_path_dire_img, format='svg', dpi=1200)
    
    return data


def Per_sequence_quality_scores (file_path,output_dire, module):
    
    """ Read Per sequence quality scores contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot quality graph.
        Data and genereted graphs are saved to output directory as png image file and plain text. 
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Per Sequence Quality Scores
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Per sequence quality scores'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per base sequence content'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "PSQS", 'Per_sequence_quality_scores.txt')
        output_path_dire_img = os.path.join(cwd, "PSQS", 'Per_sequence_quality_scores.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_sequence_quality_scores.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Per_sequence_quality_scores.png')
    
    data.insert(0,head)#insert header before writing 
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
    
    
    #load per sequence quality scores data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1], skiprows=2)
    df1 = pd.DataFrame(df, columns= ['Quality', 'Count'])
    quality  =df1.iloc[:,0]
    count = df1.iloc[:,1]
    count = count.astype('int64')
    
    plt.figure(figsize=(20,10))
    plt.plot(quality,count, c = 'b')
    plt.xticks(np.arange(0, int(max(quality)+ 2)))
    plt.grid(which='major', linestyle='--', linewidth='0.5', color='green')
    plt.title('Quality score distribution over all sequences')
    plt.xlabel("Mean Sequence Quality (Phred Score)")
    ax=plt.gca()
    ax.ticklabel_format(axis='y', style='plain', scilimits=(0,0), useOffset=None, useLocale=None, useMathText=None)
    plt.savefig(output_path_dire_img, format='png', dpi=300)
    plt.draw()
    #plt.show()
    plt.close()
    
    return data


def Per_base_sequence_content (file_path, output_dire, module):
    
    """ Read Per base sequence content contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot base sequence content graph.
        Data and genereted graphs are saved to output directory as png image file and plain text.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """
        
    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Per Base Sequence Content
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Per base sequence content'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per sequence GC content'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "PBSC", 'Per_base_sequence_content.txt')
        output_path_dire_img = os.path.join(cwd, "PBSC", 'Per_base_sequence_content.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_base_sequence_content.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Per_base_sequence_content.png')
     
    data.insert(0,head)#insert header before writing    
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
        

    #load per base sequence content data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1,2,3,4], skiprows=2)
    df1 = pd.DataFrame(df, columns= ["Base", "G", "A",	"T", "C"])
    Bases  =df1.iloc[:,1:5]
    
    
    #plot
    plt.figure(figsize=(16,10))
    plt.plot(Bases)
    plt.legend(("%G", "%A",	"%T", "%C"))
    plt.xlim(1,76)
    plt.ylim(0, 101)
    plt.yticks(np.arange(0,101, 10))
    plt.xticks(np.arange(1,76, 2))
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Position in read (bp)')
    plt.title('Sequence content across all bases')
    plt.savefig(output_path_dire_img, format='png', dpi=300)
    #plt.show()
    plt.close()
    
    return data


def Per_sequence_gc_content (file_path, output_dire, module):
    
    """ Read Per sequence GC content contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot sequence GC content graph.
        Data and genereted graphs are saved to output directory as png image file and plain text.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Per sequence GC content
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Per sequence GC content'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Per base N content'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "PSGCC", 'Per_sequence_GC_content.txt')
        output_path_dire_img = os.path.join(cwd, "PSGCC", 'Per_base_sequence_GC_content.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_base_sequence_GC_content.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Per_base_sequence_GC_content.png')
    
    data.insert(0,head)#insert header before writing
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
        

    #load Per sequence GC content data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1], skiprows=2)
    df1 = pd.DataFrame(df, columns= ["GC content", "count"])
    
    #plot
    plt.figure(figsize=(16,10))
    plt.plot(df1)
    ax=plt.gca()
    ax.ticklabel_format(axis='y', style='plain', scilimits=(0,0), useOffset=None, useLocale=None, useMathText=None)
    plt.xlim(0,100)
    plt.xticks(np.arange(0,101, 2))
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Mean GC content (%)')
    plt.title('GC distribution over all sequences')
    plt.savefig(output_path_dire_img, format='png', dpi=600)
    #plt.show()
    plt.close()
     
    return data 


def Per_base_N_content (file_path, output_dire, module): 
    
    """ Read Per base N content contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot base N content graph.
        Data and genereted graphs are saved to output directory as png image file and plain text.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open('fastqc_data2.txt', 'r')

    # list to store Per base N content
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Per base N content'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Sequence Length Distribution'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "PBNC", 'Per_base_N_content.txt')
        output_path_dire_img = os.path.join(cwd, "PBNC", 'Per_base_N_content.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Per_base_N_content.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Per_base_N_content.png')
    
    data.insert(0,head)#insert header before writing   
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
       

    #load Per base N content data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1], skiprows=2)
    df1 = pd.DataFrame(df, columns= ["Base", "Count"])
    
    
    #plot
    plt.figure(figsize=(7,4))
    plt.plot(df1["Base"], df1['Count'], c='r', linewidth=4)
    plt.ylim(0,101)
    plt.ylim(1, 76)
    plt.yticks(np.arange(0,101, 10),fontsize=8)
    plt.xticks(np.arange(1,76, 3), fontsize=8)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Position in read (bp)', fontsize=8)
    plt.title('N content across all bases', fontsize=8)
    plt.legend(("N%"))
    plt.savefig(output_path_dire_img, format='png', dpi=600)
    plt.close()
    return data 

def Sequence_Duplication_Levels (file_path, output_dire, module):
    
    """ Read Sequence Duplication Levels contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot Sequence Duplication Levels graph.
        Data and genereted graphs are saved to output directory as png image file and plain text.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Sequence Duplication Levels
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Sequence Duplication Levels'):
            read_line = 'start' 
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Overrepresented sequences'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "SDL", 'Sequence_Duplication_Levels.txt')
        output_path_dire_img = os.path.join(cwd, "SDL", 'Sequence_Duplication_Levels.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Sequence_Duplication_Levels.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Sequence_Duplication_Levels.png')
    
    data.insert(0,head)#insert header before writing        
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()

    # Get Total Deduplicated Percentage that is going to be added to title of plot
    total_per = data[1]
    total_per = total_per.split('\t')
    total_per = round(float(total_per[1]), 2)
    #load Sequence Duplication Levels data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [1,2], skiprows=3)# first col contain string to avoid this will skip it.
    df_col1 = np.loadtxt(data, dtype='str', usecols = [0], skiprows=3)#get the first col that was skipped to label x axis plot.
    
    df1 = pd.DataFrame(df, columns= ['deduplicated_%', 'Total%'])# add column names

    #plot
    plt.figure(figsize=(10,5))
    plt.plot(df1, linewidth=1.5)
    plt.ylim(0,101)
    x = np.arange(0, 16)
    plt.yticks(np.arange(0,101, 10),fontsize=8)
    plt.xticks(x,df_col1, fontsize=8)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Sequence Duplication Level', fontsize=8)
    plt.title(('Percent of seqs remaining if deduplicated ' + str(total_per) + '%'), fontsize=8)
    plt.legend(("% Deduplicated sequences", "% Total sequences" ), fontsize=6)
    plt.savefig(output_path_dire_img, format='png', dpi=600)
    plt.close()
    
    return data 

def Adapter_Content (file_path, output_dire, module): 
    
    """ Read Adapter Content contents from a fastqc file and parses to output file. 
        Returns a list of data used to plot Adapter content graph.
        Data and genereted graphs are saved to output directory as png image file and plain text.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store Adapter Content
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Adapter Content'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Kmer Content'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "AC", 'Adapter_Content.txt')
        output_path_dire_img = os.path.join(cwd, "AC", 'Adapter_Content.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Adapter_Content.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Adapter_Content.png')
    
    data.insert(0,head)#insert header before writing 
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
       

    #load Adapter Content data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1,2,3,4], skiprows=2)# first col contain string to avoid this will skip it.
    df1 = pd.DataFrame(df, columns= ['Position', "IUAdapter", "ISRAdapter", "NTS", "SSRAdapter"])# add column names
    df2 = df1.drop(["Position"], axis=1)# droping position column before plotting
    
    #plot
    plt.figure(figsize=(10,5))
    plt.plot(df2, linewidth=4)
    plt.ylim(0,101)
    plt.xlim(1,64)
    plt.yticks(np.arange(0,101, 10),fontsize=8)
    plt.xticks(np.arange(1, 64, 3), fontsize=8)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Position in read (bp)', fontsize=8)
    plt.title(('% Adapter'), fontsize=8)
    plt.legend(('Illumina Universal Adapter', 'Illumina Small RNA Adapter', 'Nextera Transposase Sequence', 'SOLID Small RNA Adapter'), fontsize=6)
    plt.savefig(output_path_dire_img, format='png', dpi=600)
    plt.close()
    
    return data 

def Overrepresented_sequences (file_path, output_dire, module): 
    
    """ Read Overrepresented sequences contents from a fastqc file and parses to output file. 
        Returns a list of data parsed to a plain txt file.
        And saved to output directory.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open(file_path, 'r')

    # list to store  Overrepresented sequences
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Overrepresented sequences'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Adapter Content'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "OS", 'Overrepresented_sequences.txt')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Overrepresented_sequences.txt')
    
    data.insert(0,head)#insert header before writing 
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
  
def Sequence_Length_Distribution (file_path, output_dire, module):
    
    """ Read Sequence Length Distribution contents from a fastqc file and parses to output file. 
        Returns a list of data parsed to a plain txt file and used to plot graph.
        Output is saved to directory.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open('fastqc_data2.txt', 'r')

    # list to store Sequence Length Distribution
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'

    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Sequence Length Distribution'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # once the netxt module is reached it break the loop and stop reading lines
        elif line.startswith('>>Sequence Duplication Levels'):
            break

        elif read_line == 'start':
            data.append(line)
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "SLD", 'Sequence_Length_Distribution.txt')
        output_path_dire_img = os.path.join(cwd, "SLD", 'Sequence_Length_Distribution.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Sequence_Length_Distribution.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Sequence_Length_Distribution.png')
    
    data.insert(0,head)#insert header before writing    
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
       
    #load Sequence Length Distribution data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [0,1], skiprows=2)
    
    #plot
    sn=sns.distplot(df, hist=True,rug=True,kde=True,norm_hist=True)
    sn.set_title('Distribution of sequence lengths of overall sequences', fontsize=15)
    sn.set_xlabel("Position in read(bp)", fontsize=15)
    sn.figure.set_figwidth(13)
    sn.figure.set_figheight(8)
    sn.figure.savefig(output_path_dire_img, format='png', dpi=600)
    return data

def Kmer_Content(file_path, output_dire, module): 
    
    
    """Read Kmer content contents from a fastqc file and parses to output file. 
        Returns a list of data parsed to a plain txt file and used to plot graph.
        Output is saved to directory.
        Acceptes three Aguments :-
                                   file_path = File name of fastqc
                                   output_dire = Output Folder | directory name
                                   module = Fastqc module name """

    # opening fastqc data
    fastqc_data = open('fastqc_data2.txt', 'r')

    # list to store Kmer Content
    data = [] 

    # creating start/stop switch to read data for specific module
    read_line = 'stop'
    for line in fastqc_data:
        #line.strip("\n\r")
        if line.startswith('>>Kmer Content'):
            read_line = 'start'
            head = line
            pass

        elif line.startswith(">>END_MODULE"):
            read_line = 'stop'

        # # once the End module is reached it break the loop and stop reading lines, as kmer is the last module in the file
        elif line.startswith('>>END_MODULE'):
            break

        elif read_line == 'start':
            data.append(line)
                
    # Here it's creating directory and the try and except will raise any error except if the file exist.  
    try:
        os.makedirs(output_dire)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
    # creating file path for both image and text file output    
    if module == "ALL":
        output_path_dire_txt = os.path.join(cwd, "KC", 'Kmer_Content.txt')
        output_path_dire_img = os.path.join(cwd, "KC", 'Kmer_Content.png')
    else:
        output_path_dire_txt = os.path.join(cwd, output_dire, 'Kmer_Content.txt')
        output_path_dire_img = os.path.join(cwd, output_dire, 'Kmer_Content.png')
    
    data.insert(0,head)#insert header before writing           
    output_path_dire = cwd + '/' + output_dire + '/' + 'Kmer_Content.txt'
    fw = open(output_path_dire_txt, 'w' )
    fw.write("\n".join(data))
    fw.close()
    fastqc_data.close()
       

    #load Kmer_Content data to numpy dataframe 
    df = np.loadtxt(data, dtype='float', usecols = [1,2,3,4], skiprows=2)# first col is kmer string will skip it.
    df_kmer = np.loadtxt(data, dtype='str', usecols = [0], skiprows=2)# get col is kmer string.
    df1 = pd.DataFrame(df, columns= ["Count","PValue","Obs/Exp Max", "Max_OE_Position"])# add column names

    tf = np.transpose(df)# transpose
    #plot
    plt.figure(figsize=(10,5))
    plt.plot(tf, linewidth=2)
    plt.ylim(0,101)
    plt.xlim(1,22)
    plt.yticks(np.arange(0,101, 10),fontsize=8)
    plt.xticks(np.arange(1, 22, 3), fontsize=8)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.xlabel('Position in read (bp)', fontsize=8)
    plt.title(('Relative enrichment over read length'), fontsize=8)
    plt.legend((df_kmer), fontsize=6)
    plt.savefig(output_path_dire_img, format='png', dpi=1200)
    plt.close()
    return data 


def run(args):
    
    ''' Assigning corresponding Arguments passed to appropraite varibles.
        Calling functions and methods defined above to accept argument passed from command line 
        Returnig output file to user difined directory for each fastqc module.
        In case user passes ALL modules it will create all folders or Directory automaticall for each module. '''
    
    filename = args.input # these match the "dest": dest="input"
    oput_direc = args.output # from dest="output"
    module = args.module_name # from dest="module_name" and default is BS
   
    
    if module == 'BS':
        get_basic_statistics (filename)
    elif module == 'PBSQ':
        Per_base_sequence_quality(filename, oput_direc, module)    
    elif module == "PTSQ":
        Per_tile_sequence_quality(filename, oput_direc, module)
    elif module == "PSQS":
        Per_sequence_quality_scores(filename, oput_direc, module)
    elif module == "PBSC":
        Per_base_sequence_content (filename, oput_direc, module)
    elif module == "PSGCC":
        Per_sequence_gc_content(filename, oput_direc, module)
    elif module == "Filter":
        filter_report_generator(filename, oput_direc, module)
    elif module == "PBNC":
        Per_base_N_content(filename, oput_direc, module)
    elif module == "SDL":
        Sequence_Duplication_Levels(filename, oput_direc, module)
    elif module == "AC":
        Adapter_Content(filename, oput_direc, module)
    elif module == "OS":
        Overrepresented_sequences(filename, oput_direc, module)
    elif module == "SLD":
        Sequence_Length_Distribution(filename, oput_direc, module)
    elif module == "KC":
        Kmer_Content(filename, oput_direc, module)
    elif module == "ALL":
        folders = ["PBSQ", "PSQS","PBSC", "PSGCC", "PBNC", 'SLD', "SDL", "OS", "AC", "KC", "PTSQ", "Filter"]
        directories(folders)
        get_basic_statistics (filename)
        Per_base_sequence_quality(filename, oput_direc, module)
        Per_tile_sequence_quality(filename, oput_direc, module)
        Per_sequence_quality_scores(filename, oput_direc, module)
        Per_base_sequence_content (filename, oput_direc, module)
        Per_sequence_gc_content(filename, oput_direc, module)
        Per_base_N_content(filename, oput_direc, module)
        Sequence_Duplication_Levels(filename, oput_direc, module)
        Adapter_Content(filename, oput_direc, module)
        Overrepresented_sequences(filename, oput_direc, module)
        Sequence_Length_Distribution(filename, oput_direc, module)
        Kmer_Content(filename, oput_direc, module)
        filter_report_generator(filename, oput_direc, module)
        
        
    

        
def main():
	parser=argparse.ArgumentParser(description="FastQC Report parser, generates appropraite graphs and summary reoprt for the modules ")
	parser.add_argument("-in",help="fastqc input filename" ,dest="input", type=str, required=True)
	parser.add_argument("-out",help="fastqc output directory (folder_)" ,dest="output", type=str, required=True)
	parser.add_argument("-module",help=('''module names and filter = 
                                     BS: Base Statistics, 
                                     PBSQ: Per Base Sequence Quality, 
                                     PTSQ: Per Tile Sequence Quality, 
                                     PSQS: Per Sequence Quality Score, 
                                     PBSC: Per Base Sequence Content, 
                                     PSGCC: Per Sequence GC Content, 
                                     PBNC: Per Base N Content, 
                                     SLD: Sequence Length Distribution, 
                                     SDL: Sequence Duplication Levels, 
                                     OS: Overrepresented Sequence, 
                                     AC: Adapter Content, 
                                     KC: Kmer Content, 
                                     ALL: All modules
                                     Filter: filter report''') ,
                     dest="module_name", type=str, choices = ["BS", "PBSQ", "PSQS","PBSC", "PSGCC", "PBNC", 'SLD', "SDL", "OS", "AC", "KC", "ALL", "Filter"], 
                     default="BS")
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)


if __name__=="__main__":
	main()