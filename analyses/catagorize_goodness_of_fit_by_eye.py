import pandas as pd
import os
import numpy as np
#from PIL import Image
import time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# from tkMessageBox import *

save_file_name='/home/corinnet/workspace/aiephys/multipatch_analysis/analyses/ML_connected_iclamp_2018_12_18.csv'
if not os.path.isfile(save_file_name): #create the file if it doesnt exist
    #load data from csv created from extract_first_pulse_fit_from_DB.py
    #df=pd.read_csv('average_fit_dynamic_weight_2018_10_29.csv')
    #df=pd.read_csv('average_fit_latency_jitter_2018_11_01.csv')
    #df=pd.read_csv('average_fit_vclamp_2018_11_15.csv') #note that this version did not have the zero weight cross-talk
    #df=pd.read_csv('/home/corinnet/workspace/aiephys/multipatch_analysis/analyses/average_fit_vclamp_2018_12_06.csv') #note that this version did not have the zero weight cross-talk
    df=pd.read_csv('/home/corinnet/workspace/aiephys/multipatch_analysis/analyses/dynamic_weight_jitter_latency2018-11-12.csv') #note that this version did not have the zero weight cross-talk
    df['uid']=df.apply(lambda row: "%.3f" % float(row.uid), axis=1)

    connected_df=df[(df.distance<1e10) & (df.boolean_connection==True)]
    print(len(connected_df))
    connected_df["good_fit"]=""
    connected_df["data_clarity"]=""
    connected_df.to_csv(save_file_name)

connected_df=pd.read_csv(save_file_name)
connected_df['uid']=connected_df['uid'].astype(str)
completed=[] 
connected_df["good_fit"]=connected_df["good_fit"].astype(str)
connected_df["data_clarity"]=connected_df["data_clarity"].astype(str)
for ii, (index, row) in enumerate(connected_df.iterrows()): 
    if row["good_fit"]!='nan': #checks for nans
        completed.append([str(row.uid), str(row.pre_cell_id), str(row.post_cell_id)]) #note that technically this isn't needed
print('completed',len(completed), 'out of', len(connected_df))
    

#show the plot and decide if it is a good fit
viewed=0
for ii, (index, row) in enumerate(connected_df.iterrows()): 
    print(ii, index, row['good_fit'])
    # if [str(row.uid), str(row.pre_cell_id), str(row.post_cell_id)] in completed:
    if row["good_fit"]!='nan': #note this might not work if this is the first file; just run again 
        print('already done passing:', row.uid, row.pre_cell_id, row.post_cell_id)
        continue
    # NOTE: I tried all sorts of ways to get the figure to show in the front and then close 
    # automatically with a time delay using Image, Tkinter, ImageTK, subprocess etc.  None of them worked.
    plt.figure(figsize=(14,14))
    plt.ion() #for what ever reason this allows plt.pause to work 
    try:
        plt.imshow(mpimg.imread(row['image_path']))
        plt.show()
        plt.pause(3)  #this allows the figure to come to front for 2 seconds
        plt.close()
        good_fit=raw_input("Is this a good fit ? \n\t'e': excellent \n\t'g': good \n\t'o': ok \n\t'b': bad \n\t'?': don't know what to think \n\t'd': double peak? \n\tIf you need to see the plot again give a different random input")
        # if the awnser was not sufficient image is open again and doesn't automatically close
        if good_fit not in ['e', 'g','o', 'b', '?', 'd']:
            plt.figure(figsize=(14,14))
            plt.ioff()
            plt.imshow(mpimg.imread(row['image_path']))
            plt.show()
            good_fit=raw_input("Is this a good fit?")
        if good_fit=='e':
            connected_df.at[index, 'good_fit']='excellent'
        elif good_fit=='g':
            connected_df.at[index, 'good_fit']='good'
        elif good_fit=='o':
            connected_df.at[index, 'good_fit']='ok'    
        elif good_fit=='b':
            connected_df.at[index, 'good_fit']='bad'
        elif good_fit=='?':
            connected_df.at[index, 'good_fit']='dont know what to think'    
        elif good_fit=='d':
            connected_df.at[index, 'good_fit']='double peak' 

        data_clarity=raw_input("How well can you see what is going on in this data ? \n\t'w': well \n\t's': sort of \n\t'n': not at all \n\t'?': don't know what to think")
        if data_clarity not in ['w', 's', 'n','?']:
            plt.figure(figsize=(14,14))
            plt.ioff()
            plt.imshow(mpimg.imread(row['image_path']))
            plt.show()
            data_clarity=raw_input("How easy is this data to fit?")
        if data_clarity=='w':
            connected_df.at[index, 'data_clarity']='well'
        elif data_clarity=='s':
            connected_df.at[index, 'data_clarity']='sort of'
        elif data_clarity=='n':
            connected_df.at[index, 'data_clarity']='not at all'    
        elif data_clarity=='?':
            connected_df.at[index, 'data_clarity']='dont know what to think'    

        

    except:
        connected_df.at[index, 'good_fit']='no image'   
        connected_df.at[index, 'data_clarity']='no image'        
        print(index, connected_df.loc[index]['good_fit'], connected_df.loc[index]['image_path'])
        plt.close()
    
    #save the file every 10 new assessments
    viewed += 1
    if viewed == 5:
        connected_df.to_csv(save_file_name)
        print('saved')
        viewed = 0

connected_df.to_csv(save_file_name)


