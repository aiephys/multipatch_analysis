import pandas as pd
import os
import numpy as np
#from PIL import Image
import time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# from tkMessageBox import *
#open csv and get connections.
#load data from csv created from extract_first_pulse_fit_from_DB.py
#df=pd.read_csv('average_fit_dynamic_weight_2018_10_29.csv')
#df=pd.read_csv('average_fit_latency_jitter_2018_11_01.csv')
#df=pd.read_csv('average_fit_vclamp_2018_11_15.csv') #note that this version did not have the zero weight cross-talk
df=pd.read_csv('/home/corinnet/workspace/aiephys/multipatch_analysis/analyses/average_fit_vclamp_2018_12_06.csv') #note that this version did not have the zero weight cross-talk
df['uid']=df.apply(lambda row: "%.3f" % float(row.uid), axis=1)

connected_df=df[(df.distance<1e10) & (df.boolean_connection==True)]
print(len(connected_df))
connected_df["good_fit"]=""
save_file_name='/home/corinnet/workspace/aiephys/multipatch_analysis/analyses/ML_connected_vclamp_2018_12_07.csv'
completed=[]
if os.path.isfile(save_file_name):
    df_done=pd.read_csv(save_file_name)
    df_done['uid']=df_done['uid'].astype(str)
    for ii, (index, row) in enumerate(df_done.iterrows()): 
        completed.append([str(row.uid), str(row.pre_cell_id), str(row.post_cell_id)])
    
    

#show the plot and decide if it is a good fit
for ii, (index, row) in enumerate(connected_df.iterrows()): 
    print(ii, index, row['good_fit'])
    if [str(row.uid), str(row.pre_cell_id), str(row.post_cell_id)] in completed:
        print('already done passing:', row.uid, row.pre_cell_id, row.post_cell_id)
        continue
    # NOTE: I tried all sorts of ways to get the figure to show in the front and then close 
    # automatically with a time delay using Image, Tkinter, ImageTK, subprocess etc.  None of them worked.
    plt.figure(figsize=(14,14))
    plt.ion() #for what ever reason this allows plt.pause to work 
    try:
        plt.imshow(mpimg.imread(row['image_path']))
        plt.show()
        plt.pause(2)  #this allows the figure to come to front for 2 seconds
        plt.close()
    
        good_fit=raw_input("Is this a good fit ('y': yes, 'n':no, '?':somewhere in between)? \nIf you need to see the plot again give a different random input")
        # if the awnser was not sufficient image is open again and doesn't automatically close
        if good_fit not in ['y', 'n', '?']:
            plt.figure(figsize=(14,14))
            plt.ioff()
            plt.imshow(mpimg.imread(row['image_path']))
            plt.show()
            good_fit=raw_input("Is this a good fit?")
        connected_df.at[index, 'good_fit']=good_fit
    except:
        connected_df.at[index, 'good_fit']='no image'        
        print(index, connected_df.loc[index]['good_fit'], connected_df.loc[index]['image_path'])
        plt.close()
    if ii%20 == 0:
        connected_df.to_csv(save_file_name)

connected_df.to_csv(save_file_name)

