#open images from DF
import os
import numpy as np

def find_image_file(path, uid_desired, pre_cell_id, post_cell_id):
    """given a uid and the pre and post cell ids find the file.
    Note that this required the image file names to be saved with a specific format.
    For example 1513976168.03_7tlx3_3tlx3_average_fit.png

    Inputs
    ------
    path: string
        path to folder
    uid_desired: string
        uid looking for 
    pre_cell_id: string
        pre synaptic cell id looking for
    post_cell_id: string
        post synaptic cell id looking for
    
    Returns
    -------
    either a string specifying the path of the image or a np.nan if an image could not be found
    """
    path, uid_desired, pre_cell_id, post_cell_id
    uid_desired=str(np.round(np.float(uid_desired), 3))
    files=os.listdir(path)
    for f in files:
        uid_f, pre, post, _, _ = f.split('_') #split the file names into parts
        #convert file name uid into comparable format
        if len(uid_f) != 14:  #if the uid is not a 10 digit number with 2 decimal places 
            if len(uid_f) == 13:  #usually the zero gets left off during the file naming
                uid_f=uid_f+'0'
            elif len(uid_f) == 12:  #usually the zero gets left off during the file naming
                uid_f=uid_f+'00'
            else:
                raise Exception("Don't know how to deal with %s uid in file name" % uid_f)  
        if uid_f == uid_desired:  # when you find the matching uids check for matching cell ids
            if pre[0] == pre_cell_id:
                if post[0] == post_cell_id:
                    return os.path.join(path, f)
    return np.nan
        

        
#print(find_image_file("/home/corinnet/workspace/DBfit_pics_any_sign_10_1ish_2018", "1536781898.381", "2", "8"))  