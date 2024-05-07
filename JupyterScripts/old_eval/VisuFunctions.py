import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import glob
import re

def plot_plate(Data, transfer, strategy, hue, hue_infos):
    sns.set_style("white")
    fig = plt.figure(str(transfer)+str(strategy), dpi=200)
    Data_ts = Data[(Data.transfer_n == transfer) & (Data.strategy == strategy) ].sort_values("row", ascending=False)
    fig = sns.scatterplot(data=Data_ts, 
                          x="col", y="row", 
                          s=100,
                          hue=hue, palette=hue_infos) #, order = plate_rows
    plt.legend(bbox_to_anchor=(0.75, -0.2), borderaxespad=0)
    fig.set_title(strategy+", Transfer "+str(transfer))
    return fig

def sort_image_pathes(image_path):
    ## Put images in right order
    image_pathes = pd.DataFrame(glob.glob(image_path+'*.png'), columns = ["image_path"])
    for i, image_p in image_pathes.iterrows():
        t = re.findall("\d+", image_p["image_path"])[-1]
        image_pathes.loc[i, "t"] = int(t)
    image_pathes = image_pathes.sort_values(by="t")
    image_pathes = image_pathes.reset_index()
    return image_pathes

def create_gif(folder_path, name):
    ## Get imagepahtes in the right order
    image_pathes = sort_image_pathes(folder_path)
    images = []
    # get all the images in the 'images for gif' folder
    for filename in image_pathes["image_path"]: # loop through all png files in the folder
        im = Image.open(filename) # open the image
        images.append(im) # add the image to the list

    # calculate the frame number of the last frame (ie the number of images)
    last_frame = (len(images)) 

    # create 10 extra copies of the last frame (to make the gif spend longer on the most recent data)
    for x in range(0, 9):
        im = images[last_frame-1]
        images.append(im)

    # save as a gif   
    images[0].save(folder_path+name+'.gif',
                   save_all=True, append_images=images[1:], optimize=True, duration=2000, loop=0)
