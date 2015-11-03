from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
import numpy as np
import tifffile

def number_cells(celldata_path,em_path,output_path):
    #get size of EM image stack to build image stack
    with tifffile.TiffFile(em_path) as imfile:        #open the image stack
	[z,y,x]=imfile.series[0].shape                #get the dimensions
    #get locations for each cell number to write to image stack
    with open(celldata_path) as celldata:             #open the cell data file
        data=celldata.readlines()                     #read the file into individual lines
        cellcenters=[]                                #create a list to hold cell center info
        zslices=[]                                    #create a list to hold z locations
        for each in data:                             #for each cell in the list
            each=each.replace('\n','')                #remove line breaks
            each=each.split(',')                      #split the line into individual values
            cent=each[len(each)-3:len(each)]          #get the last three elements (XCenter,YCenter,Zcenter coordinates)
            if cent[0]!='XCent':                      #as long as this isn't the label row
                cent=[float(n) for n in cent]           #convert the values into floating point numbers
                zslice=float(each[len(each)-1])           #get the z coordinate
                cellcenters.append(cent)              #append the cell center coordinates to the list
                zslices.append(zslice)                #append the z coordinate location to the list
                    
    i=1                                               #create a counter  
    cellnumstack=np.empty((z,y,x),np.uint8)           #create an empty numpy array to hold image information
    fnt=ImageFont.truetype('/media/feynman/3TBBackup/Scripts/CellSeeker/fonts/ARIALN.TTF',150) #set font to use for cell numbering
    while i<=z:                                       #for each image in the stack,
        a=Image.new('L',(x,y),0)                      #create blank image
        currentloc=[m for m,n in enumerate(zslices) if n==i]    #get indices for cells located on this slice
        if currentloc!=[]:                            #if there is a cell located on this slice
            for cell in currentloc:
                d=ImageDraw.Draw(a)                              #get a drawing context
                [current_x,current_y]=cellcenters[cell][0:2]     #get x and y coordinates of cell center
                #current_x=current_x                              offset for x (if necessary)
                current_y=(y-current_y)                          #offset for y measured from the top of the image
                print i,current_x,current_y                     
                d.text((current_x,current_y),str(cell+1),font=fnt, fill=255) #write current cell number on image
        cellnumstack[i-1,:,:]=np.array(a,np.uint8)
        i=i+1
    tifffile.imsave(output_path,cellnumstack)
    return cellnumstack    

