'''
Author: Paul Holcomb
Date: 09/15/2015
Project: cellseeker.py (previously Finder.py)

Purpose: Given EM data and masks of nuclei in TIFF or HDF5 format,
         find all the centers of the Nuclei.
         EM data comes from Serial Block Face Electron Microscopy (SBEM)
         and is converted into TIFF files.  These files are then processes
         through Ilastik's Pixel Classification workflow to identify all
         probable nuclei. This program finds the centers of all probable nuclei,
         displays them for verification by a proofreader, and then crops an EM
         volume for each cell related to an identified and verified nucleus.

Known Issues: ????

Usage:  python cellseeker.py

'''
import numpy as np
import glob
import sys
import os
import tifffile
from PIL import Image
from PIL import ImageTk
import h5py
import Tkinter as tk
import ttk
from tkFileDialog import askopenfilename,askdirectory
from scipy import ndimage
import subprocess
import time
from platform import system

#Main frame call

class CellSeekerGUI(tk.Tk):

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)
        self.iconbitmap('@cellseeker_icon.xbm')

        container.grid()

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.wm_title("CellSeeker: Initialization")

        self.frames = {}

        self.nucdir=tk.StringVar()
        self.emdir=tk.StringVar()
        self.xypix=tk.StringVar()
        self.zpix=tk.StringVar()
        self.celldiam=tk.StringVar()
        self.cellbuffer=tk.StringVar()
        self.progress = tk.StringVar()
        self.cellInfo = []
        self.nucarray = []
        self.emarray = []
        self.imgs = []
        self.verifiedcells=[]
        self.xextent=tk.StringVar()
        self.yextent=tk.StringVar()
        self.zextent=tk.StringVar()
        self.Upid=[]
        self.Dnid=[]
        self.Leftid=[]
        self.exdir=tk.StringVar()
        self.emcheck=tk.IntVar()
        self.rescale=[]
        self.em_xmax=[]
        self.em_ymax=[]

        for F in (StartUp, NucProcessingGUI, NucVerifyGUI, NucDisplayGUI, ExtractGUI, ExtractGUI2, EndGUI):
            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartUp)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()

        #run cell extraction if NucProcessingGUI window is visible

        if cont.__name__=="NucProcessingGUI":
            #start timer
            self.start_time=time.time()
            self.progress.set("Importing nucleus segmentation data...")
            self.wm_title("CellSeeker: Finding Cells...")
            self.update_idletasks()
            self.frames[NucProcessingGUI].update_idletasks()

            #Compute some variables necessary for later processing
            self.xextent=int(round((float(self.celldiam.get())*1000)/float(self.xypix.get())))
       	    self.yextent=int(self.xextent)
            self.zextent=int(round((float(self.celldiam.get())*1000)/float(self.zpix.get())))

            #Get list of identified cells
            self.nucarray,self.cellInfo =cellFind(self)

            #Save nucarray to a file for access later
            #tiffSave(self)

            #Get EM images and load into numpy array
            self.wm_title("CellSeeker: Getting EM...")
            self.progress.set("Finding cells finished. " + str(len(self.cellInfo)) + " possible cells found. Fetching EM images...")
            self.update_idletasks()
            self.frames[NucProcessingGUI].update_idletasks()
            self.emarray = loadImage(self.emdir.get(),0)

            #Proceed to cell verification
            print("--- %s seconds ---" % (time.time() - self.start_time))
            self.wm_title("CellSeeker: Cell Verification")
            self.progress.set("EM Images loaded. Beginning cell verification.")
            b = ttk.Button(self.frames[NucProcessingGUI], text = "Click to Continue.", command = lambda: self.show_frame(NucVerifyGUI))
            b.grid(column=0,row=1, pady=25)
            self.frames[NucProcessingGUI].grid()
            self.update_idletasks()
            self.frames[NucProcessingGUI].update_idletasks()

        if cont.__name__=="NucVerifyGUI":

            self.wm_title("CellSeeker: Cell Verification...")
            self.update_idletasks()
            self.frames[NucVerifyGUI].update_idletasks()

            self.imgs = getNucImages(self)
            self.verifiedcells = np.zeros(len(self.imgs))
            self.verifiedcells.fill(4)

            self.frames[NucDisplayGUI].tkraise()

            getVerified(self,0,'f')

        if cont.__name__=="ExtractGUI":

            self.wm_title("CellSeeker: Set Extraction Directory")
            self.update_idletasks()

        if cont.__name__=="ExtractGUI2":

            if self.emcheck.get()==1:
				xold=self.emarray.shape[0]

				with tifffile.TiffFile(self.emdir.get()) as imfile:
					[z,self.em_ymax,self.em_xmax]=imfile.series[0].shape

				self.rescale=round(self.em_xmax/float(xold))

            else:
				self.em_xmax=self.emarray.shape[0]
				self.em_ymax=self.emarray.shape[1]
				self.rescale=1

            self.wm_title("CellSeeker: EM Subvolume Extraction...")
            self.progress.set("Beginning EM subvolume extraction.")
            self.update_idletasks()

            checkdir(self.exdir.get())
            cropEM(self)

            self.frames[EndGUI].tkraise()

    def endprog(self):
        self.container.destroy()

#Sets up initial variables for CellSeeker
class StartUp(tk.Frame):
    def __init__(self, parent, controller):

        tk.Frame.__init__(self,parent)
        self.config(height=600,width=400)
        self.grid()
        self.controller=controller
        #controller.title("CellSeeker: Initialization")
        self.columnconfigure(1, minsize=300)

        #get nucleus segmentation file

        self.controller.nucdir.set('')

        nucmessage = ttk.Label(self, text = "Nucleus segmentation:")
        nucmessage.grid(column=0,row=0, sticky="E")
        nucdir = ttk.Entry(self, textvariable=self.controller.nucdir, width=40)
        nucdir.grid(column=1, columnspan=2,row=0, sticky="W")
        b = ttk.Button(self, text = "Search...", command = self.load_nuc_file)
        b.grid(column=3,row=0,sticky="W")

        #get EM image file

        self.controller.emdir.set('')
        
        emmessage = ttk.Label(self, text = "EM Images:")
        emmessage.grid(column=0,row=1, sticky="E")
        emdir = ttk.Entry(self, textvariable=self.controller.emdir, width=40)
        emdir.grid(column=1, columnspan=2,row=1, sticky="W")
        b = ttk.Button(self, text = "Search...", command = self.load_em_file)
        b.grid(column=3,row=1,sticky="W")

        #get pixel dimension

        self.controller.xypix.set("")
        
        xymessage = ttk.Label(self, text = "Pixel dimension (nm):")
        xymessage.grid(column=0,row=2, sticky="E")
        xypix = ttk.Entry(self, textvariable=self.controller.xypix, width=10)
        xypix.grid(column=1, columnspan=2, row=2, sticky="W")

        #get z-slice thickness

        self.controller.zpix.set("")
        
        zpixmessage = ttk.Label(self, text = "Z-slice thickness (nm):")
        zpixmessage.grid(column=0,row=3, sticky="E")
        zpix = ttk.Entry(self, textvariable=self.controller.zpix, width=10)
        zpix.grid(column=1, columnspan=2, row=3, sticky="W")

        #get cell diameter

        self.controller.celldiam.set("")

        cdmessage = ttk.Label(self, text = "Estimated Cell Diameter (um):")
        cdmessage.grid(column=0,row=4, sticky="E")
        cd = ttk.Entry(self, textvariable=self.controller.celldiam, width=10)
        cd.grid(column=1, columnspan=2, row=4, sticky="W")

        #get buffer (% of cell diameter to add to window size)

        self.controller.cellbuffer.set("")

        cbmessage = ttk.Label(self, text = "Buffer for Edges (% of cell diameter; 0-1):")
        cbmessage.grid(column=0,row=5, sticky="E")
        cb = ttk.Entry(self, textvariable=self.controller.cellbuffer, width=10)
        cb.grid(column=1, columnspan=2, row=5, sticky="W")

        #set all initialization variables and continue
        b = ttk.Button(self, text = "OK", command = lambda: self.controller.show_frame(NucProcessingGUI))
        b.grid(column=0, columnspan=4, row=6, pady=10)

    def load_nuc_file(self):
        fname = askopenfilename(initialdir="/home", filetypes=(("TIFF file(s)","*.tiff"),("TIFF file(s)","*.tif"),("HDF5 file","*.h5")))
        self.controller.nucdir.set(fname)

    def load_em_file(self):
        fname = askopenfilename(initialdir="/home", filetypes=(("TIFF file(s)","*.tiff"),("TIFF file(s)","*.tif"),("HDF5 file","*.h5")))
        self.controller.emdir.set(fname)

#show nucleus segmentation processing progress

class NucProcessingGUI(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid()
        self.controller=controller

        nucmessage = tk.Label(self, textvariable = self.controller.progress)
        nucmessage.grid(column=0,row=0, sticky="NSEW")

class NucVerifyGUI(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1,uniform='foo')
        self.grid_columnconfigure(1, weight=1,uniform='foo')
        self.grid_columnconfigure(2, weight=1,uniform='foo')
        self.grid()
        self.controller=controller

        nucmessage = tk.Label(self, textvariable = self.controller.progress)
        nucmessage.grid(column=0,row=0,columnspan=3,rowspan=2,sticky="NSEW")

class NucDisplayGUI(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid()
        self.controller=controller
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

        self.image1=[]
        self.image2=[]
        self.image3=[]

        self.nucmessage = tk.Label(self, text='{X,Y}')
        self.nucmessage2 = tk.Label(self,text='{X,Z}')
        self.nucmessage3 = tk.Label(self,text='{Y,Z}')
        self.backbutton=tk.Button(self, text="Back", state='disabled')
        self.approvebutton=tk.Button(self, text="Approve", state='disabled')
        self.rejectbutton=tk.Button(self, text="Reject", state='disabled')

        self.nucmessage.grid(column=0,row=0,ipadx=10,ipady=10)
        self.nucmessage2.grid(column=1,row=0,ipadx=10,ipady=10)
        self.nucmessage3.grid(column=2,row=0,ipadx=10,ipady=10)
        self.backbutton.grid(column=0,row=1)
        self.approvebutton.grid(column=1,row=1)
        self.rejectbutton.grid(column=2,row=1)
        self.update_idletasks()

class ExtractGUI(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid()
        self.controller=controller

        #get directory for subvolume export

        path=self.controller.nucdir.get()
        path=path.replace('\\','/')
        path=path[0:(path.rfind('/')-1)]
        self.controller.exdir.set(path)

        nucmessage = ttk.Label(self, text = "Choose a path for export:")
        nucmessage.grid(column=0,row=0, sticky="E")
        exdir = ttk.Entry(self, textvariable=self.controller.exdir, width=40)
        exdir.grid(column=1, row=0, sticky="W")
        searchdir = ttk.Button(self, text = "Search...", command = self.set_extract_dir)
        searchdir.grid(column=2,row=0,sticky="W")

        #add option for subvolumes from different EM volume
        emmessage = ttk.Label(self, text = "Use different EM volume for subvolumes:")
        emmessage.grid(column=0,row=1,sticky="E")
        emcheck = tk.Checkbutton(self,variable=self.controller.emcheck)
        emcheck.grid(column=1, columnspan=2,row=1, sticky="W")
        emmessage2 = ttk.Label(self, text = "New EM File:")
        emmessage2.grid(column=0,row=2,sticky="E")
        emdir2 = ttk.Entry(self, textvariable=self.controller.emdir, width=40)
        emdir2.grid(column=1,row=2, sticky="W")
        b = ttk.Button(self, text = "Search...", command = self.load_em_file)
        b.grid(column=2,row=2,sticky="W")

        #set extraction path and continue
        okbutton = ttk.Button(self, text = "OK", command = lambda: self.controller.show_frame(ExtractGUI2))
        okbutton.grid(column=0,columnspan=3,row=3, pady=10)


    def set_extract_dir(self):
        fname = askdirectory(initialdir=self.controller.exdir, title='Select directory for subvolume export...')
        self.controller.exdir.set(fname)

    def load_em_file(self):
        fname = askopenfilename(initialdir=self.controller.exdir, filetypes=(("TIFF file(s)","*.tiff"),("TIFF file(s)","*.tif"),("HDF5 file","*.h5")))
        self.controller.emdir.set(fname)
#show extraction progress

class ExtractGUI2(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid()
        self.controller=controller

        nucmessage = tk.Label(self, textvariable = self.controller.progress)
        nucmessage.grid(column=0,row=0, sticky="NSEW")

class EndGUI(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1,uniform='foo')
        self.grid()
        self.controller=controller

        #end program

        nucmessage = ttk.Label(self, text = "Export complete!")
        nucmessage.grid(column=0,row=0)

        #set extraction path and continue
        endbutton = ttk.Button(self, text = "Finish", command = lambda: self.controller.endprog)
        endbutton.grid(column=0,row=1, pady=10)


#  Find nuclei in Ilastik segmentation

def cellFind(parent):

    nucarray = loadImage(parent.nucdir.get(),1)

    nucarray=filterNuclei2D(nucarray,parent)

    #create structured element to fill holes
    structel=ndimage.morphology.generate_binary_structure(3,3)

    #fill holes
    currentslice="Dilating..."
    parent.progress.set(currentslice)
    parent.frames[NucProcessingGUI].update_idletasks()

    nucarray=ndimage.morphology.binary_dilation(nucarray,structel,2)

    currentslice="Filling holes..."
    parent.progress.set(currentslice)
    parent.frames[NucProcessingGUI].update_idletasks()

    nucarray=ndimage.morphology.binary_fill_holes(nucarray,structel)

    currentslice="Eroding..."
    parent.progress.set(currentslice)
    parent.frames[NucProcessingGUI].update_idletasks()

    nucarray=ndimage.morphology.binary_erosion(nucarray,structel,2)

    #filter nuclei in 3D to combine connected components
    nucarray,sizes,centers=filterNuclei3D(nucarray,parent)
    buffer = float(parent.cellbuffer.get())/2

    sizes,centers,extents=filterPartialCells(nucarray,sizes,centers,buffer,parent)

    voxelvol=(float(parent.xypix.get()) ** 2) * float(parent.zpix.get())
    cellinfo=[]
    cellNum=0
    while cellNum!=len(sizes):
        cellinfo.append([cellNum+1,centers[cellNum],sizes[cellNum]*voxelvol,extents[cellNum]])
        cellNum += 1

    return nucarray,cellinfo

def filterNuclei2D(nucarray,parent):

    #2D filter in Z

    fullmask=np.zeros(nucarray.shape)
    slicecount=0

    while slicecount!=fullmask.shape[2]:
        currentslice="Processing Z-slice " + str(slicecount+1) + " of " + str(nucarray.shape[2])
        parent.progress.set(currentslice)
        parent.frames[NucProcessingGUI].update_idletasks()

        im=nucarray[:,:,slicecount]
        label_im, nb_labels = ndimage.label(im)
        sizes=ndimage.sum(im,label_im,range(nb_labels+1))
        #check for artifacts
        sortsize=np.copy(sizes)
        sortsize.sort()
        sortsize=sortsize[::-1]
        if len(sortsize)>1:
            if sortsize[1]>(0.4*sortsize[0]):
                maxsize=sortsize[1]
            else:
                maxsize=sortsize[0]
        else:
            maxsize=sortsize[0]
        mask_size = sizes < 0.4*maxsize
        mask_size2 = sizes > maxsize
        mask_size =np.add(mask_size,mask_size2)
        remove_pixel = mask_size[label_im]
        label_im[remove_pixel]=0

        #relabel the remaining regions
        labels=np.unique(label_im)
        label_im=np.searchsorted(labels,label_im)

        #save the result to our fullmask array
        fullmask[:,:,slicecount]=label_im
        slicecount += 1

    #2D filter in Y

    fullmask2=np.zeros(nucarray.shape)
    slicecount=0

    while slicecount!=fullmask.shape[1]:
        currentslice="Processing Y-slice " + str(slicecount+1) + " of " + str(nucarray.shape[1])
        parent.progress.set(currentslice)
        parent.frames[NucProcessingGUI].update_idletasks()

        im2=nucarray[:,slicecount,:]
        label_im2, nb_labels2 = ndimage.label(im2)
        sizes2=ndimage.sum(im2,label_im2,range(nb_labels2+1))
        #check for artifacts
        sortsize2=np.copy(sizes2)
        sortsize2.sort()
        sortsize2=sortsize2[::-1]
        if len(sortsize2)>1:
            if sortsize2[1]>(0.4*sortsize2[0]):
                maxsize2=sortsize2[1]
            else:
                maxsize2=sortsize2[0]
        else:
            maxsize2=sortsize2[0]
        mask_size3 = sizes2 < 0.4*maxsize2
        mask_size4 = sizes2 > maxsize2
        mask_size3 =np.add(mask_size3,mask_size4)
        remove_pixel2 = mask_size3[label_im2]
        label_im2[remove_pixel2]=0

        #relabel the remaining regions
        labels2=np.unique(label_im2)
        label_im2=np.searchsorted(labels2,label_im2)

        #save the result to our fullmask array
        fullmask2[:,slicecount,:]=label_im2
        slicecount += 1

    fullmask=np.add(fullmask,fullmask2)
    fullmask[fullmask>0]=1

    return fullmask

def filterNuclei3D(nucarray,parent):

    ThreeDmsg="Performing 3D filtering...please wait."
    parent.progress.set(ThreeDmsg)
    parent.frames[NucProcessingGUI].update_idletasks()

    #3D filter
    label_im2, nb_labels2 = ndimage.label(nucarray)
    labels2=np.array(range(1,nb_labels2+1))
    sizes2=ndimage.sum(nucarray,label_im2,range(1,nb_labels2+1))
    com = np.round(ndimage.measurements.center_of_mass(nucarray,label_im2,labels2))
    return label_im2,sizes2,com

def filterPartialCells(nucarray,sizes,centers,buffer,parent):

    ThreeDmsg="Filtering partial cells...please wait."
    parent.progress.set(ThreeDmsg)
    parent.frames[NucProcessingGUI].update_idletasks()

    xextent=parent.xextent
    yextent=parent.yextent
    zextent=parent.zextent

    keep=np.zeros(len(centers))
    extentlist=[]

    count=0
    for each in centers:
        #calculate min and max extents in x,y,z
        xmin=each[0]-(xextent/2)
        xmax=each[0]+(xextent/2)
        ymin=each[1]-(yextent/2)
        ymax=each[1]+(yextent/2)
        zmin=each[2]-(zextent/2)
        zmax=each[2]+(zextent/2)

        #check to see if any values exceed the data volume
        if xmin<0 or ymin<0 or zmin<0 or xmax>nucarray.shape[0] or ymax>nucarray.shape[1] or zmax>nucarray.shape[2]:
            keep[count]=0
        else:
            keep[count]=1

            #add buffer specified by user
            xmin=xmin-round(buffer*xextent/2)
            xmax=xmax+round(buffer*xextent/2)
            ymin=ymin-round(buffer*yextent/2)
            ymax=ymax-round(buffer*yextent/2)
            zmin=zmin-round(buffer*zextent/2)
            zmax=zmax-round(buffer*zextent/2)

            if xmin<0:
                xmin=0
            if ymin<0:
                ymin=0
            if zmin<0:
                zmin=0
            if xmax>nucarray.shape[0]:
                xmax=nucarray.shape[0]
            if ymax>nucarray.shape[1]:
                ymax=nucarray.shape[1]
            if zmax>nucarray.shape[2]:
                zmax=nucarray.shape[2]


            extentlist.append([xmin,xmax,ymin,ymax,zmin,zmax])

        count += 1

    sizes2=sizes[np.nonzero(keep)]
    centers2=centers[np.nonzero(keep)]

    return sizes2,centers2,extentlist

# checks to see if a directory exists, and if not, create it
def checkdir(path):
    p = os.path.dirname(path)
    if not os.path.exists(p):
        os.makedirs(p)

#display nuclei with EM image for verification

def loadImage(path,binarize):
	if os.path.isdir(path)==True:
		if path[len(path)-1]=="/":
		    path=path[:len(path)-2]
		a=glob.glob(path+'/*.tiff')
		b=glob.glob(path+'/*.tif')
		if a==[] and b==[]:
		    sys.exit('No usable files found. Please check path.')
		elif a!=[] and b==[]:
		    files=a
                elif a==[] and b!=[]:
                    files=b
		elif a!=[] and b!=[]:
		    sys.exit('Mix of .tiff and .tif files found. Please standardize and re-run.')

		count=0
		for each in files:
		    im=Image.open(each)                                       #open first image to load into numpy array
		    temparray=np.array(im)
		    temparray2=convertto8bit(temparray)#convert image to numpy array
                    temparray2=temparray2.astype(np.uint8)                      #convert array into unsigned 8-bit integer format
		    temparray2=np.expand_dims(temparray2,axis=2)                #add a third dimension to store additional images
		    if count==0:                                              #if this is the first time through the loop
		        imagearray2=temparray2                                  #initialize imagearray using our first image
		        count += 1                                            #increment the counter
		    else:                                                     #if this is any time through the loop but the first
		        imagearray2=np.append(imagearray2,temparray2,axis=2)     #append the current image to the stack of images

    	else:
		path=path.replace('\\','/')
		extension = path[path.rfind('.')+1:len(path)]
		if extension=="tiff" or extension=="tif":

		    with tifffile.TiffFile(path) as imfile:
		        imagearray=imfile.asarray()
		        imagearray=np.swapaxes(imagearray,0,2)
		        imagearray2=convertto8bit(imagearray)
		        imagearray2=imagearray2.astype(np.uint8)

                elif extension=="h5":

                    grab = h5py.File(path,'r')
                    img = grab.get(grab.keys()[0])
                    img = np.array(img)
                    grab.close()
		    img = np.swapaxes(img,0,2)
		    imagearray = img.squeeze()
		    imagearray2=convertto8bit(imagearray)
		    imagearray2=imagearray2.astype(np.uint8)

		    if binarize==1:
                        imagearray2[imagearray2>1]=0

                else:
                    sys.exit('No usable files found. Please check path.')

	return imagearray2

def tiffSave(parent,array,path,cellname,cellNumber):
    array=np.swapaxes(array,0,1)
    numints=len(str(array.shape[2]))
    fullname=path+cellname
    w=0
    while w!=array.shape[2]:
        currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Exporting image " + str(w+1) + " of " + str(array.shape[2])
       	parent.progress.set(currentslice)
        parent.frames[ExtractGUI2].update_idletasks()
        intdiff=numints-len(str(w))
        ints=0
        intstr=""
        while ints!=intdiff:
            intstr=intstr+"0"
            ints += 1
        intstr=intstr+str(w)
        tempfilename=fullname+"_slice"+intstr+".tiff"
        tempfile=array[:,:,w].astype("uint8")
        tifffile.imsave(tempfilename,tempfile)
        w += 1

    currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Stacking TIFFs"
    parent.progress.set(currentslice)
    parent.frames[ExtractGUI2].update_idletasks()

    mptiffstring="convert " + path + "*.tiff " + fullname +"_EM.tiff"
    subprocess.call(mptiffstring,shell=True)

    currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Cleaning up TIFFs"
    parent.progress.set(currentslice)
    parent.frames[ExtractGUI2].update_idletasks()

    tempfiles=glob.glob(fullname+"_slice*.tiff")
    for f in tempfiles:
        os.remove(f)

def histeq(im,nbr_bins): #from http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html

   #get image histogram
   imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = 255 * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   im2 = np.interp(im.flatten(),bins[:-1],cdf)

   return im2.reshape(im.shape), cdf

def convertto8bit(arr):
    imagemax=arr.max()
    eightbit=(2**8)-1
    sixteenbit=(2**16)-1
    thirtytwobit=(2**32)-1
    sixtyfourbit=(2**64)-1
    if imagemax>eightbit and imagemax<=sixteenbit:          #16 bit image
        arr=arr*(float(eightbit)/sixteenbit)                #convert to 8 bit
    elif imagemax>sixteenbit and imagemax<=thirtytwobit:    #32 bit image
        arr=arr*(float(eightbit)/thirtytwobit)              #convert to 8 bit
    elif imagemax>thirtytwobit and imagemax<=sixtyfourbit:  #64 bit image
        arr=arr*(float(eightbit)/sixtyfourbit)              #convert to 8 bit

    return arr

def Precarve(parent,path,cellNumber):
    currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Converting to H5"
    parent.progress.set(currentslice)
    parent.frames[ExtractGUI2].update_idletasks()

    #Get system OS to write correct output for Ilastik
    ostype=system()

    #Create filename for Ilastik project file
    projfile=path[:path.rfind('.')] + ".ilp"

    #Convert to HDF5

    #Create filename for hdf5 file
    hdf5file=path[:path.rfind('.')] + "_Input.h5"

    #Create temporary project for file conversion
    tempproj=path[:path.rfind('.')] + "_tempproj.ilp"

    if ostype=="Windows":
  	  ilastikcmd="ilastik.exe --headless --new_project=" + tempproj + " --workflow=DataConversionWorkflow --input_axes=zyx --output_format=hdf5 " + path
    elif ostype=="Linux":
  	  ilastikcmd="run_ilastik.sh --headless --new_project=" + tempproj + " --workflow=DataConversionWorkflow --input_axes=zyx --output_format=hdf5 " + path
    elif ostype=="darwin":
  	  ilastikcmd="./ilastik.app/Contents/MacOS/ilastik --headless --new_project=" + tempproj + " --workflow=DataConversionWorkflow --input_axes=zyx --output_format=hdf5 " + path

    print ilastikcmd

    with open('h5conversionout.txt',"w") as h5log:
        subprocess.call(ilastikcmd, shell=True, stdout=h5log)

    #Create Ilastik project file and preprocess data

    currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Preprocessing for carving..."
    parent.progress.set(currentslice)
    parent.frames[ExtractGUI2].update_idletasks()

    if ostype=="Windows":
        ilastikcmd2="ilastik.exe --headless --new_project=" + projfile + " --workflow=Carving --run-preprocessing --preprocessing-sigma=1.2 --preprocessing-filter=step-edges " + hdf5file
    elif ostype=="Linux":
        ilastikcmd2="run_ilastik.sh --headless --new_project=" + projfile + " --workflow=Carving --run-preprocessing --preprocessing-sigma=1.2 --preprocessing-filter=step-edges " + hdf5file
    elif ostype=="darwin":
        ilastikcmd2="./ilastik.app/Contents/MacOS/ilastik --headless --new_project=" + projfile + " --workflow=Carving --run-preprocessing --preprocessing-sigma=1.2 --preprocessing-filter=step-edges " + hdf5file

    with open('preprocessout.txt',"w") as pplog:
        subprocess.call(ilastikcmd2, shell=True, stdout=pplog)

def getNucImages(parent):
        #Define local variables
        nucarray=parent.nucarray
        emarray=parent.emarray
        cellInfo=parent.cellInfo
        imgs = []
        count=0
	while count!=len(cellInfo):
           currentslice="Extracting nucleus image " + str(count+1) + " of " + str(len(cellInfo))
	   parent.progress.set(currentslice)
	   parent.frames[NucVerifyGUI].update_idletasks()
	   temparray=[]

           x = cellInfo[count][1][0]
	   y = cellInfo[count][1][1]
	   z = cellInfo[count][1][2]

	   xmin=int(round(cellInfo[count][3][0]))
	   xmax=int(round(cellInfo[count][3][1]))
	   ymin=int(round(cellInfo[count][3][2]))
	   ymax=int(round(cellInfo[count][3][3]))
	   zmin=int(round(cellInfo[count][3][4]))
	   zmax=int(round(cellInfo[count][3][5]))

	   if count==0:
            print str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " " + str(zmin) + " " + str(zmax)

	   #for x,y
	   src = np.copy(nucarray[xmin:xmax,ymin:ymax,z])
	   src2 = np.copy(emarray[xmin:xmax,ymin:ymax,z])

	   mask_nuc1 = np.ma.masked_not_equal(src, src[x-xmin,y-ymin],copy=True).mask
	   mask_nuc2 = np.ma.masked_equal(src, src[x-xmin,y-ymin],copy=True).mask
           rgbimage1=np.zeros([src.shape[0],src.shape[1],4],dtype='uint8')
	   foralpha1=rgbimage1[:,:,3]
	   foralpha1[:,:]=140

           if np.sum(mask_nuc2)>(0.5*src.size):
               src[mask_nuc1]=254
               src[mask_nuc2]=0
               foralpha1[mask_nuc2]=0
           else:
               src[mask_nuc1]=0
               src[mask_nuc2]=254
               foralpha1[mask_nuc1]=0

           rgbimage1[:,:,0] = src
	   new_nuc_image1=Image.fromarray(rgbimage1,mode="RGBA")

           rgbimage2=np.empty([src2.shape[0],src2.shape[1],4],dtype='uint8')
           src2,cdf=histeq(src2,254)
           rgbimage2[:,:,0] = np.copy(src2)
	   rgbimage2[:,:,1] = rgbimage2[:,:,2] = rgbimage2[:,:,0]
	   rgbimage2[:,:,3] = 254
	   new_em_image1=Image.fromarray(rgbimage2,mode="RGBA")

	   new_image1=Image.alpha_composite(new_em_image1,new_nuc_image1)
	   temparray.append(new_image1)

	   #for x,z
	   src3 = np.copy(nucarray[xmin:xmax,y,zmin:zmax])
	   src4 = np.copy(emarray[xmin:xmax,y,zmin:zmax])

	   mask_nuc3 = np.ma.masked_not_equal(src3, src3[x-xmin,z-zmin],copy=True).mask
	   mask_nuc4 = np.ma.masked_equal(src3, src3[x-xmin,z-zmin],copy=True).mask
           rgbimage3=np.zeros([src3.shape[0],src3.shape[1],4],dtype='uint8')
	   foralpha2=rgbimage3[:,:,3]
	   foralpha2[:,:]=140

           if np.sum(mask_nuc4)>(0.5*src3.size):
               src3[mask_nuc3]=254
               src3[mask_nuc4]=0
               foralpha2[mask_nuc4]=0
           else:
               src3[mask_nuc3]=0
               src3[mask_nuc4]=254
               foralpha2[mask_nuc3]=0

           rgbimage3[:,:,0] = src3
	   new_nuc_image2=Image.fromarray(rgbimage3,mode="RGBA")

           rgbimage4=np.empty([src4.shape[0],src4.shape[1],4],dtype='uint8')
           src4,cdf=histeq(src4,254)
           rgbimage4[:,:,0] = np.copy(src4)
	   rgbimage4[:,:,1] = rgbimage4[:,:,2] = rgbimage4[:,:,0]
	   rgbimage4[:,:,3] = 254
	   new_em_image2=Image.fromarray(rgbimage4,mode="RGBA")

	   new_image2=Image.alpha_composite(new_em_image2,new_nuc_image2)
	   temparray.append(new_image2)

	   #for y,z
	   src5 = np.copy(nucarray[x,ymin:ymax,zmin:zmax])
	   src6 = np.copy(emarray[x,ymin:ymax,zmin:zmax])

	   mask_nuc5 = np.ma.masked_not_equal(src5, src5[y-ymin,z-zmin],copy=True).mask
	   mask_nuc6 = np.ma.masked_equal(src5, src5[y-ymin,z-zmin],copy=True).mask
           rgbimage5=np.zeros([src5.shape[0],src5.shape[1],4],dtype='uint8')
	   foralpha3=rgbimage5[:,:,3]
	   foralpha3[:,:]=140

           if np.sum(mask_nuc6)>(0.5*src5.size):
               src5[mask_nuc5]=254
               src5[mask_nuc6]=0
               foralpha3[mask_nuc6]=0
           else:
               src5[mask_nuc5]=0
               src5[mask_nuc6]=254
               foralpha3[mask_nuc5]=0

	   rgbimage5[:,:,0] = src5
	   new_nuc_image3=Image.fromarray(rgbimage5,mode="RGBA")

           rgbimage6=np.empty([src6.shape[0],src6.shape[1],4],dtype='uint8')
           src6,cdf=histeq(src6,254)
           rgbimage6[:,:,0] = np.copy(src6)
	   rgbimage6[:,:,1] = rgbimage6[:,:,2] = rgbimage6[:,:,0]
	   rgbimage6[:,:,3] = 254
	   new_em_image3=Image.fromarray(rgbimage6,mode="RGBA")

	   new_image3=Image.alpha_composite(new_em_image3,new_nuc_image3)
	   temparray.append(new_image3)

           imgs.append(temparray)
           count += 1
        return imgs

def getVerified(parent,current,lastvalue):      #For some reason, the hotkey code doesn't work yet
    #initialize variables
    if len(parent.verifiedcells)>1:
        if parent.verifiedcells[current]!=4:
            if parent.verifiedcells[current]==1:
                parent.wm_title("CellSeeker: Cell Verification...Accepted")
            else:
                parent.wm_title("CellSeeker: Cell Verification...Rejected")
            parent.update_idletasks()
            parent.frames[NucVerifyGUI].update_idletasks()


    if lastvalue!=2 and lastvalue!=3 and lastvalue!='f':
        parent.verifiedcells[current-1]=lastvalue

    if lastvalue==3:
        parent.frames[NucDisplayGUI].nucmessage.grid_forget()
        parent.frames[NucDisplayGUI].nucmessage2.grid_forget()
        parent.frames[NucDisplayGUI].nucmessage3.grid_forget()
        parent.frames[NucDisplayGUI].approvebutton.grid_forget()
        parent.frames[NucDisplayGUI].rejectbutton.grid_forget()
        parent.frames[NucDisplayGUI].update_idletasks()

    #if the last image hasn't been verified
    if current!=len(parent.cellInfo):
        #initialize key bindings and buttons based on image position
        if parent.Leftid!=[]:
            parent.unbind('<Up>',parent.Leftid)
        parent.Upid=parent.bind('<Up>',lambda: getVerified(parent,current+1,1))
        parent.Dnid=parent.bind('<Down>',lambda: getVerified(parent,current+1,0))
        parent.frames[NucDisplayGUI].approvebutton.config(command= lambda: getVerified(parent,current+1,1), state='normal')
        parent.frames[NucDisplayGUI].rejectbutton.config(command= lambda: getVerified(parent,current+1,0), text="Reject", state='normal')
        if current!=0:
            parent.Leftid=parent.bind('<Left>',lambda: getVerified(parent,current-1,2))
            parent.frames[NucDisplayGUI].backbutton.config(command= lambda: getVerified(parent,current-1,2), state='normal')
        else:
            parent.frames[NucDisplayGUI].backbutton.config(state='disabled')

        #get images
        XY=parent.imgs[current][0]
        XZ=parent.imgs[current][1]
        YZ=parent.imgs[current][2]

        #transform images into ImageTk format for display
        nucleus1=ImageTk.PhotoImage(XY)
        nucleus2=ImageTk.PhotoImage(XZ)
        nucleus3=ImageTk.PhotoImage(YZ)

        #for x,z
        parent.frames[NucDisplayGUI].nucmessage.config(image = nucleus1,text='{X,Y}',compound='bottom')#,width=grid1width)
        parent.frames[NucDisplayGUI].image1 = nucleus1

        parent.frames[NucDisplayGUI].nucmessage2.config(image = nucleus2,text='{X,Z}',compound='bottom')#,width=grid2width)
        parent.frames[NucDisplayGUI].image2 = nucleus2

        #for x,z
        parent.frames[NucDisplayGUI].nucmessage3.config(image = nucleus3,text='{Y,Z}',compound='bottom')#,width=grid3width)
        parent.frames[NucDisplayGUI].image3 = nucleus3

        parent.frames[NucDisplayGUI].nucmessage.grid(column=0,row=0,ipadx=10,ipady=10)
        parent.frames[NucDisplayGUI].nucmessage2.grid(column=1,row=0,ipadx=10,ipady=10)
        parent.frames[NucDisplayGUI].nucmessage3.grid(column=2,row=0,ipadx=10,ipady=10)
        parent.frames[NucDisplayGUI].backbutton.grid(column=0,row=1)
        parent.frames[NucDisplayGUI].approvebutton.grid(column=1,row=1)
        parent.frames[NucDisplayGUI].rejectbutton.grid(column=2,row=1)
        parent.frames[NucDisplayGUI].grid()
        parent.frames[NucDisplayGUI].update_idletasks()

    else:
        parent.verifiedcells[current-1]=lastvalue
        parent.frames[NucDisplayGUI].nucmessage.grid_forget()
        parent.frames[NucDisplayGUI].nucmessage2.grid_forget()
        parent.frames[NucDisplayGUI].nucmessage3.grid_forget()
        parent.frames[NucDisplayGUI].approvebutton.grid_forget()
        parent.frames[NucDisplayGUI].rejectbutton.grid_forget()
        parent.frames[NucDisplayGUI].image1 = []
        parent.frames[NucDisplayGUI].update_idletasks()

        parent.frames[NucDisplayGUI].nucmessage.config(image='', compound='none', text = "Nucleus verification complete. Ready to extract " + str(int(sum(parent.verifiedcells))) + " cells. Continue?")
        parent.frames[NucDisplayGUI].nucmessage.grid(column=0,row=0,columnspan=3,sticky="NSEW")

        parent.bind('<Left>',lambda: getVerified(parent,current-1,3))
        parent.unbind('<Up>',parent.Upid)
        parent.unbind('<Down>',parent.Dnid)
        parent.frames[NucDisplayGUI].backbutton.config(text="Back", command= lambda: getVerified(parent,current-1,3))
        parent.frames[NucDisplayGUI].backbutton.grid(column=0,row=1)



        parent.frames[NucDisplayGUI].rejectbutton.config(text="Begin Extraction", command= lambda: parent.show_frame(ExtractGUI))
        parent.frames[NucDisplayGUI].rejectbutton.grid(column=2,row=1)
        parent.frames[NucDisplayGUI].grid()
        parent.frames[NucDisplayGUI].update_idletasks()

def cropEM(parent):

        #Define local variables
        cellInfo=parent.cellInfo
        cellParts=['Nucleus','CellBody','Axon','Dendrites','Inputs','AdditionalStructures']
        cellNumber=1
        count=0

        #Extract EM subvolumes based on verified cells
        for each in parent.verifiedcells:
            if each==1:
                currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Building directory structure"
       	        parent.progress.set(currentslice)
       	        parent.frames[ExtractGUI2].update_idletasks()

                cellName='Cell_' + str(cellNumber)
                dirname=parent.exdir.get() + '/ExtractedVolumes/' + cellName +'/'
                for part in cellParts:
                    tempname=dirname + part + '/'
                    checkdir(tempname)

                currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Saving location information"
       	        parent.progress.set(currentslice)
       	        parent.frames[ExtractGUI2].update_idletasks()

                savepath=dirname+'CellBody/'
                xmin=int(round((cellInfo[count][3][0])*parent.rescale))
    	        xmax=int(round((cellInfo[count][3][1])*parent.rescale))
    	        ymin=parent.em_ymax-int(round((cellInfo[count][3][3])*parent.rescale)) #these have to be reversed due to a change in where newstack
    	        ymax=parent.em_ymax-int(round((cellInfo[count][3][2])*parent.rescale)) #thinks the origin for y is (top versus bottom of image)
    	        zmin=int(round(cellInfo[count][3][4]))
    	        zmax=int(round(cellInfo[count][3][5]))

    	        xcen=round(parent.em_xmax/2)
    	        ycen=round(parent.em_ymax/2)

    	        xsize=xmax-xmin
    	        ysize=ymax-ymin

    	        xoff=int(round((xmin+(xsize/2))-xcen))
    	        yoff=int(round((ymin+(ysize/2))-ycen))

    	        zstart=zmin-1
    	        zend=zmax-1


    	        cellcoords=savepath+'coordinates.txt'
    	        with open(cellcoords,'w') as textfile:
    	            textfile.write("xmin: " + str(xmin) + "\n")
    	            textfile.write("xmax: " + str(xmax) + "\n")
    	            textfile.write("ymin: " + str(ymin) + "\n")
    	            textfile.write("ymax: " + str(ymax) + "\n")
    	            textfile.write("zmin: " + str(zmin) + "\n")
    	            textfile.write("zmax: " + str(zmax) + "\n")

    	        currentslice="Exporting cell " + str(cellNumber) + " of " + str(int(np.sum(parent.verifiedcells))) + ": Saving location information"
    	        parent.progress.set(currentslice)
    	        parent.frames[ExtractGUI2].update_idletasks()


    	        #newstack is implemented using code from the IMOD software
    	        #produced by the Boulder Laboratory for 3-D Electron Microscopy
    	        #of Cells, University of Colorado: Boulder (bio3d.colorado.edu/imod)
    	        #under the GNU GPL (29 June 2007)

    	        myPath = os.path.realpath(".")
    	        outpath=savepath+cellName+"_EM.tiff"
    	        mptiffstring=myPath + "/newstack -siz " + str(xsize) + "," + str(ysize) +" -off " + str(xoff) + "," + str(yoff) + " -ori -sec " + str(zstart) + "-" + str(zend) + " " + parent.emdir.get() + " -fo TIF " + outpath
                print mptiffstring
    	        subprocess.call(mptiffstring,shell=True)
    	        Precarve(parent,outpath,cellNumber)
    	        cellNumber += 1

    	    count += 1



def main():
        #rewrite so that this is shown in Help menu (F1)
        """
        if len(sys.argv) < 1:
       	    print 'CellSeeker -- by Paul Holcomb'
       	    print 'Usage: python cellseeker.py [Nucleus Segmentation filename or directory] [EM filename or directory] [Pixel Dimension in nm] \
       	           [z-slice thickness in nm] [Scale factor for dataset (ie. bin1 -> bin 6 would be 6)] [Estimated cell diameter in um] \
       	           [Number of cell to extract (for testing)]'
       	    print 'Accepts TIFF or HDF5 image stacks or directory of individual TIFF images in a series (expects axes formatted as z,y,x as from Ilastik)'
       	    print 'Outputs processed subvolumes of EM for carving in Ilastik and masks for identified nuclei.'
        """
        cellseeker=CellSeekerGUI()
        cellseeker.eval('tk::PlaceWindow %s center' % cellseeker.winfo_pathname(cellseeker.winfo_id()))
        cellseeker.mainloop()

if __name__ == "__main__":
	main()
