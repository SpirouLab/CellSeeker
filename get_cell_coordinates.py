#get_cell_coordinates.py
#
#Combines all cell coordinates into a single CSV file for easy use in
#spreadsheet software such as Excel
#Additionally, calculates the center of each cell volume 
#
#Input:     numcells - Number of cells to gather coordinates for, numbered
#                  sequentially from one
#
#           input_path - Path to the ExtractedVolumes folder created by CellSeeker
#
#           output_path - Path to place the output file
#
#
#Output:    allcells.txt - CSV text file, one cell per line
#

def get_cell_coordinates(numcells,input_path,output_path):
    #number of exported cells (counting from 1)
    numcells=52
    
    #create a CSV file to write values to
    allcells=open('/media/feynman/3TBBackup/Scripts/CellSeeker/P2/ExtractedVolumes/allcells.txt','w')
    
    #initialize cell counter
    i=1
    
    #create string with titles
    cellstring='Cell Number,XMin,XMax,YMin,YMax,ZMin,ZMax,XCent,YCent,ZCent\n'
    
    #for each text file in all exported cells
    while i<=numcells:
    
        #add cell number to CSV string
        cellstring=cellstring + str(i) +','
    
        #format string for text file import
        cellopen='/media/feynman/3TBBackup/Scripts/CellSeeker/P2/ExtractedVolumes/Cell_' + str(i) + '/CellBody/coordinates.txt'
    
        #use with to put the text file into a list without having to close it
        with open(cellopen) as file:
            data=file.readlines()
            
            #get coordinate values and calculate centers
            xmin=int(data[0].split()[1])
            xmax=int(data[1].split()[1])
            ymin=int(data[2].split()[1])
            ymax=int(data[3].split()[1])
            zmin=int(data[4].split()[1])
            zmax=int(data[5].split()[1])
            xcent=round(float(xmax-xmin)/2)+xmin
            ycent=round(float(ymax-ymin)/2)+ymin
            zcent=round(float(zmax-zmin)/2)+zmin
            
            #append the values to the CSV string and add a line break (\n)
            cellstring=cellstring + str(xmin) + ',' + str(xmax) + ',' + str(ymin) + ',' + str(ymax) + ',' + str(zmin) + ',' + str(zmax) + ',' + str(xcent) + ',' + str(ycent)+ ',' + str(zcent) + '\n' 
            
        #increment the cell counter
        i=i+1
        
    #once all of the coordinates have been placed in the string,
    #write the string to the CSV file
    allcells.write(cellstring)
    allcells.close()