# CellSeeker
<p>CellSeeker takes EM images and nuclear segmentation masks from Ilastik and produces EM subvolumes for individual cells as well as Ilastik carving project files. A menu-driven interface provides ease of use for inputting parameters and reviewing results.</p>

<h3>If you use CellSeeker, please cite:</h3>

<p>Holcomb PS, Morehead M, Doretto G, Chen P, Berg S, Plaza S, Spirou G. <i>“Rapid and Semi-Automatic Extraction of Neuronal Cell Bodies and Nuclei from Electron Microscopy Image Stacks.”</i> Methods Mol Biol., 2016 Jun; 1427: 277-290.</p>

<p>Sommer C, Strähle C, Köthe U, Hamprecht FA. <i>"ilastik: Interactive Learning and Segmentation Toolkit"</i> Eighth IEEE International Symposium on Biomedical Imaging (ISBI) Proceedings, (2011); 230-233.</p>

# Requirements

<h3>Software:</h3>
<p>Python 2.7<br>
ImageMagick<br>
ilastik 1.1.6 or later <br>
<b>NOTE:</b> CellSeeker expects that the Ilastik bin folder has been added to your path so that Ilastik can be run from the command prompt.</p>

<h3>Python dependencies:</h3>
<p>SciPy (https://www.scipy.org/)<br>
Pillow (https://pypi.python.org/pypi/Pillow/3.0.0)<br>
H5Py (http://www.h5py.org/)</p>

# Usage

<h3>To run:</h3>
<p>Open a terminal or command window, navigate to the folder containing the cellseeker.py file, and type:<br><br>
<i>python cellseeker.py</i><br><br></p>

<h3>Initialization:</h3>
<p>CellSeeker will first present an initialization window, where the user will input the following parameters:<br>
<br>
<b>Nucleus Segmentation:</b> the location of the nucleus segmentation in TIFF (series or multipage) or HDF5 format<br>
<b>EM Images:</b> the location of the electron microscopy images that correspond to the nucleus segmentation in TIFF (series or multipage) or HDF5 format<br>
<b>Pixel Dimension:</b> the size of the image pixels in either x or y in nanometers. Square pixels are assumed.<br>
<b>Z-slice thickness:</b> the distance between image planes in the depth (or z) dimension in nanometers.<br>
<b>Estimated Cell Diameter:</b> the expected diameter of the cells being searched for in microns (µm). This parameter may need to be adjusted in subsequent runs for optimal results.<br>
<b>Buffer for Edges:</b> the decimal percentage (0-1) of the cell diameter to enlarge the image window in order to provide a buffer and ensure that cells are completely captured. The larger this value, the more likely the complete cell body will be output, but at the cost of larger file sizes.</p><br>

<h3>Filtering and Cell Identification:</h3>
<p>After initialization, CellSeeker imports the provided EM and nucleus segmentation data.  A series of step-wise filtering and morphological operations are then performed to clean the segmentation data and identify individual nuclei as follows:<br>
<ul>
<li>2D size filtering in Z, removing components below a threshold of 40% of the largest connected component in each image</li>
<li>2D size filtering in y, removing components below a threshold of 40% of the largest connected component in each image</li>
<li>Binary ilation using a 3x3x3 structured element (two passes)</li>
<li>Hole filling</li>
<li>Binary erosion using a 3x3x3 structured element (two passes)</li>
<li>3D connected component labeling to identify nuclei and locate centers of mass</li>
<li>Partial cell filtering using cell diameter + buffer to keep only cells completely enclosed in volume <i>(currently hard-coded as on, may eventually be converted to on/off option in initialization)</i></li>
</ul></p><br>

<h3>Cell Verification:</h3>
<br>
<p>After cells with nuclei have been identified, each cell will be presented in three cross-sectional views (x,y,z) through the center of the image volume to be produced. Users have the option of approving or rejecting each cell individually using the buttons under each image set. This ensures that only good cells are output and artifacts are removed.</p><br>

<h3>Cell Output</h3>
<p>Following the verification of all found cells, CellSeeker runs a multi-step output process consisting of the following:<br>
<ul>
<li>Creation of the following directory structure in a user-defined location to store identified cell images and Ilastik projects:
<ul>
<li><b>ExtractedVolumes</b></li>
<ul>
<li>Cell01</li>
<ul>
<li><i>Nucleus</i></li>
<li><i>CellBody</i></li>
<li><i>Axon</i></li>
<li><i>Dendrites</i></li>
<li><i>Inputs</i></li>
<li><i>AdditionalStructures</i></li>
</ul>
<li>Cell02</li>
<ul>
<li>...</li>
</ul>
</ul>
</ul>
<li>Calculation of cell image stack origin within the coordinate system the whole EM volume (for placement of segmentation later). Origin is saved as CellBody/coordinates.txt.</li>
<li>Extraction of EM image multipage TIFF containing the current cell to CellBody/cell#_EM.tiff. (Currently accomplished using the <i>newstack</i> script written for the IMOD software produced by the Boulder Laboratory for 3-D Electron Microscopy of Cells, University of Colorado: Boulder (bio3d.colorado.edu/imod) under the GNU GPL (29 June 2007).)</li>
<li>Preprocessing of EM data for use in the "carving" workflow in Ilastik, which produces a project file saved to CellBody/cell#_EM.ilp. Preprocessing is accomplished using a headless Ilastik instance with the carving workflow and preprocessing settings of a step-edges filter for edge identification and a sigma size of 1.2 pixels (currently hard-coded, may be changed to user-controlled in future update).</li>
</ul></p>
<p><b>Note:</b> Folders other than the CellBody folder will be empty, and are produced simply for convenience for future segmentation.</p>

# Output Usage

<p>CellSeeker produces both a multipage TIFF stack of EM images which can be read by many segmentation softwares (Ilastik, IMOD, Seg3D, FIJI, etc.) and an Ilastik project file that can be used immediately in the Ilastik carving workflow.</p>






