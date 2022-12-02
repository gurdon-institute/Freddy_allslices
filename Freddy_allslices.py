import math as maths
from os.path import expanduser

from ij import IJ, ImagePlus, ImageStack, WindowManager, Prefs
from ij.gui import Roi, TextRoi, ShapeRoi, Overlay
from ij.process import StackStatistics, Blitter, ImageProcessor, ByteProcessor, ShortProcessor, AutoThresholder, FloodFiller
from ij.plugin import Duplicator
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.measure import ResultsTable

from java.awt import Color
from java.io import File
from javax.swing import JFileChooser, JOptionPane


minA = 8.0 # µm²
maxA = 300.0

sigma = 2.5	# µm
k = 1.4



def DoG(ip, sigma, k):
	proc = ip.duplicate()
	sub = ip.duplicate()
	proc.blurGaussian(sigma)
	sub.blurGaussian(k*sigma)
	proc.copyBits(sub, 0,0, Blitter.SUBTRACT)
	return proc


def fillHoles(mask):
	width = mask.getWidth()
	height = mask.getHeight()
	ff = FloodFiller(mask)
	mask.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if mask.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if mask.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if mask.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if mask.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if mask.get(i)==127:
		    mask.set(i, 0)
		else:
		    mask.set(i, 255)


def watershed(ip, tol):
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, tol, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)


def getMask(ip, method):
	stats = ip.getStatistics()
	hist = ip.getHistogram(256)
	thresh = AutoThresholder().getThreshold( AutoThresholder.Method.Huang, hist )
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	ip.threshold(int(thresh))
	mask = ip.convertToByte(False)
	return mask


def getRois(mask):
	rois = []
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	return rois


def run(imp, all, allC, show):
	cal = imp.getCalibration()
	W = imp.getWidth()
	H = imp.getHeight()
	C = imp.getNChannels()
	Z = imp.getNSlices()

	sigmaPx = sigma/cal.pixelWidth

	ol = Overlay()
	z0 = 1 if all else imp.getSlice()
	z1 = Z+1 if all else z0+1
	for z in range(z0,z1):
		masks = [None for c in range(C+1)]
		c0 = 1 if allC else 4
		for c in range(c0,C+1):
			ip = imp.getStack().getProcessor(imp.getStackIndex(c,z,1))
			proc = DoG(ip, sigmaPx, k)
			mask = getMask(proc, AutoThresholder.Method.Otsu)
			mask.dilate()
			mask.erode()
			
			fillHoles(mask)
			watershed(mask, 0.5)
			masks[c] = mask.duplicate()
			channelrois = getRois(mask)

		if allC:	# Rois from union of all channel masks
			maskc2c3c4 = ByteProcessor(W, H)
			maskc2c3c4.copyBits(masks[2], 0,0, Blitter.ADD)
			maskc2c3c4.copyBits(masks[3], 0,0, Blitter.OR)
			maskc2c3c4.copyBits(masks[4], 0,0, Blitter.OR)
			watershed(maskc2c3c4, 0.5)
			
			for i in range(5):
				maskc2c3c4.dilate()
			for i in range(5):
				maskc2c3c4.erode()
			
			nucrois = getRois(maskc2c3c4)

		else:	# Rois from C4 only
			nucrois = getRois(masks[4])
		
		indexmask = ShortProcessor(W,H)
	
		
		for i,roi in enumerate(nucrois):
			roiA = roi.getStatistics().area * cal.pixelWidth * cal.pixelHeight
			if roiA >= minA and roiA <= maxA:
				if show:
					indexmask.setValue(i/float(len(nucrois))*pow(2,16))	 # make index masks
					indexmask.fill(roi)
				row = rt.getCounter()
				roi.setPosition(0,z,0)
				ol.add(roi)
				area = roi.getStatistics().area
				rect = roi.getBounds()
				perim = roi.getLength()
				feret = roi.getFeretValues()
				circ = 4*maths.pi*(area/(perim*perim))
				rt.setValue("Image", row, imp.getTitle())
				rt.setValue("X", row, rect.x+rect.width/2.0)
				rt.setValue("Y", row, rect.y+rect.height/2.0)
				rt.setValue("Z", row, z)
				rt.setValue("Area", row, area*cal.pixelWidth*cal.pixelHeight)
				rt.setValue("Max Feret", row, feret[0]*cal.pixelWidth)
				rt.setValue("Min Feret", row, feret[2]*cal.pixelWidth)
				rt.setValue("Feret Ratio", row, feret[2]/feret[0])
				rt.setValue("Circularity", row, circ)
				for c in range(1,C+1):
					ip = imp.getStack().getProcessor(imp.getStackIndex(c,z,1))
					ip.setRoi(roi)
					roiStats = ip.getStatistics()
					
					rt.setValue("C"+str(c)+"Mean", row, roiStats.mean)
					rt.setValue("C"+str(c)+"Min", row, roiStats.min)
					rt.setValue("C"+str(c)+"Max", row, roiStats.max)
					rt.setValue("C"+str(c)+"StdDev", row, roiStats.stdDev)
			else:
				roi.setStrokeColor(Color.RED)
				ol.add(roi)
		if show:
			imp.setOverlay(ol)
			indexmaskimp = ImagePlus("C2 C3 C4 union index mask Z"+str(z), indexmask)
			indexmaskimp.setCalibration(cal)
			indexmaskimp.show()

allC = False

rt = ResultsTable()
if WindowManager.getImageCount()>0:
	image = IJ.getImage()
	op = JOptionPane.showConfirmDialog(None, "Run on all slices?")
	allSlices = op == JOptionPane.YES_OPTION
	run(image, allSlices, allC, True)
	rt.show(image.getTitle()+" Results")
else:
	path = Prefs.get("Freddy_allslices.path", expanduser("~"))
	fc = JFileChooser()
	fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
	fc.setDialogTitle("Directory...")
	if(fc.showDialog(None, "OK")==JFileChooser.APPROVE_OPTION):
		path = fc.getSelectedFile().getAbsolutePath()
		Prefs.set("Freddy_allslices.path", path)
		d = File(path)
		for f in d.listFiles():
			if f.getName().endswith(".tif"):
				image = IJ.openImage(f.getAbsolutePath())
				run(image, True, allC, False)
	rt.show(path+" Results")
