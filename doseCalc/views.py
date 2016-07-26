from django.shortcuts import get_object_or_404, render, render_to_response
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.template import loader, RequestContext
from django.core.urlresolvers import reverse
from django.views import generic

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings

from .models import Setup
from .forms import CT
import os
import zipfile
import io

# Create your views here.
#def detail(request):
#	return HttpResponse("Dose calculation app.")

class IndexView(generic.TemplateView):
    model = Setup
    template_name = 'doseCalc/index.html'
	
	
	
def submit(request):	
	Setup.beamWidth[0] = float(request.POST.get('beamsize_x'))
	Setup.beamWidth[1] = float(request.POST.get('beamsize_y'))
	Setup.energy = float(request.POST.get('energy'))*12
	Setup.nHistories = int(request.POST.get('nhistories'))
	
	Setup.voxSize[0] = float(request.POST.get('voxsize_x'))
	Setup.voxSize[1] = float(request.POST.get('voxsize_y'))
	Setup.voxSize[2] = float(request.POST.get('voxsize_z'))
	Setup.resolution[0] = float(request.POST.get('resolution_x'))
	Setup.resolution[1] = float(request.POST.get('resolution_y'))
	Setup.resolution[2] = float(request.POST.get('resolution_z'))
	
	Setup.scoringQuantity = request.POST.get('scoring_quantity')
	
	Setup.CTFile = request.FILES.get('CT')
	if Setup.CTFile is not None:
		Setup.CTFile.name = "ctvolume.dat"
		path = default_storage.save('doseCalc/goCMC/'+Setup.CTFile.name, ContentFile(Setup.CTFile.read()))
	Setup.material = request.POST.get('material')
	
	return HttpResponseRedirect(reverse('doseCalc:calculating'))
	
def calculating(request):
	os.chdir("C:/Users/zhen/mysite1/mysite/doseCalc/goCMC")
	
	config = open("carbon_config", "w+")
	#config = open("test.txt", "w+")
	config.write('''// directory for macro cross section:
input/carbon.crossSection
// directory for mass stopping power ratio:
input/carbon.mspr
// directory for restricted stopping power in water:
input/carbon.rspw
// directory for nuclear data
input/nuclear
// beam parameters: width, energy, number of carbons, source center\n''')
	config.write("%f\t%f\n" % (Setup.beamWidth[0], Setup.beamWidth[1]))
	config.write("%f\n" % Setup.energy)
	config.write("%d\n" % Setup.nHistories)
	config.write("%f\t%f\t%f\n" % (0,0,-50))
	config.write("// phantom parameters: vox size, iso center, width, hight, depth\n")
	config.write("%f\t%f\t%f\n" % (Setup.voxSize[0],Setup.voxSize[1],Setup.voxSize[2]))
	config.write("%f\t%f\t%f\n" % (0,0,0))
	config.write("%d\t%d\t%d\n" % (Setup.resolution[0],Setup.resolution[1],Setup.resolution[2]))
	config.write("// CTFile or standard material:\n")
	if Setup.CTFile is not None:
		config.write(Setup.CTFile.name)
	else:
		config.write(Setup.material)
	config.write("\n// scoring quantity:\n")
	config.write(Setup.scoringQuantity)
	config.write('''
// output directory
Output/''')
	config.close()
	
	os.system('CLTransport.exe')
	
    # Files (local path) to put in the .zip
    # FIXME: Change this (get paths from DB etc)
	filenames = ["Output/"+Setup.scoringQuantity+"_mean.bin", "Output/"+Setup.scoringQuantity+"_std.bin"]

    # Folder name in ZIP archive which contains the above files
    # E.g [thearchive.zip]/somefiles/file2.txt
    # FIXME: Set this to something better
	zip_subdir = "Output"
	zip_filename = "%s.zip" % zip_subdir

    # Open StringIO to grab in-memory ZIP contents
	s = io.BytesIO()

    # The zip compressor
	zf = zipfile.ZipFile(s, "w")

	for fpath in filenames:
		# Calculate path for file in zip
		fdir, fname = os.path.split(fpath)
		zip_path = os.path.join(zip_subdir, fname)

        # Add file, at correct path
		zf.write(fpath, zip_path)

    # Must close zip for all contents to be written
	zf.close()

	# Grab ZIP file from in-memory, make response with correct MIME-type
	resp = HttpResponse(s.getvalue(), content_type = "application/x-zip-compressed")
	# ..and correct content-disposition
	resp['Content-Disposition'] = 'attachment; filename=%s' % zip_filename
	if Setup.CTFile is not None:
		os.remove(Setup.CTFile.name)
	os.chdir("C:/Users/zhen/mysite1/mysite")
	return resp
	

	
