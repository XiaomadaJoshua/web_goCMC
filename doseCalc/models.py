from django.db import models

# Create your models here.

class Setup(models.Model):
	"""Class for setting up simulation parameters."""
	beamWidth = [0, 0]
	energy = 0
	nHistories = 0
	
	voxSize = [0, 0, 0]
	resolution = [0, 0, 0]
	scoringQuantity = ""
	CTFile = ""
	material=""