# **************************************************************************
# *
# * Authors:         Jose Luis Vilas (jvilas@cnb.csic.es)
#                    Javier Vargas (jvargas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, isfile
from shutil import copyfile
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import Viewer
import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)


class XmippProtMultiRefAlignability(ProtAnalysis3D):
    """    
    Performs soft alignment validation of a set of particles confronting them
    against a given 3DEM map. This protocol produces particle alignment
    precision and accuracy parameters.
    """
    _label = 'multireference aligneability'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='SetOfVolumes,Volume',
                      label="Input volumes",  
                      help='Select the input volumes.')     
                
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignment',
                      label="Input particles", important=True,
                      help='Select the input projection images.')
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('angularSampling', FloatParam, default=5,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('numOrientations', FloatParam, default=10,
                      label="Number of Orientations for particle",  
                      help='Parameter to define the number of most similar volume \n' 
                      '    projected images for each projection image')
        
        form.addParam('phaseFlipped', BooleanParam, default=False,
                      label="Is the data already phase flipped?",
                      help='In case the data has been already phase flipped select True')

        form.addParam('doNotUseWeights', BooleanParam, default=False,
                      label="Do not use the weights",
                      help='Do not use the weights in the clustering calculation')
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):   
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.inputParticles.get().getObjId())
        deps = [] # store volumes steps id to use as dependencies for last step
        commonParams    = self._getCommonParams()
        commonParamsRef = self._getCommonParamsRef()
        sym = self.symmetryGroup.get()

        for i, vol in enumerate(self._iterInputVols()):
            volName = getImageLocation(vol)
            volDir = self._getVolDir(i+1)
            
            sigStepId = self._insertFunctionStep('significantStep', 
                                                 volName, volDir,
                                                 'exp_particles.xmd',
                                                 commonParams, 
                                                 prerequisites=[convertId])
            
            phanProjStepId = self._insertFunctionStep('phantomProject', 
                                                 volName, 
                                                 prerequisites=[convertId])

            sigStepId = self._insertFunctionStep('significantStep', 
                                                 volName, volDir,
                                                 'ref_particles.xmd',
                                                 commonParamsRef, 
                                                 prerequisites=[phanProjStepId])

            volStepId = self._insertFunctionStep('validationStep', 
                                                 volName, volDir,
                                                 sym,
                                                 prerequisites=[sigStepId])
            
            deps.append(volStepId)
            
        self._insertFunctionStep('createOutputStep', 
                                 prerequisites=deps)
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        a = self.inputParticles.get() 
        print a.printAll() 
        writeSetOfParticles(self.inputParticles.get(), 
                            self._getPath('input_particles.xmd'))
                    
    def _getCommonParams(self):
        params =  '  -i %s' % self._getPath('input_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --alpha0 0.05 --alphaF 0.05'
        params += ' --angularSampling %0.3f' % self.angularSampling.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %d' % self.numOrientations.get()
        return params
    
    def _getCommonParamsRef(self):
        params =  '  -i %s' % self._getPath('reference_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --alpha0 0.05 --alphaF 0.05'
        params += ' --angularSampling %0.3f' % self.angularSampling.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %d' % self.numOrientations.get()
        return params
    
    def phantomProject(self,volName):
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()         
        pathParticles = self._getPath('input_particles.xmd')
        Nx,Ny,Nz = self.inputParticles.get().getDim()
        R = -int(Nx/2)

        f = open(self._getExtraPath('params'),'w') 
        if (self.phaseFlipped.get()):
            doPhaseFlip = 1
        else:
            doPhaseFlip = 0  
        print doPhaseFlip             
        f.write("""# XMIPP_STAR_1 *
#
data_block1
_dimensions2D '%d %d'
_projAngleFile %s
_ctfPhaseFlipped %d
_applyShift 0
_noisePixelLevel   '0 0'""" % (Nx, Ny, pathParticles,self.phaseFlipped.get()))
        f.close()
        param =     ' -i %s' % volName
        param +=    ' --params %s' % self._getExtraPath('params')
        param +=    ' -o %s' % self._getPath('reference_particles.xmd')
        param +=    ' --sampling_rate % 0.3f' % self.inputParticles.get().getSamplingRate()
                
        #while (~isfile(self._getExtraPath('params'))):
        #    print 'No created'
        
        self.runJob('xmipp_phantom_project', 
                    param, numberOfMpi=1,numberOfThreads=1)

        param =     ' -i %s' % self._getPath('reference_particles.stk')
        param +=    ' --mask circular %d' % R
        self.runJob('xmipp_transform_mask', 
                    param, numberOfMpi=1,numberOfThreads=1)
                
    def numberOfProjections(self,sym,angularSampling):
        #http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/HowManyReferenceProjections
        dict = {0.5: 165016, 1.0: 41258, 1.5: 18338, 2.0: 10318, 2.5: 6600,
                3.0: 4586, 3.5: 3367,4.0: 2586, 4.5: 2042, 5.0: 1652, 5.5: 1367,
                6.0: 1148, 6.5: 977, 7.0: 843, 7.5: 732, 8.0: 643, 8.5: 569,
                9.0: 510, 9.5: 460, 10.0: 412, 11.0: 340, 12.0: 288, 13.0: 245,
                14.0: 211, 15.0: 184
                }
    
    def significantStep(self, volName, volDir, anglesPath, params):
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        makePath(volDir)  
        params += '  --initvolumes %s' % volName  
        params += ' --odir %s' % volDir
        params += ' --iter %d' % 1
        self.runJob('xmipp_reconstruct_significant', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        copyfile(volDir+'/angles_iter001_00.xmd', self._getExtraPath(anglesPath))
        
    def validationStep(self, volName,volDir,sym):
        makePath(volDir)  
        aFile = self._getExtraPath('exp_particles.xmd')
        aFileRef =self._getExtraPath('ref_particles.xmd')
        params = '  --volume %s' % volName  
        params += '  --angles_file %s' % aFile
        params += '  --angles_file_ref %s' % aFileRef
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
        
        if self.doNotUseWeights:
            params += ' --dontUseWeights'
            
        self.runJob('xmipp_multireference_aligneability',
                    params, numberOfMpi=1, numberOfThreads=1)
        
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        imgSet = self.inputParticles.get()
        for i, vol in enumerate(self._iterInputVols()):
            volume = vol.clone()               
            volDir = self._getVolDir(i+1)
            volPrefix = 'vol%03d_' % (i+1)
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            outImgSet = self._createSetOfParticles(volPrefix)
            
            outImgSet.copyInfo(imgSet)

            outImgSet.copyItems(imgSet,
                                updateItemCallback=self._setWeight,
                                itemDataIterator=md.iterRows(clusterMd,
                                                             sortByLabel=md.MDL_ITEM_ID))
                        
            mdValidatoin = md.MetaData(validationMd)
            weight = mdValidatoin.getValue(md.MDL_WEIGHT, mdValidatoin.firstObject())
            volume.weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set            
            outputVols.append(volume)
            self._defineOutputs(outputParticles=outImgSet)
        
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        #self._defineTransformRelation(self.inputVolumes.get(), volume)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolumes.get() and not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs

    def _summary(self):
        summary = ["Input particles:  %s" % self.inputParticles.get().getNameId()]
        summary.append("-----------------")

        if self.inputVolumes.get():
            for i, vol in enumerate(self._iterInputVols()):
                summary.append("Input volume(s)_%d: [%s]" % (i+1,vol))
        summary.append("-----------------")

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            for i, vol in enumerate(self._iterInputVols()):
                VolPrefix = 'vol%03d_' % (i+1)
                mdVal = md.MetaData(self._getExtraPath(VolPrefix+'validation.xmd'))                
                weight = mdVal.getValue(md.MDL_WEIGHT, mdVal.firstObject())
                summary.append("Output volume(s)_%d : %s" % (i+1,self.outputVolumes.getNameId()))
                summary.append("Quality parameter_%d : %f" % (i+1,weight))
                summary.append("-----------------")        
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using '
                            'the approach [Vargas2014a] with angular sampling '
                            'of %f and number of orientations of %f' % (self.angularSampling.get(),
                                                                        self.numOrientations.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getTmpPath('vol%03d' % volIndex)
    
    def _iterInputVols(self):
        """ In this function we will encapsulate the logic
        to iterate through the input volumes.
        This give the flexibility of having Volumes, SetOfVolumes or 
        a combination of them as input and the protocol code
        remain the same.
        """
        inputVols = self.inputVolumes.get()
        
        if isinstance(inputVols, Volume):
            yield inputVols
        else:
            for vol in inputVols:
                yield vol
        
    def _defineMetadataRootName(self, mdrootname,volId):
        
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn
    
    def _setWeight(self, item, row):  
        item._xmipp_weightClusterability = Float(row.getValue(md.MDL_VOLUME_SCORE1))
