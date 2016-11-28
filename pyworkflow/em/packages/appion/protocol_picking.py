# **************************************************************************
# *
# * Authors:     Laura del Cano (ldelcano@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, basename

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.em.convert import ImageHandler
from pyworkflow.utils.properties import Message

from convert import readSetOfCoordinates



class DogPickerProtPicking(ProtParticlePicking):
    """ Pick particles in a set of micrographs using Appion dogpicker. """

    _label = 'dogpicker'
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):

        ProtParticlePicking._defineParams(self, form)
        form.addParam('diameter', params.IntParam, default=100,
                   label='Diameter of particle')
        form.addParam('invert', params.BooleanParam, default=False,
                   label='Invert',
                      help = "Invert image before picking, DoG normally picks "
                             "white particles.")
        form.addParam('threshold', params.FloatParam, default=0.5,
                      label='Threshold',
                      help = "Threshold in standard deviations above the mean, "
                             "e.g. --thresh=0.7")
        form.addParam('extraParams', params.StringParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='Additional parameters for dogpicker: \n'
                      '  --numberSizes, --sizeRange, --threshold,...')

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        inputMics = self.getInputMicrographs()
        sampling = inputMics.getSamplingRate()

        self._params = {}
        # diameter must be passed in Angstrongs and therefore should be
        # converted to pixels
        self._params['apix'] = sampling
        self._params['diam'] = self.diameter.get() * sampling
        self._params['thresh'] = self.threshold.get()

        args = ""
        for par, val in self._params.iteritems():
            args += " --%s=%s" % (par, str(val))

        if self.invert:
            args += " --invert"

        args += " " + self.extraParams.get('')

        # Store all steps ids, final step createOutput depends on all of them
        deps = []
        ih = ImageHandler()

        for mic in self.inputMicrographs.get():
            # Create micrograph folder
            micName = mic.getFileName()
            micDir = self._getTmpPath(pwutils.removeBaseExt(micName))
            pwutils.makePath(micDir)

            # If needed convert micrograph to mrc format, otherwise link it
            if pwutils.getExt(micName) != ".mrc":
                fnMicBase = pwutils.replaceBaseExt(micName, 'mrc')
                inputMic = join(micDir, fnMicBase)
                ih.convert(mic.getLocation(), inputMic)
            else:
                inputMic = join(micDir, basename(micName))
                pwutils.createLink(micName, inputMic)

            # Insert step to execute program
            stepId = self._insertFunctionStep('executeDogpickerStep',
                                              inputMic, args)
            deps.append(stepId)


        self._insertFinalSteps(deps)


    def _insertFinalSteps(self, deps):
        # Insert step to create output objects
        self._insertFunctionStep('createOutputStep', prerequisites=deps)


    #--------------------------- STEPS functions -------------------------------
    def executeDogpickerStep(self, inputMic, args):

        # Program to execute and it arguments
        program = "ApDogPicker.py"
        outputFile = self._getExtraPath(pwutils.replaceBaseExt(inputMic, "txt"))

        args += " --image=%s --outfile=%s" % (inputMic, outputFile)

        # Run the command with formatted parameters

        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)


    def createOutputStep(self):
        coordSet = self._createSetOfCoordinates(self.getInputMicrographs())
        self.readSetOfCoordinates(self._getExtraPath(), coordSet)
        coordSet.setBoxSize(self.diameter.get())
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(self.inputMicrographs, coordSet)

    
    #--------------------------- UTILS functions -------------------------------
    def getFiles(self):
        return self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMicrographs.get(), coordSet)

    def _summary(self):
        summary = []
        summary.append("Number of input micrographs: %d"
                       % self.getInputMicrographs().getSize())

        if self.getOutputsSize() > 0:
            summary.append("Number of particles picked: %d"
                           % self.getCoords().getSize())
            summary.append("Particle size: %d" % self.getCoords().getBoxSize())
            summary.append("Threshold: %0.2f" % self.threshold)
            if self.extraParams.hasValue():
                summary.append("And other parameters: %s" % self.extraParams)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)

        return summary

    def _citations(self):
        return ['Voss2009']

    def _methods(self):
        methodsMsgs = []

        if self.getInputMicrographs() is None:
            return ['Input micrographs not available yet.']

        methodsMsgs.append("Input micrographs %s of size %d."
                           % (self.getObjectTag(self.getInputMicrographs()),
                              self.getInputMicrographs().getSize()))

        if self.getOutputsSize() > 0:
            output = self.getCoords()
            methodsMsgs.append('%s: User picked %d particles with a particle '
                               'size of %d and threshold %0.2f.'
                               % (self.getObjectTag(output), output.getSize(),
                                  output.getBoxSize(), self.threshold))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs