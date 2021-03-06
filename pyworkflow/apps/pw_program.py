#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
Launch main project window 
"""

import sys
import pyworkflow.utils as pwutils



if __name__ == '__main__':

    program = sys.argv[1]
    params = ' '.join('"%s"' % x for x in sys.argv[2:])
    
    env = None
    
    if program.startswith('xmipp'):
        import pyworkflow.em.packages.xmipp3 as xmipp3
        env = xmipp3.getEnviron()
    if program.startswith('relion'):
        import pyworkflow.em.packages.relion as relion
        env = relion.getEnviron()        
    elif (program.startswith('e2') or 
          program.startswith('sx')):
        import pyworkflow.em.packages.eman2 as eman2
        env = eman2.getEnviron()
    elif program.startswith('b'):
        import pyworkflow.em.packages.bsoft as bsoft
        env = bsoft.getEnviron()
    
    pwutils.runJob(None, program, params, env=env)
         
