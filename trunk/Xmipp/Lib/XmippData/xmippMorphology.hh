/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Pedro A. de Alarc�n     (pedro@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or   
 * (at your option) any later version.                                 
 *                                                                     
 * This program is distributed in the hope that it will be useful,     
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 * GNU General Public License for more details.                        
 *                                                                     
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
#ifndef _XMIPP_MORPHOLOGY_HH
#  define _XMIPP_MORPHOLOGY_HH

#include "xmippMatrices2D.hh"
#include "xmippMatrices3D.hh"

/**@name Mathematical morphology */
//@{
/**@name 2D processes.
    The neighbourhood must be 4 or 8.
    
    Count is the number of pixels
    meeting the condition so that the operation is performed, by default,
    0. Ie, if more than 0 pixels meet the condition then the corresponding
    operation is applied.
    
    Size is the size of the structuring element (box).
    
    The output image must be already resized to the desired shape*/
//@{
/** Dilate.
    See the group documentation for the parameter meanings */
void dilate2D(matrix2D<double> &in,matrix2D<double> &out, int neig,
   int count, int size);
/** Erode.
    See the group documentation for the parameter meanings */
void erode2D(matrix2D<double> &in,matrix2D<double> &out, int neig,
   int count, int size);
/** Closing=Dilation+Erosion */
void closing2D(matrix2D<double> &in,matrix2D<double> &out, int neig,
   int count, int size);
/** Opening=Erosion+Dilation */
void opening2D(matrix2D<double> &in,matrix2D<double> &out, int neig,
   int count, int size);
//@}

/**@name 3D processes.
    The neighbourhood must be 6, 18 or 26.
    
    Count is the number of voxels
    meeting the condition so that the operation is performed, by default,
    0. Ie, if more than 0 voxels meet the condition then the corresponding
    operation is applied.
    
    Size is the size of the structuring element (box).

    The output image must be already resized to the desired shape*/
//@{
/** Dilate.
    See the group documentation for the parameter meanings */
void dilate3D(matrix3D<double> &in,matrix3D<double> &out, int neig,
   int count, int size);
/** Erode.
    See the group documentation for the parameter meanings */
void erode3D(matrix3D<double> &in,matrix3D<double> &out, int neig,
   int count, int size);
/** Closing=Dilation+Erosion */
void closing3D(matrix3D<double> &in,matrix3D<double> &out, int neig,
   int count, int size);
/** Opening=Erosion+Dilation */
void opening3D(matrix3D<double> &in,matrix3D<double> &out, int neig,
   int count, int size);
//@}
//@}
#endif
