/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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

/* Includes ---------------------------------------------------------------- */
#include "QtColorLabel.hh"
#include "qcolor.h"


/* Constructor ------------------------------------------------------------- */
QtColorLabel::QtColorLabel() {
   __col = new QColor[__numColors];
   __col[0]  = Qt::red;         
   __col[1]  = Qt::green;      
   __col[2]  = Qt::blue;
   __col[3]  = Qt::black;       
   __col[4]  = Qt::cyan;       
   __col[5]  = Qt::magenta;
   __col[6]  = Qt::yellow;      
   __col[7]  = Qt::darkRed;    
   __col[8]  = Qt::darkGreen;
   __col[9]  = Qt::darkBlue;    
   __col[10] = Qt::darkCyan;   
   __col[11] = Qt::darkMagenta;
   __col[12] = Qt::darkYellow;  
   __col[13] = Qt::darkGray;   
   __col[14] = Qt::lightGray;
}

QtColorLabel::~QtColorLabel() {
   delete __col;
}

/* Get the color for the i family ------------------------------------------ */
QColor QtColorLabel::col( int i ) {
   return( __col[ i % __numColors ] );
}
