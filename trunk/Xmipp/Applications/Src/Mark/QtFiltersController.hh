/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef __QT_FILTERS_CONTROLLER_HH__
#define __QT_FILTERS_CONTROLLER_HH__

/* Includes ---------------------------------------------------------------- */
#include <vector>
#include "qwidget.h"
#include <XmippData/xmippMicrograph.hh>

/* Forward declarations ---------------------------------------------------- */
class QImage;
class QtFilter;
class QDialog;
class QListBox;

/* Filters generic class --------------------------------------------------- */
class QtFiltersController : public QWidget {
   Q_OBJECT   
   
private:
   vector<QtFilter*>  __filterList;
   QDialog           *__addFilterDialog;
   QListBox          *__listFilters;
   const Micrograph  *__M;
   
   enum filters {
      invertContrastFilter,
      enhanceContrastFilter,
      substractBackgroundFilter,
      removeOutlierFilter,
      lowpassFilter,
      highpassFilter
   };
   
public:
   // Constructor
   QtFiltersController( QWidget * _parent, const Micrograph *_M );
   ~QtFiltersController();

   // Apply the filters list
   void applyFilters( QImage *_img );

public slots:
   void slotAddFilter();
   void slotAddFilter( int _f );
   void slotCleanFilters();
};

#endif
