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
#include "QtWidgetMicrograph.hh"
#include "QtFilterMenu.hh"
#include "QtImageMicrograph.hh"
#include "QtImageOverviewMicrograph.hh"
#include <XmippData/xmippMicrograph.hh>
#include <XmippData/xmippArgs.hh>

/* Constructor ------------------------------------------------------------- */
QtWidgetMicrograph::QtWidgetMicrograph( QtMainWidgetMark *_mainWidget, 
                                        QtFiltersController *_f, 
                                        Micrograph *_m ) : 
   QWidget( (QWidget*) _mainWidget ) {
   __filtersController = _f;
   __m              = NULL;
   __activeFamily   = -1;
   __tilted         = FALSE;

   __gridLayout     = new QVBoxLayout( this );
   __menuBar        = new QMenuBar( this );
      __menuBar->setSeparator( QMenuBar::InWindowsStyle );
   __mImage         = new QtImageMicrograph( 0 );
   __mImageOverview = new QtImageOverviewMicrograph( this );
   __file_menu      = NULL;
   
   __mImage->setFiltersController( _f );
   __mImageOverview->setFiltersController( _f );
        
   connect( __mImageOverview, SIGNAL(signalSetCoords(int, int)),
            __mImage, SLOT(slotSetCoords(int, int)) );
   connect( __mImage, SIGNAL(signalSetCoords(int, int)),
            __mImageOverview, SLOT(slotSetCoords(int, int)) );   
   connect( __mImage, SIGNAL(signalSetWidthHeight(int, int)),
            __mImageOverview, SLOT(slotSetWidthHeight(int, int)) );
   connect( __mImage, SIGNAL(signalRepaint( void )),
            __mImageOverview, SLOT(slotRepaint( void )) );
   connect( __mImageOverview, SIGNAL(signalRepaint( void )),
            __mImage, SLOT(slotRepaint( void )) );
   connect( __mImage, SIGNAL(signalRepaint( void )),
            __mImage, SLOT(slotRepaint( void )) );
   connect( __mImageOverview, SIGNAL(signalRepaint( void )),
            __mImageOverview, SLOT(slotRepaint( void )) );
   connect( __mImage, SIGNAL(signalAddCoordOther( int, int, int )),
            this, SLOT(slotDrawEllipse( int,int,int )) );
      
   connect( this, SIGNAL(signalActiveFamily(int)),
            __mImage, SLOT(slotActiveFamily(int)) );
   connect( this, SIGNAL(signalActiveFamily(int)),
            __mImageOverview, SLOT(slotActiveFamily(int)) );
   
   QAccel *ctrl = new QAccel( this );
   ctrl->connectItem( ctrl->insertItem(Key_G+CTRL, 200),
                             this, SLOT(slotChangeContrast(void)) );

   setMicrograph( _m );
   
   __mImage->show();
   __gridLayout->setMenuBar( __menuBar );
   __gridLayout->addWidget( __mImageOverview );
   
   openMenus();   
}

QtWidgetMicrograph::~QtWidgetMicrograph() {
   delete __mImage;
   delete __mImageOverview;
   delete __menuBar;
   delete __gridLayout;
}

/* Set Micrograph ---------------------------------------------------------- */
void QtWidgetMicrograph::setMicrograph( Micrograph *_m ) {
   if ( _m != NULL ) {
      __m = _m;
      __mImage->setMicrograph( _m );
      __mImageOverview->setMicrograph( _m );
   }   
}

/* Open menus -------------------------------------------------------------- */
void QtWidgetMicrograph::openMenus() {
   __file_menu = new QtFileMenu( this );
   connect( __mImage, SIGNAL(signalAddCoordOther(int,int,int)),
            __file_menu, SLOT(slotCoordChange()) );   

   QtFilterMenu *filterMenu = new QtFilterMenu( this );
   
   addMenuItem( "&File", (QtPopupMenuMark *)(__file_menu) );
   addMenuItem( "F&ilters", (QtPopupMenuMark *)(filterMenu) );
   
   connect( __file_menu, SIGNAL(signalAddFamily(const char *)),
            this, SLOT(slotAddFamily(const char*)) );
   
   connect( (QObject*)filterMenu, SIGNAL(signalAddFilter()),
            (QObject*)__filtersController, SLOT(slotAddFilter()) );   
   connect( (QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            (QObject*)__filtersController, SLOT(slotCleanFilters()) );
   connect( (QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            this, SLOT(slotRepaint()) );

   // *** Add your own menus
}

void QtWidgetMicrograph::changeContrast(int _mingray, int _maxgray, float _gamma) {
   __mingray=_mingray;
   __maxgray=_maxgray;
   __gamma  =_gamma;
   __mImage->changeContrast(_mingray,_maxgray,_gamma);
   __mImageOverview->changeContrast(_mingray,_maxgray,_gamma);
}

void QtWidgetMicrograph::repaint( int t ) {
   __mImage->repaint( FALSE );
   __mImageOverview->repaint( FALSE );
}

void QtWidgetMicrograph::slotDrawEllipse(int _x, int _y, int _f) {
   __mImage->drawEllipse(_x,_y,_f);
   __mImageOverview->drawEllipse(_x,_y,_f);
}

void QtWidgetMicrograph::slotDrawLastEllipse(int _x, int _y, int _f) {
   __mImage->drawEllipse(_x,_y,_f);
   __mImageOverview->drawEllipse(_x,_y,_f);
   __mImage->drawLastEllipse(_x,_y,_f);
}

/* Active family ----------------------------------------------------------- */
void QtWidgetMicrograph::slotActiveFamily( int _f ) {
   __activeFamily = _f;
   emit signalActiveFamily( _f );
}

void QtWidgetMicrograph::slotAddFamily( const char *_familyName ) {
   emit signalAddFamily( _familyName );
}

void QtWidgetMicrograph::slotDeleteMarkOther( int _coord ) {
   __m->coord(_coord).valid = false;
   repaint();
}

void QtWidgetMicrograph::slotChangeFamilyOther( int _coord, int _f ) {
   __m->coord(_coord).label = _f;
   repaint();
}

void QtWidgetMicrograph::slotQuit() {
   __file_menu->slotQuit();
}

void QtWidgetMicrograph::slotChangeContrast() {
   AdjustContrastWidget *adjustContrast=new
      AdjustContrastWidget(0,255,1.0F,this,
        0,"new window", WDestructiveClose);
   adjustContrast->show();
}

/* AdjustContrastWidget ---------------------------------------------------- */
// Constructor
AdjustContrastWidget::AdjustContrastWidget(int min, int max, float gamma,
   QtWidgetMicrograph *_qtwidgetmicrograph,
   QWidget *parent, const char *name, int wflags):
   QWidget(parent,name,wflags) {
   __qtwidgetmicrograph=_qtwidgetmicrograph;

    // Set this window caption
    setCaption( "Adjust Contrast" );

    // Create a layout to position the widgets
    QBoxLayout *Layout = new QVBoxLayout( this, 10 );

    // Create a grid layout to hold most of the widgets
    QGridLayout *grid = new QGridLayout( 3,3 );
    Layout->addLayout( grid, 5 );

    // Minimum
    QLabel     *label_min= new QLabel( this, "label" );    
    label_min->setFont( QFont("times",12,QFont::Bold) );
    label_min->setText( "Minimum" );
    label_min->setFixedSize(label_min->sizeHint());
    grid->addWidget( label_min, 0, 0, AlignCenter );

    __scroll_min= new QScrollBar(0, 255, 1, 1, min,
       QScrollBar::Horizontal,this,"scroll");
    __scroll_min->setFixedWidth(100); 
    __scroll_min->setFixedHeight(15); 
    grid->addWidget( __scroll_min, 0, 1, AlignCenter );
    connect( __scroll_min, SIGNAL(valueChanged(int)),
       SLOT(scrollValueChanged(int)) );

    __label_min= new QLabel( this, "label" );	 
    __label_min->setFont( QFont("courier",14) );
    __label_min->setText( ItoA(min,3) );
    __label_min->setFixedSize(__label_min->sizeHint());
    grid->addWidget( __label_min, 0, 2, AlignCenter );

    // Maximum
    QLabel     *label_max= new QLabel( this, "label" );    
    label_max->setFont( QFont("times",12,QFont::Bold) );
    label_max->setText( "Maximum" );
    label_max->setFixedSize(label_max->sizeHint());
    grid->addWidget( label_max, 1, 0, AlignCenter );

    __scroll_max= new QScrollBar(0, 255, 1, 1, max,
       QScrollBar::Horizontal,this,"scroll");
    __scroll_max->setFixedWidth(100); 
    __scroll_max->setFixedHeight(15); 
    grid->addWidget( __scroll_max, 1, 1, AlignCenter );
    connect( __scroll_max, SIGNAL(valueChanged(int)),
       SLOT(scrollValueChanged(int)) );

    __label_max= new QLabel( this, "label" );	 
    __label_max->setFont( QFont("courier",14) );
    __label_max->setText( ItoA(max,3) );
    __label_max->setFixedSize(__label_max->sizeHint());
    grid->addWidget( __label_max, 1, 2, AlignCenter );

    // Gamma
    QLabel     *label_gamma= new QLabel( this, "label" );    
    label_gamma->setFont( QFont("times",12,QFont::Bold) );
    label_gamma->setText( "Gamma" );
    label_gamma->setFixedSize(label_gamma->sizeHint());
    grid->addWidget( label_gamma, 2, 0, AlignCenter );

    __scroll_gamma= new QScrollBar(0, 40, 1, 1, (int)(10*gamma),
       QScrollBar::Horizontal,this,"scroll");
    __scroll_gamma->setFixedWidth(100); 
    __scroll_gamma->setFixedHeight(15); 
    grid->addWidget( __scroll_gamma, 2, 1, AlignCenter );
    connect( __scroll_gamma, SIGNAL(valueChanged(int)),
       SLOT(scrollValueChanged(int)) );

    __label_gamma= new QLabel( this, "label" );    
    __label_gamma->setFont( QFont("courier",14) );
    __label_gamma->setText( FtoA(gamma,3,2) );
    __label_gamma->setFixedSize(__label_gamma->sizeHint());
    grid->addWidget( __label_gamma, 2, 2, AlignCenter );

}

// One of the sliders changed ----------------------------------------------
void AdjustContrastWidget::scrollValueChanged(int new_val) {
    __label_min  ->setText( ItoA(__scroll_min  ->value(),3) );
    __label_max  ->setText( ItoA(__scroll_max  ->value(),3) );
    __label_gamma->setText( FtoA((__scroll_gamma->value())/10.0,3,2) );
    __qtwidgetmicrograph->changeContrast(__scroll_min->value(),
       __scroll_max->value(),__scroll_gamma->value()/10.0);
}
