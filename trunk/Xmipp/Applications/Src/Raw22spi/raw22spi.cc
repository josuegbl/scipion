/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <XmippData/xmippTypes.hh>

/* Prototypes -============================================================= */
void Usage (char **argv);
void raw22spi (const FileName &, const FileName &, char, int, int, int); 

int main (int argc, char *argv[]) {
/* Input Parameters ======================================================== */
FileName       fn_in;    // input file
FileName       fn_out;   // output file
FileName       sel_file; // selection file
string 	       sel_ext;  // extension for output files in selection file.	
char	       raw_type = 'b';
int            Zdim, Ydim, Xdim;
int 	       size_arg;

Zdim=Ydim=Xdim=0;
/* Parameters ============================================================== */
   try {
       if (check_param(argc, argv, "-i")){
         fn_in = get_param(argc, argv, "-i");
	 if (check_param(argc, argv, "-sel")||check_param(argc, argv, "-oext")) {
         EXIT_ERROR(1,"Raw22spi: -i option is not compatible with -sel or -oext");
	 }
       } 
       if (check_param(argc, argv, "-o")){
         fn_out = get_param(argc, argv, "-o");
	 if (check_param(argc, argv, "-sel")||check_param(argc, argv, "-oext")) {
         EXIT_ERROR(1,"Raw22spi: -o option is not compatible with -sel or -oext");
	 }
       } 
       if (check_param(argc, argv, "-sel")){
         sel_file = get_param(argc, argv, "-sel");
	 sel_ext  = get_param(argc, argv, "-oext");
         if (check_param(argc, argv, "-i")||check_param(argc, argv, "-o")){ 
          /*error cause -sel is not compatible with -i -o*/ 
	 EXIT_ERROR(1,"Raw22spi: -sel option is not compatible with -i or -o");
         }
       }

       if (check_param(argc, argv, "-f")) raw_type='f'; 
       
       if (check_param(argc, argv, "-s")){ 
         size_arg = position_param(argc, argv, "-s");
         if (size_arg+3>=argc) EXIT_ERROR(1,"Not enough parameters behind -s");
         Zdim= AtoI(argv[size_arg+1]);
         Ydim= AtoI(argv[size_arg+2]);
         Xdim= AtoI(argv[size_arg+3]);
       }
       if (argc == 1) {Usage(argv);}
   }
   catch (Xmipp_error XE) {cout << XE; Usage(argv);}

try {
/* Perform conversion ====================================================== */

 /* input is a sel file*/
 if (check_param(argc, argv,"-sel") && check_param(argc, argv,"-oext")){
   SelFile SF(sel_file);
   FileName SF_out_name;
   SF_out_name = sel_file.without_extension().add_prefix("out_");
   SF_out_name += sel_file.get_extension();
   SelFile SF_out;
   SF_out.clear();
   while (!SF.eof()) {
      SelLine line= SF.current();
      if (line.Is_data()) { //The SelLine is not a comment
       FileName in_name = line.get_text();
       short label = line.get_label();
       FileName out_name = in_name.without_extension();
       out_name.add_extension(sel_ext);
       SF_out.insert(out_name,(SelLine::Label)label);
       raw22spi(in_name,out_name,raw_type,Zdim,Ydim,Xdim);
      }
      else if (line.Is_comment()){
        SF_out.insert(line);
      }
      SF.next();
   }  // while 
   SF_out.write(SF_out_name); //write output sel file
 }
 //input/output are single files
 
 else if (check_param(argc, argv,"-i") && check_param(argc, argv,"-o")){
   raw22spi(fn_in,fn_out,raw_type,Zdim,Ydim,Xdim);
 }
 
} catch (Xmipp_error XE) {cout << XE;}
   exit(0);
} //main

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage (char **argv) {
  printf (
     "Usage: %s [Purpose and Parameters]"
     "\nPurpose: Convert from a 2d/3d raw images to Xmipp ones (and viceversa)"
     "\n        Input/Output can be either a single file or a set of them "
     "\n        (specified in a 'sel' file)"           
     "\nParameter Values: (note space before value)"
     "\nI/O parameters"
     "\nESPECIFIC PARAMETERS FOR SINGLE-FILE CONVERSION"
     "\n    -i    file_in        input raw or Xmipp file"
     "\n    -o    file_out       output Xmipp or raw file"
     "\nESPECIFIC PARAMETERS FOR SEL-FILE CONVERSION"     
     "\n    -sel  input_file     input sel file"
     "\n    -oext input_file     extension for the output files if the input"
     "\n			 files were specified in a sel file"
     "\nGENERAL PARAMETERS"
     "\n   [-f]		         raw file is read/written in float format "
     "\n                         (byte by default)."      
     "\n    -s Zdim Ydim Xdim    Z,Y,X dimensions for input files."
     "\n			 For 2D raw images set the Zdim to 1\n"
     ,argv[0]);
}
void raw22spi (const FileName &fn_in, const FileName &fn_out,
     char raw_type, int Zdim,  int Ydim, int Xdim) {

 VolumeXmipp    Vx;
 ImageXmipp     Ix;

  // Volume Xmipp --> Raw Volume
   if (Is_VolumeXmipp(fn_in)) {
      Vx.read(fn_in);
      switch (raw_type) {
         case 'b': ((Volume)Vx).write(fn_out,Vx.reversed(),VBYTE); break;
         case 'f': ((Volume)Vx).write(fn_out,Vx.reversed(),VFLOAT); break;
      }
   // Image Xmipp --> Raw Image
   } else if (Is_ImageXmipp(fn_in)) {
     Ix.read(fn_in);
      switch (raw_type) {
         case 'b': ((Image)Ix).write(fn_out,Ix.reversed(),IBYTE); break;
         case 'f': ((Image)Ix).write(fn_out,Ix.reversed(),IFLOAT); break;
      }
   }   
   // Raw image --> Xmipp Image
   else if (Zdim==1) {
      Image *I=&Ix; // This is a trick for the compiler
      switch (raw_type) {
         case 'b': I->read(fn_in,Ydim,Xdim,FALSE,IBYTE); break;
         case 'f': I->read(fn_in,Ydim,Xdim,FALSE,IFLOAT); break;
      }
      Ix.write(fn_out);
   // Raw Volume --> Xmipp Volume
   } else {
      Volume *V=&Vx; // This is a trick for the compiler
      switch (raw_type) {
         case 'b': V->read(fn_in,Zdim,Ydim,Xdim,FALSE,VBYTE); break;
         case 'f': V->read(fn_in,Zdim,Ydim,Xdim,FALSE,VFLOAT); break;
      }
      Vx.write(fn_out);
   }
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Raw22spi {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Raw22spi/Help/raw22spi.html";
      help="Converts one or several images and volumes from Raw format to
            Xmipp and ivceversa.";
      OPEN MENU menu_raw22spi;
      COMMAND LINES {
	+ single: xmipp_raw22spi -i $FILE_IN -o $FILE_OUT
                      [-f] [-s $ZDIM $YDIM $XDIM]
        + multiple: xmipp_raw22spi -sel $SELFILE_IN -osel $OEXT
                      OPT(-f) OPT(-s)
      } 
      PARAMETER DEFINITIONS {
        $FILE_IN {
	   label="Input File";
           help="Raw or Xmipp";
	   type=file existing;
	}
        $FILE_OUT {
	   label="Output File";
           help="Xmipp or Raw";
	   type=file;
	}
        OPT(-f) {
           label="The Raw file is a float raw";
           help="Otherwise, it is a byte raw";
        }
        $SELFILE_IN {
	   label="Input SelFile";
           help="All images to be converted";
	   type=file existing;
	}
        $OEXT {
	   label="Output Extension";
           help="This extension substitutes the one in the input files";
	   type=text;
	}
        OPT(-s) {
           label="Raw dimensions";
           help="Only needed when input raw";
        }
           $ZDIM {
              label="Z Dimension";
              type=natural;
           }
           $YDIM {
              label="Y Dimension";
              type=natural;
           }
           $XDIM {
              label="X Dimension";
              type=natural;
           }
      }
   }

   MENU menu_raw22spi {
      "Single conversion"
      $FILE_IN
      $FILE_OUT
      "Multiple conversions"
      $SELFILE_IN
      $OEXT
      "Size specifications"
      OPT(-f)
      OPT(-s)
   }
*/
