/*
	rwSPIDER.h
	Header file for reading and writing SPIDER files
	Format: 3D image file format for the SPIDER package
	Author: Bernard Heymann
	Created: 19990410 	Modified: 20010928
*/

#ifndef RWSPIDER_H
#define RWSPIDER_H

#define SPIDERSIZE 1024	// Minimum size of the SPIDER header (variable)
struct SPIDERhead {         	// file header for SPIDER data
    	float nslice;   	//  0	    	slices in volume (image = 1)
    	float nrow;   	    	//  1	    	rows per slice
    	float irec;   	    	//  2	    	# records in file (unused)
    	float nhistrec;	    	//  3	    	(obsolete)
    	float iform;  	    	//  4	    	file type specifier
    	float imami;  	    	//  5	    	max/min flag (=1 if calculated)
    	float fmax; 	    	//  6	    	maximum
    	float fmin; 	    	//  7	    	minimum
    	float av;   	    	//  8	    	average
    	float sig;  	    	//  9	    	standard deviation (=-1 if not calculated)
    	float ihist;  	    	// 10	    	(obsolete)
    	float nsam;   	    	// 11	    	pixels per row
    	float labrec;	    	// 12	    	# records in header
    	float iangle; 	    	// 13	    	flag: tilt angles filled
    	float phi;  	    	// 14	    	tilt angles
    	float theta;	    	// 15
    	float gamma;	    	// 16	    	(=psi)
    	float xoff; 	    	// 17	    	translation
    	float yoff; 	    	// 18
    	float zoff; 	    	// 19
    	float scale;	    	// 20	    	scaling
    	float labbyt; 	    	// 21	    	# bytes in header
    	float lenbyt; 	    	// 22	    	record length in bytes (row length)
    	float istack; 	    	// 23	    	indicates stack of images
    	float inuse;  	    	// 24	    	indicates this image in stack is used (not used)
    	float maxim;  	    	// 25	    	max image in stack used
    	float imgnum;  	    	// 26	    	number of current image
    	float unused[2];	// 27-28    	(unused)
    	float kangle; 	    	// 29	    	flag: additional angles set
    	float phi1; 	    	// 30	    	additional angles
    	float theta1;	    	// 31
    	float psi1; 	    	// 32
    	float phi2; 	    	// 33
    	float theta2;	    	// 34
    	float psi2; 	    	// 35

        double fGeo_matrix[3][3]; // x9 = 72 bytes: Geometric info
        float fAngle1; // angle info

        float fr1;
        float fr2; // lift up cosine mask parameters

        /** Fraga 23/05/97  For Radon transforms **/
        float RTflag; // 1=RT, 2=FFT(RT)
        float Astart;
        float Aend;
        float Ainc;
        float Rsigma; // 4*7 = 28 bytes
        float Tstart;
        float Tend;
        float Tinc; // 4*3 = 12, 12+28 = 40B

        /** Sjors Scheres 17/12/04 **/
        float Weight; // For Maximum-Likelihood refinement
        float Flip;   // 0=no flipping operation (false), 1=flipping (true)

        char fNada2[576]; // empty 700-76-40=624-40-8= 576 bytes

    	char cdat[12];  	// 211-213  	creation date
    	char ctim[8];		// 214-215  	creation time
    	char ctit[160]; 	// 216-255  	title
} ;


/************************************************************************
@Function: readSPIDER
@Description:
	Reading a SPIDER image file format.
@Algorithm:
	A 3D multi-image format used in electron microscopy.
	Header size:				1024 bytes (not same as data offset!).
	Data offset:				sizeof(float)*x_size*ceil(1024/x_size)
	File format extensions:  	.spi
	The identifier is a 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	File type and third dimension values
								must be less than 256*256.
	Data type: 					only float.
	Transform type: 			Hermitian
								The x-dimension contains the x-size
								of the full transform
	A multi-image file has a global header followed by a header and data
	for each sub-image.
@Arguments:
	Bimage*	p			the image structure.
	int img_select		image selection in multi-image file (-1 = all images).
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
int  readSPIDER(int img_select)
{
#define DEBUG
    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL ) return(-1);
	
    SPIDERhead*	header = (SPIDERhead *) balloc(sizeof(SPIDERhead));
    if ( fread( header, SPIDERSIZE, 1, fimg ) < 1 ) return(-2);
	
    // Determine byte order and swap bytes if from different-endian machine
    char*   	b = (char *) header;
    int     	i, j, swap = 0;
    int     	extent = SPIDERSIZE - 180;  // exclude char bytes from swapping
    if ( ( fabs(header->nslice) > SWAPTRIG ) || ( fabs(header->iform) > SWAPTRIG ) ||
         ( fabs(header->nslice) < 1 ) ) 
    {
#ifdef DEBUG
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif
    	swap = 1;
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
    
    // Map the parameters
    XSIZE(data) = px = (int) header->nsam;
    YSIZE(data) = py = (int) header->nrow;
    ZSIZE(data) = pz = (int) header->nslice;
    datatype = Float;
    transform = NoTransform;
    if ( header->iform < 0 ) 
    {
        transform = Hermitian;
        datatype = ComplexFloat;
        XSIZE(data) -= 2;
    }
    offset = (int) header->labbyt;
    min = header->fmin;
    max = header->fmax;
    avg = header->av;
    std = header->sig;
    ux = uy = uz = header->scale;
		
    size_t header_size = offset;
    size_t image_size = header_size + XSIZE(data)*YSIZE(data)*ZSIZE(data)*sizeof(float);
    size_t pad = offset;
    
    if ( header->istack > 0 ) 
        NSIZE(data) = (int) header->maxim;

    unsigned long 		imgstart = 0;
    unsigned long 		imgend = NSIZE(data);
    char*			hend;

    if ( img_select > -1 ) 
    {
        if ( img_select >= (long)NSIZE(data) ) 
            img_select = NSIZE(data) - 1;
        imgstart = img_select;
        imgend = img_select + 1;
        NSIZE(data) = 1;
        i = img_select;
    }
	
    image = new SubImage [NSIZE(data)];
    image->xoff = header->xoff;
    image->yoff = header->yoff;
    image->zoff = header->zoff;
    image->rot  = header->phi;
    image->tilt = header->theta;
    image->psi  = header->gamma;
    image->weight = header->Weight;
    image->flip = header->Flip;

    if ( header->istack > 0 ) {
        offset += offset;
        for ( i=imgstart; i<imgend; i++ ) {
            fseek( fimg, header_size + i*image_size, SEEK_SET );
            if ( fread( header, SPIDERSIZE, 1, fimg ) < 1 ) return(-3);
            hend = (char *) header + extent;
            if ( swap ) 
                for ( b = (char *) header; b<hend; b+=4 ) 
                    swapbytes(b, 4);
            j = ( pNSIZE(data) > 1 )? j = i: 0;
            image[j].xoff = header->xoff;
            image[j].yoff = header->yoff;
            image[j].zoff = header->zoff;
            image[j].rot  = header->phi;
            image[j].tilt = header->theta;
            image[j].psi  = header->gamma;
            image[j].weight = header->Weight;
            image[j].flip = header->Flip;
        }
    }

    bfree(header, sizeof(SPIDERhead));
	
#ifdef DEBUG
    std::cerr<<"DEBUG readSPIDER: header_size = "<<header_size<<" image_size = "<<image_size<<std::endl;
    std::cerr<<"DEBUG readSPIDER: img_select= "<<img_select<<" n= "<<NSIZE(data)<<" pad = "<<pad<<std::endl;
#endif
	
    readData(fimg, img_select, swap, pad );
    
    fclose(fimg);
		
    return(0);

}
/************************************************************************
@Function: writeSPIDER
@Description:
	Writing a SPIDER image file format.
@Algorithm:
	A 3D image format used in electron microscopy.
@Arguments:
	Bimage*				the image structure.
@Returns:
	int					error code (<0 means failure).
**************************************************************************/
/*
int 	writeSPIDER()
{
    if (typeid(T) == typeid(double))
    {
        float* fdata = (float *) balloc(datasize*sizeof(float));
        cast2Datatype
    }
    else if (typeid(T) == typeid(std::complex<double>))
    {

    }

	if ( datatype < ComplexShort ) img_to_float(p);
	else img_to_complex_float(p);
	
	if ( transform != NoTransform )
            img_convert_fourier(p, Hermitian);


	
#ifdef DEBUG
        printf("DEBUG writeSPIDER: Writing Spider file\n");
#endif            

    float		lenbyt = sizeof(float)*x;		// Record length (in bytes)
    float		labrec = floor(SPIDERSIZE/lenbyt);	// # header records
    if ( fmod(SPIDERSIZE,lenbyt) != 0 ) labrec++;
    float		labbyt = labrec*lenbyt; 		// Size of header in bytes
    offset = (int) labbyt;
	
    SPIDERhead*	header = (SPIDERhead *) balloc((int)labbyt*sizeof(char));
	
	// Map the parameters
    header->lenbyt = lenbyt;					// Record length (in bytes)
    header->labrec = labrec;					// # header records
    header->labbyt = labbyt; 					// Size of header in bytes

    header->irec = labrec + floor((x*y*z*sizeof(float))/lenbyt + 0.999999);	// Total # records
    header->nsam = x;
    header->nrow = y;
    header->nslice = z;

    // If a transform, then the physical storage in x is only half+1
    unsigned long		xstore = x;
    if ( transform == Hermitian ) {
        xstore = x/2 + 1;
        header->nsam = 2*xstore;
    }
	
#ifdef DEBUG
    printf("DEBUG writeSPIDER: Size: %g %g %g\n", header->nsam, header->nrow, header->nslice);
#endif
    if ( z < 2 ) 
    {
        if ( transform == NoTransform )
            header->iform = 1;				 // 2D image
        else
            header->iform = -12 + (int)header->nsam%2;   // 2D Fourier transform
    } else 
    {
        if ( transform == NoTransform )
            header->iform = 3;				 // 3D volume
        else
            header->iform = -22 + (int)header->nsam%2;   // 3D Fourier transform
    }
    header->imami = 1;
    header->fmin = min;
    header->fmax = max;
    header->av = avg;
    header->sig = std;
	
    // For multi-image files
    if (n > 1 ) 
    {
        header->istack = 2;
        header->inuse = -1;
        header->maxim = n;
    }
	
    header->xoff = 0.;
    header->yoff = 0.;
    header->zoff = 0.;
    header->phi = 0.;
    header->theta = 0.;
    header->gamma = 0.;

    unsigned long datatypesize = gettypesize(datatype);
    unsigned long datasize = xstore*y*z*datatypesize;
    unsigned long datasize_n = xstore*y*z;
	
#ifdef DEBUG
    printf("DEBUG writeSPIDER: Header size: %g\n", header->labbyt);
    printf("DEBUG writeSPIDER: Header records and record length: %g %g\n", header->labrec, header->lenbyt);
    printf("DEBUG writeSPIDER: Data type size: %ld\n", datatypesize);
    printf("DEBUG writeSPIDER: Data size: %ld\n", datasize);
    printf("DEBUG writeSPIDER: Data offset: %ld\n", offset);
    printf("DEBUG writeSPIDER: File %s\n", filename.c_str());
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "w") ) == NULL ) return(-1);
    fwrite( header, offset, 1, fimg );
	
    if ( verbose & VERB_DEBUG )
        printf("DEBUG writeSPIDER: File header written\n");
	
    float* fdata = (float *) balloc(datasize);
    if ( n == 1 ) {
        if ( dataflag ) 
        {
            castPage2Datatype(data, fdata, Float, datasize_n);
            fwrite( fdata, datasize, 1, fimg );
        }
    } else {
        header->istack = 0;
        header->inuse = 0;
        header->maxim = 0;
        for ( size_t i=0; i<n; i++ ) 
        {
            header->imgnum = i + 1;
            header->xoff = 0.;
            header->yoff = 0.;
            header->zoff = 0.;
            header->phi = 0.;
            header->theta = 0.;
            header->gamma = 0.;
            fwrite( header, offset, 1, fimg );
            castPage2Datatype(data + i*datasize, fdata, Float, datasize_n);
            fwrite( fdata, datasize, 1, fimg );
        }
    }
	
    fclose(fimg);
    bfree(fdata, datasize);
    bfree(header, (int)labbyt*sizeof(char));
	
    return(0);
}
*/
#endif

