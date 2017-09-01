#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std ;

#define D2R M_PI/180.
void uvec (double *in, float *out) ;

int main (int argc, char *argv[]) {
	int nextx, i, j, is, ns, nl,xloc, yloc ; 
	unsigned short attbyte = 0 ;
	unsigned int count=0 ;
	char ifname [420] , header[80];
	double x, y, z, lradius, degPerPixel, baseRad ;
	double lat, lon, val, *outdat, *xyzout ;
	float norm[3], triang[9], triang1[9] ;
	short *idata ;
	float  vexag, startlat ;
	float minx, maxx, miny, maxy, minz, maxz ;
	float minel, maxel, scalefac ;
	scalefac = 50. ;
	FILE *fin = fopen ("earth_dem_quartdeg", "r") ;


	for (i=0; i<80; i++) header[i] = 0;
	
	// we get min max of x y z for scaling to mm
	minx = 1.E9 ;
	miny = 1.E9 ;
	minz = 1.E9 ;
	maxx = -1.E9 ;
	maxy = -1.E9 ;
	maxz = -1.E9 ;
	minel = 1.E9 ;
	maxel = -1.E9 ;
	 

	ns = 1440 ;
	// for full globe nl /=1 else nl/=2
	nl = 720/2 ;
	baseRad = 6371. ;
	vexag = 17. ;

	xyzout = new double [3*ns * nl] ;
	idata = new  short [ns * nl] ;
	// for southern, comment out for northern
	//fseek (fin, ns*2*(nl-1), SEEK_SET) ; 
	fread ((char *) idata, 2, ns*nl, fin) ;
	fclose (fin) ;

//	open the stl file
	//FILE *fout =fopen("global_dem_s.stl", "w") ;
	FILE *fout =fopen("earth_dem_s_50.stl", "w") ;
	fwrite (header, 1, 80, fout) ;

	degPerPixel = 360. / ns ;
	startlat = 90. - (.5*degPerPixel) ;

// first go through dem and convert lat,lon, elev, to xyz
	cout << "converting lat lon to xyz" << endl ;
	for (i=0; i<nl; i+=1) {
		// mod for northern hem
		// lat = (startlat - (i * degPerPixel)) * D2R ;
		// for southern
		lat = (startlat - ((nl+i-1) * degPerPixel)) * D2R ;
		for (j=0; j<ns ; j+=1) {
			lon = (j * degPerPixel) * D2R ;
			val = idata[i*ns+j]/1000. ;
			lradius = val * vexag + baseRad ;
			x = lradius * cos(lat) * cos(lon) ;
			y = lradius * cos(lat) * sin(lon) ;
			z = lradius * sin (lat) ;
			x /= scalefac ;
			y /= scalefac ;
			z /= scalefac ;

			if (x > maxx) maxx = x ;
			if (y > maxx) maxy = y ;
			if (z > maxz) maxz = z ;
			if (x < minx) minx = x ;
			if (y < miny) miny = y ;
			if (z < minz) minz = z ;
			if (lradius < minel) minel= lradius  ;
			if (lradius > maxel) maxel= lradius  ;
			xyzout[i*ns*3+j*3]=x ;
			xyzout[i*ns*3+j*3+1]=y ;
			xyzout[i*ns*3+j*3+2]=z ;

		}
	}
	cout << "min and max x " << minx<< " " << maxx << endl ;
	cout << "min and max y " << miny<< " " << maxy << endl ;
	cout << "min and max z " << minz<< " " << maxz << endl ;
	cout << "min and max el " << minel<< " " << maxel << endl ;


	// now count the number of triangles
	count = (nl-1) * (ns) * 2 ;
	cout << "Number of triangles is " << count << endl ;
	// write the count to file
	fwrite ((char *) &count,4, 1,fout) ;

// create the triangles ;
	norm[0] =0 ;
	norm[1] =0 ;
	norm[2] =0 ;

	int numT = 0;
	for (i=0;i<nl-1; i++) {
		for (j=0; j<ns; j++) {
			//uvec (&xyzout[i*ns*3+j*3], &norm[0]) ;
			nextx = j + 1 ;
			if (nextx >= ns) nextx = 0 ;
			//if (j%2) {
			for (is=0; is<3; is++) triang[is]=xyzout[i*ns*3+j*3+is] ;
			for (is=0; is<3; is++) triang[3+is]=xyzout[(i)*ns*3+(nextx)*3+is] ;
			for (is=0; is<3; is++) triang[6+is]=xyzout[(i+1)*ns*3+j*3+is] ;
			fwrite ((char *) &norm[0], 4, 3, fout) ;
			fwrite ((char *) &triang[0], 4, 9, fout) ;
			fwrite ((char *)&attbyte, 2, 1, fout) ;
			numT++ ;
			
		//	} else {
			uvec (&xyzout[i*ns*3+nextx*3], &norm[0]) ;
			for (is=0; is<3; is++) triang[is]=xyzout[(i+1)*ns*3+j*3+is] ;
			for (is=0; is<3; is++) triang[is+3]=xyzout[(i+1)*ns*3+(nextx)*3+is] ;
			for (is=0; is<3; is++) triang[is+6]=xyzout[(i)*ns*3+(nextx)*3+is] ;
		//	}
			fwrite ((char *) &norm[0], 4, 3, fout) ;
			fwrite ((char *) &triang[0], 4, 9, fout) ;
			fwrite ((char *)&attbyte, 2, 1, fout) ;
			numT++ ;
				
		}	

	}
	cout << "Num tri written is "<< numT << endl ;
//FILE *fout =fopen ("globedat", "w") ;
	//fwrite (outdat,4,nptsTot, fout) ;
	fclose (fout) ;
	delete [] xyzout ;
	delete [] idata ;

}

			



void uvec (double *in, float *out) {
	int i ;
	double total=0 ;
	for (i=0;i<3; i++) 
		total += (in[i] * in[i]) ;

	total = sqrt(total) ;
	for (i=0; i<3; i++) 
		out[i] = in[i] / total ;
}

