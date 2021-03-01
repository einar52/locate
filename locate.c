/*  © Copyright  2008, 2021 Einar Kjartansson */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "locate.h"
#include "sgeom.h"

/*********
Parameters for fit:
origin time
lat
lon
m 
*********/

#define ERadius 6382.75 /* average earth radius of curvature at 65 N lat */
#define VEL 6.58

typedef struct { double x,y ; } Proj ;
typedef struct {
	Vec3 pos ;
	double  pgv,f0,fBTime ;
	char name[8] ;
	double x,y,trav ;
	double mag,wmag ;
	int status, band, nPick ;
	double pgvP,f0P ;
} StaEvent ;
StaEvent *data ;

typedef struct {
	short sta1,sta2,sta3, u ;
	float x,y,t,l1Sum ;
} StaTriplet ;
StaTriplet minTri ;
double magnitude ,sigmaMag,medianF0,medianF0P,locQ ;
char locStatus ;

time_t uBaseTime ;
int event,date,year,nSta,nLog ;
int locMode  = 3 ;
char month[10],timeHMS[20] ;
Vec3 nPole, pCenter ;
Mat33 projection, inverseProjection ;
double *aMatrix, *bVector ;
FILE *outFile ;

Proj v32ProjX( Vec3 v )
{	
/* Azimuthal equidistant projection */
	double r , azi ;
	Proj res ;
	azi = v3Vertex(nPole,pCenter,v) ;
	r = ERadius * v3Angle2(pCenter,v) ;
	res.x = r * sin(azi)  ;
	res.y = r * cos(azi)  ;
	return res ;
}
Vec3 proj2V3X( Proj pp )
/* inverse Azimuthal equidistant projection */
{	
	double angle, azi ;
	angle = sqrt( pp.x * pp.x + pp.y * pp.y) /ERadius ;
	azi = atan2(pp.x,pp.y) ;
	return v3DistAzi(pCenter,angle,azi) ;
}
Proj v32Proj( Vec3 v )
/* orthograpic projection */
{
	Vec3 pv ;
	Proj p ;
	pv = v3MMult(&projection,v) ;
	p.x = pv.x ;
	p.y = pv.y ;
	return p ;
}
Vec3 proj2V3( Proj pp )
/* inverse orthograpic projection */
{
	double z2 ;
	Vec3  pv ;
	z2 = ERadius*ERadius - pp.x*pp.x - pp.y*pp.y ;
	pv.x = pp.x ;
	pv.y = pp.y ;
	pv.z = sqrt(z2) ;
	return  v3MMult( &inverseProjection,pv) ;
}
void xy2LatLon( double x, double y, double *lat, double *lon)
{
        LatLon ll ;
        Vec3 vv ;
        Proj pp ;
        pp.x = x ; pp.y = y ;
        vv = proj2V3(pp) ;
        ll = v32LatLon(vv) ;
        *lon = ll.lon ; *lat = ll.lat ;
}

double locateLinear( int *staList, int n, double *x )
/*	compute location using a linear traveltime ( constant vel ) 
	staList lists indecies into data, n is number of stations to use.
	n  >= 4
*/
#define T0  0.63
{
	double *a,*b ;
	double xi,xj,yi,yj,ti,tj ;
	double sum,d,r ;
	Proj loc ;
	Vec3 v ;
	int i,j ;
	StaEvent *si,*sj ;
	a = aMatrix ;
	b = bVector ;
	for( i = 0 ; i < n ; i++ ) {
		j = (i+1) % n ;	
		si = data + staList[i] ;
		sj = data + staList[j] ;
		xi = si->x ; xj = sj->x ;
		yi = si->y ; yj = sj->y ;
		ti = VEL * si->fBTime ; tj = VEL * sj->fBTime  ;
		a[i]      = 2.0*(xj - xi) ;
		a[i+n]    = 2.0*(yj - yi) ;
		a[i+2*n]  = 2.0*(ti - tj) ;
		b[i] = xj*xj - xi*xi + yj*yj - yi*yi + ti*ti - tj*tj ;
	}
/*	for( i = 0 ; i < n ; i++) 
		fprintf(outFile,"%10g %10g %10g    %10g\n",a[i],a[i+n],a[i+2*n],b[i]) ; */
	golubC(a,x,b,n,3) ;
/*	fprintf(outFile,"x = [%10g %10g %10g]\n",x[0],x[1],x[2] ) ; */
	for( i = 0 ; i < n ; i++) {
		si = data + staList[i] ;
		xi = si-> x - x[0] ; yi = si->y - x[1] ; ti=si->fBTime ;
		r = sqrt( xi*xi + yi*yi)  ;
		tj =  (r+x[2]) /VEL ;
		d = ti - tj ;
		sum += d*d ;
		fprintf(outFile,"%2d %2d %s %8.3f %8.3f %8.3f %8.3f\n",i,staList[i],si->name,r,ti,tj,d) ;
	}
	loc.x = *x ;
	loc.y = x[1] ;
	v = proj2V3(loc) ;
	v3Print(outFile,v,"location ") ;
	x[2] /= VEL ;
	return sqrt(sum/n-2) ;
}
void printIntArr(char *str, int *a, int n)
{
	int i ;
	if( NULL == outFile) return ;
	fprintf(outFile,"%s :",str) ;
	for( i = 0 ; i < n ; i++) fprintf(outFile," %2d",a[i]) ;
	fprintf(outFile,"\n") ;
		
}
int compDouble(const void *p1, const void *p2 )
{
	double d ;
	d = *(double*) p2 - *(double*)p1 ;
	if( d < 0.0 ) return 1 ;
	if( d > 0.0 ) return -1 ;
	return 0 ;
}
double testTeleseis( double *x )
{
	int i,n,j ;
	StaEvent *p ;
	double *a, *b , tt,d,sum,sigma,apparentV ;
	n = 0 ;
	for( i = 0 ; i < nSta ; i++)  if( data[i].status == 0 ) n++ ;
	if( n < 4 ) return(0) ;
	a = aMatrix ; b = bVector ;
	j = 0 ;
	for( i = 0 ; i < nSta ; i++) {
		p = data + i ;	
		if(p->status) continue ;
		a[j] =   p->x ;
		a[j+n] = p->y ;
		a[j+2*n] = 1.0 ;
		b[j] = p->fBTime ;
		j++ ;
	}
	golubC(a,x,b,n,3) ;
	apparentV = 1.0/sqrt(x[0] * x[0] + x[1]*x[1]) ;
	fprintf(outFile,"Test Teleseis:  [%12.7f %12.7f %12.7f]\n", x[0], x[1], x[2] ) ;
	sum = 0.0 ;	
	for( i = 0 ; i < nSta ; i++) {
		p = data + i ;	
		if(p->status) continue ;
		tt = p->x*x[0] + p->y*x[1] + x[2] ;
		d = p->fBTime - tt ;
		sum += d*d ;
		fprintf(outFile,"%s %7.3f %7.3f %7.3f\n",p->name,p->fBTime,tt,d) ;
	}
	sigma = sqrt(sum/(n-3) ) ;
	fprintf(outFile," sigma = %6.3f   Vel =%8.3f\n",sigma,apparentV) ;
	return sigma ;
}
void spatialCheck()
{
	int i,n,nn,j,ij,iTest ;
	int *s1,*s2,*ss ;
	StaEvent **sa ;
	double *dist,d,*dist2, test ;
	n = 0 ;
	for( i = 0 ; i < nSta ; i++)  if( data[i].status == 0 ) n++ ;
	if( n < 4 ) return ;
	nn = ( n*n - n ) / 2 ;
	sa = (StaEvent **) calloc(n,sizeof(*sa)) ;
	s1 = (int*) calloc(nn,sizeof(int)) ;
	s2 = (int*) calloc(nn,sizeof(int)) ;
	ss = (int*) calloc(n, sizeof(int)) ;
	dist  = (double*) calloc(nn,sizeof(double)) ;
	dist2 = (double*) calloc(nn,sizeof(double)) ;
	j = 0 ;
	for( i = 0 ; i < nSta ; i++) if( data[i].status == 0 ) sa[j++] = data+i ;
	ij = 0 ;
	for( i = 0 ; i < n - 1 ; i++ ) {
		for( j = i+1 ; j < n ; j++ ) {
			d = ERadius*v3Angle2( sa[i]->pos,sa[j]->pos ) ; 
			if(outFile) fprintf(outFile,"%s %s %6.1f\n",sa[i]->name,sa[j]->name,d) ;
			s1[ij] = i   ; s2[ij] = j ;
			dist[ij] = d ; dist2[ij++] = d ;
		}
	}
	qsort(dist2,nn,sizeof(double),compDouble) ;
	if( outFile) {
		for(j= 0 ; j < nn ; j++) fprintf(outFile,"%5.1f ",dist2[j] ) ;
		fprintf(outFile,"\n") ;
	}
	iTest = nn/2-1 ;
	test = 3.0*dist2[iTest] ;  
	if(outFile) fprintf(outFile,"iTest = %d  test = %5.1f\n",iTest,test) ;
	for( i = 0 ; i < nn ; i++) if( dist[i] > test ) {
		ss[s1[i]]++ ; 
		ss[s2[i]]++ ;
	}
	j = 0 ; for( i = 1 ; i < n ; i++) if( ss[i] >= ss[j] ) j=i ;
	while( ss[j] > 0 ) {
		printIntArr("ss ",ss,n); printIntArr("s1 ",s1,nn); printIntArr("s2 ",s2,nn);
		sa[j]->status = 100 ;
		if(outFile) fprintf(outFile,"Spatial %s\n",sa[j]->name) ;
		for( i = 0 ; i < nn ; i++) {
			if(( s1[i] == j) || ( s2[i] == j )) {
				ss[s1[i]]-- ;
				ss[s2[i]]-- ;
			}
		}
		j = 0 ; for( i = 1 ; i < n ; i++) if( ss[i] >= ss[j] ) j=i ;
	}
	free(dist2) ; free(dist) ; free(ss) ; free(s2) ; free(s1) ; free(sa) ;
}
void markImpossible()
{
	int n,nn,i,j,ij,*s1,*s2,*ss ;
	StaEvent *p1,*p2 ;
	nn = (nSta*nSta-nSta)/2 ;
	double dist,dt,dtm ;
	s1 = (int*) calloc(nn,sizeof(int)) ;
	s2 = (int*) calloc(nn,sizeof(int)) ;
	ss = (int*) calloc(nSta,sizeof(*ss)) ;
	ij = 0 ;
	for( i = 0 ; i < nSta ; i++) {
		p1 = data+i ; 
		if( p1->nPick < 3 ) p1->status = 1 ;
	}
	for( i = 0 ; i < nSta-1 ; i++) {
		p1 = data + i  ;
		if(outFile) fprintf(outFile," %s %d %d %d\n",p1->name,i,p1->status,p1->nPick) ;
		if( p1->status == 0 ) for( j = i+1 ; j < nSta ; j++ ) {
			p2 = data + j ;
			if( p2->status == 0 ) {
				dtm = ERadius * v3Angle2(p1->pos,p2->pos)/VEL ;
				dt = p1->fBTime - p2->fBTime ;
				if(outFile) fprintf(outFile,"%s %s %8.2f %8.2f\n",p1->name,p2->name,dt,dtm ) ; 
				if( (fabs(dt) - dtm) > 0.5) {
					ss[i]+=2 ; ss[j]+=3 ;
					s1[ij] = i ;
					s2[ij++] = j ;
				}
			}
		}
	}
	nn = ij ;
	printIntArr( "ss" ,ss,nSta); printIntArr( "s1" ,s1,nn); printIntArr( "s2" ,s2,nn);
	j = 0 ; 
	for( i = 1 ; i < nSta ; i++)  if( ss[i] >= ss[j] ) j = i ; 
	while( ss[j] > 1 ) { 
		p1 = j + data ;
		p1-> status = ss[j]+1 ;
		if(outFile) fprintf(outFile,"%s\n",p1->name) ;
		ij = 0 ;
		while (ij < nn ) {
			if( (s1[ij] == j ) || ( s2[ij] == j )) {
				ss[s1[ij]] -= 2 ;
				ss[s2[ij]] -= 3 ;
				nn-- ;
				s1[ij] = s1[nn] ;
				s2[ij] = s2[nn] ;
			}
			else ij++ ;
		}	
	printIntArr( "ss" ,ss,nSta); printIntArr( "s1" ,s1,nn); printIntArr( "s2" ,s2,nn);
		j = 0 ;
		for( i = 1 ; i < nSta ; i++)  if( ss[i] >= ss[j] ) j = i ; 
	}
	free(ss) ; free(s2); free(s1) ;
}
void printX(char *text, double *x, int m )
{
	int i ;
	if( NULL == outFile) return ;
	fprintf(outFile,"%s [",text) ;
	for( i = 0 ; i < m ; i++) fprintf(outFile,"%12.7f",x[i]) ;
	fprintf(outFile,"]\n") ;
}
double locateNL( int *staList, int n0, int m, double *x) 
{
#define DELTA 0.01 
	double *a,*b ;
	double xd[4] ;
	double xi,xj,yi,yj,ti,t0,tj ;
	double sum,d,r,dt,w,rr,rMax,dMax ;
	int i,j,n,iMax ;
	StaEvent *si,*sj ;
	n = n0 ;
	if( m == 4 ) n = n0*2 ;
	a = aMatrix ;
	b = bVector ;
	for( i = 0 ; i < n0 ; i++ ) {
		si = data + staList[i] ;
		xi = x[0] - si->x  ; 
		yi = x[1] - si->y  ; 
		ti = si->fBTime -x[2] ; 
		r = sqrt( xi*xi + yi*yi ) ;
		t0 = travelTime(r) ;
		dt = (travelTime(r+ DELTA )-t0)/ DELTA ;
		w = si->pgv ;
		w = sqrt(w) ;
		a[i]     = xi*w*dt/r ;
		a[i+n]   = yi*w*dt/r ;
		a[i+2*n] = w ;
		b[i] = w*(ti-t0)  ; 
		if( m == 4 ) {
			a[i+7*n0] = w ;
			a[i+n0] =  a[i]/(2.3*r) ;
			a[i+3*n0] =  a[i+n]/(2.3*r) ;
			b[i+n0] = w*(log10(si->pgv) - x[3] - 1.12 + 1.63*log10(r)) ;
		}
	}
	if (outFile && (m == 4) ) for( i = 0 ; i < n ; i++) 
		fprintf(outFile,"%12.5f %12.5f %12.5f %12.5f  %12.5f\n",a[i],a[i+n],a[i+2*n],a[i+3*n],b[i]) ;
	golubC(a,xd,b,n,m) ;
	printX("xd = ",xd,m) ;
	r = sqrt(xd[0]*xd[0] + xd[1]*xd[1] ) ;
	rMax = 10.0 ;
	if( r > rMax ) for( i = 0 ; i < m ; i++) xd[i] *= rMax/r ;
	for( i = 0 ; i < m ; i++) x[i] += xd[i] ;
	printX("x  = ",x,m) ;
	return r ;
}
void checkLoc( int *staList, int n, int m, double *x, double deltaT) 
{
	double *a,*b ;
	double xd[4] ;
	double xi,xj,yi,yj,ti,t0,tj ;
	double sum,d,dt,w,rr,rMax,dMax ;
	int i,j,iMax ;
	StaEvent *si,*sj ;
	Proj loc ;
	Vec3 v ;
	if( outFile == NULL) return ;
	dMax = 0.0 ;
	sum = 0.0 ;
	for( i = 0 ; i < n ; i++) {
		si = data + staList[i] ;
		xi = si-> x - x[0] ; yi = si->y - x[1] ; ti=si->fBTime ;
		rr = sqrt( xi*xi + yi*yi)  ;
		tj = travelTime(rr) + x[2] ;
		d = ti - tj ;
		if( fabs(d) > dMax ) { dMax = fabs(d) ; iMax = i ; }
		sum += d*d ;
		fprintf(outFile,"%2d %2d %s %8.3f %8.3f %8.3f %8.3f\n",i,staList[i],si->name,rr,ti,tj,d) ;
	}
	if(( n > 4 ) && (dMax > deltaT )) data[staList[iMax]].status= 100 ;
	loc.x = *x ;
	loc.y = x[1] ;
	v = proj2V3(loc) ;
	v3Print(outFile,v,"location ") ;
	fprintf(outFile,"n =%2d sd =%8.2f\n",n,sqrt(sum/(n-3))) ;
}
int makeStaList(int *list) 
{
	int i,n ;
	n = 0 ;
	for( i = 0 ; i < nSta ; i++) {
		if( data[i].status == 0  ) {
			list[n++] = i ;
		}
	}
	return n ;
}
int compMagI( const void *p1, const void *p2 ) 
{
	StaEvent *s1, *s2 ;
	double d ;
	s1 = data + *(int *)p1 ; 
	s2 = data + *(int *)p2 ;
	d = s2->mag - s1->mag ;
        if( d < 0.0 ) return 1 ;
        if( d > 0.0 ) return -1 ;
}
double magWeight( double mag, double r )
{
	double logPGV,w,den,x  ;
	logPGV = mag + 1.12 - 1.63*log10(r) ;
	if( logPGV > 4.0 ) return 0.0 ;
	x = 4.0 - logPGV ;
	den = x + 0.4  ;
	w = 16.0*x/(den*den) ;
	w = w / ( 1.0 + 0.01*r ) ;
/*	if( r > 150.0 ) w = 0 ; */
	return w ;
}
double doMagnitude(double *x, int n , int *list)
{
	int *index ; 
	int i,j ;
	StaEvent *si ;
	double xi,yi,r,w,sumw,sw ;
	index = (int * ) calloc(nSta,sizeof(*index)) ;
	for( j = 0 ; j < 3 ; j++) {
	   sumw = 0.0 ;
	   for( i = 0 ; i < nSta ; i++ ) {
		index[i] = i ;
		si = data + i ;
		xi = x[0] - si->x ;
		yi = x[1] - si->y ;
		r = sqrt( xi*xi + yi*yi ) ;
		si->mag = log10(si->pgv) - 1.12 + 1.63*log10(r) ;
		if(j)  w = magWeight(magnitude,r) ;
		else w = 1.0 ;
		si->wmag = w ;
		sumw += w ;
		if( outFile) fprintf(outFile,"%s %6.1f %9.1f %6.2f %6.2f\n",si->name,r,si->pgv,si->mag,si->wmag) ; 
	   }
	   qsort(index,nSta,sizeof(int),compMagI) ;
	   sw = sumw *0.5 ;
	   i = 0 ;
           while( sw > 0.0 ) sw -= data[index[i++]].wmag  ;
	   magnitude  = data[index[i]].mag ;
	   if(outFile) fprintf(outFile,"mag = %8.2f\n",magnitude) ;
	}
/*	for( i = 0 ; i < nSta ; i++) { fprintf(outFile, " %7.2f",data[index[i]].mag ) ; } */
	
	free( index ) ;
}
void setInitX( double *x, int m, int *list ) 
{
	StaEvent *si,*sj ;
	x[2] = 0.0 ; x[m-1] = 0.0 ;
	si = data+list[1] ;
	sj = data+list[2] ;
	x[0] = 0.5*( si->x + sj->x ) ;
	x[1] = 0.5*( si->y + sj->y ) ;
}
double distd(double x, double y, double r)
{
	return sqrt(x*x + y*y) - sqrt(r*r) ;  
	return x*x + y*y - r*r ; 
}
int solveQuad(double dx, double dy, double x, double y, double r, double *r1, double *r2)
{	/* solve quadradic equation for location based on 3 stations - 2 solutions */
	double a,b,c,ba,d,dd ;
	a = 1.0 - dx*dx - dy*dy  ; 
	b = r + dx*x + dy*y  ;
	c = r*r - x*x - y*y ;
	dd = b*b - c*a ;
	if( fabs(a) < 0.01 ) return 0 ;
	if( dd <= 0.0 ) return 0  ;
	ba = b/a ;
	d = sqrt(dd)/a ;
	*r1 = ba - d ;
	*r2 = ba + d ;
	if(outFile) fprintf(outFile,"a=%g b=%g ba=%g c=%g d=%g r1=%g r2=%g\n", a,b,ba,c,d,*r1,*r2) ;
	return 2 ;
}
int locate3(int i1, int i2, int i3, double *xv1, double *xv2 ) 
{
	StaEvent *s1,*s2,*s3 ;
	double x1,x2,x3,y1,y2,y3,r1,r2,r3,rMin ;
	double a1,a2,b1,b2,c1,c2,cc1,cc2,det,det1 ;
	double x0,y0,dx,dy ;
	int i,status ;
	double x,y,r,d1,d2,d3,rr1,rr2 ;
	s1 = data + i1 ;
	s2 = data + i2 ;
	s3 = data + i3 ;
	x1 = s1->x ; x2 = s2->x ; x3 = s3->x ;
	y1 = s1->y ; y2 = s2->y ; y3 = s3->y ;
	r1 = VEL* s1->fBTime ; r2 = VEL*s2->fBTime ; r3 = VEL*s3->fBTime ;
	rMin = r1 ; 
	if( rMin < r2 ) rMin = r2 ;
	if( rMin < r3 ) rMin = r3 ;
	a1 = 2.0*(x2-x1) ;
	b1 = 2.0*(y2-y1) ;
	a2 = 2.0*(x3-x2) ;
	b2 = 2.0*(y3-y2) ;
	c1 = x1*x1 - x2*x2 + y1*y1 - y2*y2 - r1*r1 + r2*r2 ;
	c2 = x2*x2 - x3*x3 + y2*y2 - y3*y3 - r2*r2 + r3*r3 ;
	cc1 = c1 + 2.0*(r1-r2) ;
	cc2 = c2 + 2.0*(r2-r3) ;
	det = a1*b2 - a2*b1 ;
	if( fabs(det) < 0.01 ) return 0 ;
	det1 = 1.0/det ;
	x0 = ( c2*b1 - c1*b2 ) * det1 ;
	y0 = ( c1*a2 - c2*a1 ) * det1 ;
	dx = ( cc2*b1 - cc1*b2 ) * det1 - x0 ;
	dy = ( cc1*a2 - cc2*a1 ) * det1 - y0 ;
	status = solveQuad(dx,dy,x0-x1,y0-y1,r1,&rr1,&rr2) ; 
	if(outFile) fprintf(outFile,"det=%g dx=%6.2f dy=%6.2f len = %6.2f\n",det,dx,dy,sqrt(dx*dx+dy*dy)) ;
	if( status == 0 ) return 0 ;
	if( rr1 >= rMin ) return 0 ;
	xv1[0] = x0 + rr1*dx ; xv1[1] = y0 + rr1*dy ; xv1[2] = rr1/ VEL ;
	if( rr2 >= rMin ) return 1 ;
	xv2[0] = x0 + rr2*dx ; xv2[1] = y0 + rr2*dy ; xv2[2] = rr2/ VEL ;
	return 2 ;
	for( i = 0 ; i < 50 ; i++) {
		r = i*0.1+14 ;
		x = x0 + r*dx ;
		y = y0 + r*dy ;
		d1 = distd(x1-x,y1-y,r1-r) ;
		d2 = distd(x2-x,y2-y,r2-r) ;
		d3 = distd(x3-x,y3-y,r3-r) ;
		fprintf(outFile,"%3d %6.2f %6.2f %6.2f %8.2f %8.2f %8.2f\n",i,x,y,r,d1,d2,d3) ;
	}
}
double distSta( int i, double x, double y )
{
	double xx,yy,d ;
	xx = data[i].x -x ;
	yy = data[i].y -y ;
	d = sqrt(xx*xx+yy*yy) ;
	if(outFile) fprintf(outFile,"%s %6.1f  ",data[i].name,d) ;
	return d ;
}
int loc3Test(int s1, int s2, int s3, double *x, StaTriplet *minTri, int *staList,int n )
{
	int list[3],i ;
	double rr,ra,sum,sumw, xi,yi,ti,tj,w ;
	StaEvent *si ;
	list[0] = s1 ; list[1] = s2 ; list[2] = s3 ;
	rr = locateNL(list,3,3,x) ;
	rr = locateNL(list,3,3,x) ;
	checkLoc(list,3,3,x,100) ;
	if( rr > 1.0 ) return 0 ;	/* no or slow convergence */
	sum = 0.0 ; sumw = 0.0 ;
	for( i = 0 ; i < n ; i++) {
		si = data + staList[i] ;
		xi = si-> x - x[0] ; yi = si->y - x[1] ; 
		ti=si->fBTime ;
		ra = sqrt( xi*xi + yi*yi)  ;
		tj = travelTime(ra) + x[2] ;
		w = 50.0/ra ;
		sum += w * fabs(tj-ti) ;
		sumw += w ;
	}
	sum = sum / sumw ; 
	if(outFile) fprintf(outFile,"sum= %7.1f  sumw= %8.0f\n",sum,sumw) ; 
	minTri->u++ ;
	if( sum < minTri->l1Sum ) {
		if(outFile) fprintf(outFile,"sum=%g minTri->l1Sum=%g\n",sum,minTri->l1Sum) ;
		minTri->sta1 = s1 ; minTri->sta2 = s2 ; minTri->sta3 = s3 ; 
		minTri->x = x[0] ; minTri->y = x[1] ; minTri->t = x[2] ;
		minTri->l1Sum = sum ;
	}
	return 1 ;
}
int locateInit3( int *list, int n, double *x )
{
	int i1,i2,i3,s1,s2,s3,nn,nSol ;
	double xa[3],xb[3] ;
	minTri.l1Sum = 1.e8 ;
	nSol = 0 ;
	for( i1 = 0 ; i1 < n-2 ; i1++) {
	  for(i2 = i1+1; i2 < n-1 ; i2++) {
	    for(i3 = i2+1; i3 < n  ; i3++ ) {
		s1 = list[i1] ; s2 = list[i2] ; s3 = list[i3] ;
		if(outFile) fprintf(outFile,"nSol= %d %s %s %s\n", 
			nSol,data[s1].name,data[s2].name,data[s3].name) ;
		nn = locate3(s1,s2,s3,xa,xb) ;
		if( nn > 0 ) nSol += loc3Test(s1,s2,s3,xa,&minTri,list,n ) ;
		if( nn > 1 ) nSol += loc3Test(s1,s2,s3,xb,&minTri,list,n ) ;
	    }
	  }
	}
	if( nSol == 0 ) return 0 ;
	x[0] = minTri.x ;
	x[1] = minTri.y ;
	x[2] = minTri.t ;
	if(outFile) fprintf(outFile,"FINAL nSol=%d %s %s %s\n", nSol,data[minTri.sta1].name, 
		data[minTri.sta2].name, data[minTri.sta3].name) ;
	return nSol ;
}
int locate(int flags)
{
	int list[200],n,nSol ;
	static double x[3] ;
	double r ;
	if( 2 & flags ) markImpossible() ;
	spatialCheck() ;
	n = makeStaList(list) ;
	if( n < 4 ) return 0 ;
	nSol = locateInit3(list,n,x) ;
	if( nSol < 1 ) return 0 ;
	checkLoc(list,n,3,x,200.0) ;
	r = doMagnitude(x,n,list) ;
	return n ;
}
void locateOld()
{
	int list[200],i,j,n,nn,m ;
	static double x[4] ;
	double d,zl,r ;
	m = 3 ;
	j = 0 ;
	markImpossible() ; 
	spatialCheck() ;
	testTeleseis(x) ;
	n = makeStaList(list) ;
	if( n < 4 ) return ;
/*	zl = locateLinear(list,n,x) ;  */
	if ( locMode == 0 ) setInitX(x,m,list) ;
	if ( locMode == 3 ) { 
		locateInit3(list,n,x) ;
	}
	do {
		checkLoc(list,n,m,x,5.0 ) ;
		nn = n ;
		n = makeStaList(list) ;
	} while( nn != n )  ;
	j = 0 ;
	d = 11.0 ;
	while( (d > 10.0) && ( j < 20)) {
		d = locateNL(list,n,m,x) ; j++ ;
		checkLoc(list,n,m,x,4.0 ) ;
		n = makeStaList(list) ;
		}
	d = locateNL(list,n,m,x) ;
	checkLoc(list,n,m,x,3.0 ) ;
	n = makeStaList(list) ;
	j = 0 ;
	while( (d > 0.001) && ( j < 10 ) ) {
		n = makeStaList(list) ;
		checkLoc(list,n,m,x,3.0 ) ;
		d = locateNL(list,n,m,x) ;
		j++ ;
	}
	r = doMagnitude(x,n,list) ;
}
void initProjection( Vec3 center )
{
	nPole.z = 1.0 ;
	makeRotation(&projection,center,nPole,ERadius) ;
	transposeMat33(&inverseProjection,&projection,
		1.0/(ERadius*ERadius)) ;
}
void fixCenter()
{
	int i ;
	double xx,yy,zz,sx,sy,w	 ;
	Proj pp ;
	Vec3 *p,v ;
	v.x = 0.0 ; v.y = 0.0 ; v.z = 0.0 ;
	for( i = 0 ; i < nSta ; i++) {
		p = &(data[i].pos) ;
		w = data[i].pgv ;
		w = sqrt(w) ; 	/* weigh by sqrt of peak ground velocity */
		v.x += w* p->x ; v.y += w* p->y ; v.z += w* p->z ; 
	}
	pCenter = v3Normalize(v) ;
	initProjection(pCenter) ;
	sx = 0.0 ; sy = 0.0 ;
	for( i = 0 ; i < nSta ; i++) {
		p = &(data[i].pos) ;
		pp =  v32Proj( *p ) ;
/*		v3Print(outFile,*p," *p ") ;
	fprintf(outFile,"%s %8.2f %8.2f \n",data[i].name, pp.x,pp.y) ; */
		data[i].x = pp.x ;
		data[i].y = pp.y ;
		sx += pp.x ; sy += pp.y ;
	}
	if(outFile) fprintf(outFile,"sx,sy = %8g, %8g\n",sx/nSta, sy/nSta ) ;
}


int a2month( char *name) 
{                   /* hash routine to convert month name to number */
	static int val[] = 
{0,0,0,3,0,0,4,11,0,10,0,9,0,12,0,0,6,0,0,0,0,0,0,0,0,5,1,0,0,7,8,0,0,0,0,0,0,0,0,0,0,0,0,0,2 } ;
	int *ip ;
	ip = ( int *) name ;
	return val[*ip % 45] ;
}
void setTime() 
{	
	static struct tm t ;
	t.tm_sec = 10*(timeHMS[6]-'0') + timeHMS[7] - '0' ;
	t.tm_min = 10*(timeHMS[3]-'0') + timeHMS[4] - '0' ;
	t.tm_hour = 10*(timeHMS[0]-'0') + timeHMS[1] - '0' ;
	t.tm_mday = date ;
	t.tm_mon =  a2month(month)-1 ;
	t.tm_year = year - 1900 ;
	uBaseTime = mktime(&t) ;
/*printf(" %d %d %d %d %d %d %s\n",t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,uBaseTime,ctime(&uBaseTime)) ; */
}
void readEvent( char *file )
{
	char b1[200],b2[20],b3[20],b4[20],b5[20],b6[20],b7[20],b8[20],b9[20] ;
	int j1,j2,j3,j4,j5 ;
	int ns,i ;
	FILE *fd ;
	double lat,lon ;
	StaEvent  *s ;
	fd = fopen(file,"r") ;
	if(outFile) fprintf(outFile," read %s fd=%x\n",file,fd) ;
	ns = fscanf(fd,"%d %d %d %d %s %s",&event,&j1,&j2,&j3,b1,month ) ;
	ns = fscanf(fd,"%d %s %d %d %d ",&date,timeHMS,&year,&nLog,&nSta ) ;
	if(outFile) fprintf(outFile,"ns=%d month=%s nSta=%d\n",ns,month,nSta) ;
	data = ( StaEvent*) calloc(nSta,sizeof(StaEvent)) ;
	for( i = 0 ; i < nSta ; i++) {
		s = data+i ;
		ns=fscanf(fd,"%lf %lf %s %lf %s %lf %d %d %s %d %s %s %s",
			&lon,&lat,b2,&(s->pgv),s->name,&(s->f0),
			&(s->band),&j3,b4,&(s->nPick),b5,b6,b7) ;
		s->pos = v3New(lat,lon) ;
		s->fBTime = j3*0.01 ;
		if( *b6 != '-') {
			s->f0P = atof(b6) ;
			s->pgvP = atof(b7) ;
		} else { s->f0P = 0.0 ; s->pgvP = 0.0 ; }
		
	}
	setTime() ;
	fixCenter() ;
	aMatrix = (double *) calloc(8*nSta,sizeof(aMatrix)) ;
	bVector = (double *) calloc(2*nSta,sizeof(bVector)) ;
}
void printData()
{
	int i ;
	StaEvent *s ;
	char b1[40] ;
	if( outFile == NULL) return ;
	for( i = 0 ; i < nSta ; i++) {
		s = data + i ;
		fprintf(outFile,"%s %8g %9g %6g %3d (%10.3f %10.3f )",
			s->name,s->fBTime,s->pgv,s->f0,s->nPick,s->x,s->y) ;
		v3Print(outFile,s->pos," ") ;
	}
	v3Print(outFile,pCenter,"pCenter =") ;
}
void formatTest(char *ss, int test, double val) 
{
	if ( test ) sprintf(ss,"%6.2f",val) ;
	else sprintf(ss,"%6s","- ") ;
	return ;
}
void getQuality()
{
	int i,n ;
	StaEvent *si ;
	double xi,yi,rr,ti,tj,res,magSum,dm, *fArr ,*fPArr ;
	int n3,n10,n1,np ;
	fArr = ( double *) calloc(nSta,sizeof(fArr)) ;
	fPArr = ( double *) calloc(nSta,sizeof(fPArr)) ;
	n3 = 0 ; n1 = 0 ; n10 = 0 ;
	magSum = 0 ; 
	n = 0 ;
	np = 0 ;
	for( i = 0 ; i < nSta ; i++ ) {
		si = data + i ;
		xi = si->x -minTri.x ; yi = si->y -minTri.y ; 
		rr = sqrt(xi*xi+yi*yi) ;
		ti = si->fBTime - minTri.t ;
		tj = travelTime(rr) ;
		res = fabs(ti-tj) ;
		if( si->status == 0 ) {
			if( res < 0.1 ) n1++ ;
			if( res < 0.3 ) n3++ ;
			if( res < 1.0 ) n10++ ;
			dm = si->mag - magnitude ;
			magSum += dm*dm ;
			fArr[n] = si->f0 ;
			n++ ;
			if( si->f0P > 0.0 )  fPArr[np++] = si->f0P ; 
		}
	}
	qsort(fArr,n,sizeof(double),compDouble) ; 
	medianF0  = 0.5* (fArr[n/2] + fArr[(n-1)/2]) ;
	sigmaMag = sqrt(magSum/n) ;
	if(np > 1 ) {
		qsort(fPArr,np,sizeof(double),compDouble) ;
		medianF0P = 0.5* (fPArr[np/2] + fPArr[(np-1)/2]) ;
	}
/*
	locQ = (1.0-1.0/(n3-1.5))*(1.0-1.0/(n1-1.5))*(1.0-1.0/(n10-1.4)) ;
	locQ /=  0.7 + 2.0*sigmaMag ;
	locQ = 10.0*sqrt(locQ) ;
*/
	locStatus = 'C' ;
	if( n1> 4 ) locStatus = 'B' ;
	if(( n3 > 5 ) && (n3 > nSta/5)) locStatus = 'B' ;
	if( n1 > 6 ) locStatus = 'A' ;
	if( n10 < 4 ) locStatus = 'D' ;
	if( n10+3 < nSta/2 ) locStatus = 'D' ;
	if( sigmaMag > 0.6 ) locStatus = 'M' ;
	if( (medianF0 < 0.77) && (locStatus > 'B' )) locStatus = 'T' ;
	n1 = n1 - 3 ; if( n1 > 9 ) n1 = 9 ;
	n3 = n3 - 3 ; if( n3 > 9 ) n3 = 9 ;
	n10 = n10 - 3 ; if( n10 > 9 ) n10 = 9 ;
	locQ = n1+0.1*n3+0.01*n10 ;
	free(fArr) ;
}
void printReport(int list)
{
	double lat,lon ;
	char timeStr[30] ;
	time_t uT ;
	int ss,i ;
	double xi,yi,ti,rr,tj, fTest ;
	char s1[16],s2[16],s3[16],s4[16];
	StaEvent *si ;
	xy2LatLon(minTri.x, minTri.y, &lat,&lon) ;
	uT = uBaseTime - 3600 + ( minTri.t + 3600.0 ) ;
	strftime(timeStr,30,"%Y%m%d %H%M%S",gmtime(&uT)) ;
	ss = (minTri.t + 3600.0) * 10 ;
	ss = ss % 10 ;
	fTest = 2.0*log10(medianF0)+magnitude ;
	if( minTri.u <= 2 )  return ;
  	printf("_%03d %s.%1d %7.3f %8.3f 4.000 %4.2f ",event,timeStr,ss,lat,lon,magnitude) ;
	printf("%c %03.0f %5.2f %5.2f %4.1f",locStatus,100.0*locQ,sigmaMag,medianF0,fTest ) ;
	printf(" %2d %5.2f %3d %2d\n",nSta,medianF0P,minTri.u,minTri.sta3) ;
/*		ada      1.4  160.9  2.63   0 -12.30  24.94 -37.24 */
	if( list == 0 ) return ;
	printf("\n         pgv   dist   mag   f0   st fBreak     t     dt\n") ;
	for( i = 0 ; i < nSta ; i++) {
		si = data + i ;
		xi = si->x -minTri.x ; yi = si->y -minTri.y ; 
		rr = sqrt(xi*xi+yi*yi) ;
		ti = si->fBTime - minTri.t ;
		tj = travelTime(rr) ;
		sprintf(s1,"%6.2f",ti) ;
		sprintf(s2,"%6.2f",tj) ;
		sprintf(s3,"%6.2f",ti-tj) ;
		formatTest(s1,1-si->status,ti) ;
		formatTest(s2,1,tj) ;
		formatTest(s3,0==si->status,ti-tj) ;
		printf("%s %8.1f %6.1f %5.2f %5.2f %3d",si->name,si->pgv,rr,si->mag,si->f0,si->status) ;
		printf(" %s %s %s\n",s1,s2,s3) ;
		
	}
}
void testProjection()
{
	Vec3 c,a,b,d ;
	Proj pa,pb,pc ;
	initProjection(v3New(65,-19)) ;
	a = v3New(65,-20) ; pa = v32Proj(a) ;
	fprintf(outFile,"%10.4f %10.4f ",pa.x,pa.y) ;
	v3Print(outFile,a," a ") ;
	a = v3New(65,-21) ; pa = v32Proj(a) ;
	fprintf(outFile,"%10.4f %10.4f ",pa.x,pa.y) ;
	v3Print(outFile,a," a ") ;
	d = proj2V3(pa) ; v3Print(outFile,d," d ") ;
	a = v3New(66,-19) ; pa = v32Proj(a) ;
	fprintf(outFile,"%10.4f %10.4f ",pa.x,pa.y) ;
	v3Print(outFile,a," a ") ;
	d = proj2V3(pa) ; v3Print(outFile,d," d ") ;
	
}
void testt()
{
	char b1[40],b2[40],b3[40] ;
	Mat33 a,b ;
	Vec3 vv,v2,v3 ;
	return ;
	pCenter = v3New(65,-19) ;
	vv = v3DistAzi( pCenter,100.0/ERadius,-3.14159/4.0) ;
	v32A(vv,b1) ; fprintf(outFile,"_%s_\n",b1) ;
	makeRotation(&a,pCenter,nPole,1.0) ;
	transposeMat33(&b,&a,1.0) ;
	v3Print(outFile,nPole,"nPole  ") ;
	v3Print(outFile,pCenter,"pCenter") ;
	vv = v3MMult(&a,pCenter) ;
	v3Print(outFile,vv,"vv     ") ;
	v2 = v3MMult(&b,vv) ;
	v3Print(outFile,v2,"v2     ") ;
	v2 = v3MMult(&b,nPole) ;
	v3Print(outFile,v2,"v2     ") ;
	vv = v3MMult(&b,v2) ;
	v3Print(outFile,vv,"vv     ") ;
	testProjection() ;
}
void test3() 
{
	int list[200],n,i ;
	double xa[3],xb[3] ;
	for( i = 0 ; i < 3 ; i++) list[i] = i ;
	n = locate3(0,1,2,xa,xb) ;
	printX("xa =",xa,3) ;
	if( n == 0 ) exit(0) ;
	locateNL(list,3,3,xa) ;
	locateNL(list,3,3,xa) ;
	if( n == 1 ) exit(0) ;
	printX("xb =",xb,3) ;
	locateNL(list,3,3,xb) ;
	locateNL(list,3,3,xb) ;
	exit(0) ;
}
int main( int ac, char **av )
{
	int cc ;
	extern char *optarg ;
	static int flags = 3, groundVelocity ;
	while( EOF != (cc = getopt(ac,av,"i:m:l:3vLg"))) {
		switch(cc) {
			case 'i' : readEvent(optarg) ; break ;
			case 'l' : outFile = fopen(optarg,"w") ; break ;
			case '3' : test3() ; break ;
			case 'm' : locMode = atoi(optarg) ; break ;
			case 'v' : flags = flags & -3 ; break ;
			case 'L' : flags = flags & -2 ; break ;
			case 'g' : groundVelocity = 1 ; break ;
			default  : ;
		}
	}
	if(outFile) fprintf(outFile,"flags=%d\n",flags) ;
	locate( flags ) ;
	getQuality() ;
	printReport(flags & 1) ;
/*	printData() ; */
	testt() ;
	return 0 ;
}
