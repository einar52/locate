/*---------------------------------------------------------------------------- 
 * © Copyright 1998, 2005, 2020 Einar Kjartansson
 *
 *
 *   MODULE_ID : v3_sgeom.c
 *   TITLE     : 
 *   AUTHOR    : Einar Kjartansson
 *   DATE      : Jan 02 2001
 *
 *   PURPOSE   : Spherical geometry tools, including code determine direction
 *   			of definition of airspace reservation.
 *   Task : 1208

 REVISION_HISTORY

 *--------------------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include "sgeom.h"
/*#include "asdm.h" */
#define SRAD 3600
	char b1[30],b2[30],b3[30],b4[30] ;

char *deg2SecN( double deg , char *buffer )
{
	int sec,s,m,d ;
	char ns ;
	if( deg > 0.0 ){
		sec = rint(deg*SRAD) ;
		ns = 'N' ;
	}else {	sec = rint(-deg*SRAD) ;
		ns = 'S' ;
	}
	d = sec/3600 ;
	m = (sec/60)%60 ;
	s = sec%60 ;
	sprintf(buffer,"%c%02d°%02d\'%02d\"",ns,d,m,s) ;
	return(buffer) ;
}
char *deg2SecE( double deg , char *buffer )
{
	int sec,s,m,d ;
	char ns ;
	if( deg > 0.0 ){
		sec = rint(deg*SRAD) ;
		ns = 'E' ;
	}else {	sec = rint(-deg*SRAD) ;
		ns = 'W' ;
	}
	d = sec/3600 ;
	m = (sec/60)%60 ;
	s = sec%60 ;
	sprintf(buffer,"%c%03d°%02d\'%02d\"",ns,d,m,s) ;
	return(buffer) ;
}

char *deg2AN( double deg , char *buffer )
{
	int sec,s,m,d ;
	char ns ;
	if( deg > 0.0 ){
		sec = rint(deg*SRAD) ;
		ns = 'N' ;
	}else {	sec = rint(-deg*SRAD) ;
		ns = 'S' ;
	}
	d = sec/3600 ;
	m = (sec/60)%60 ;
	s = sec%60 ;
	if( s) sprintf(buffer,"%02d%02d%02d%c",d,m,s,ns) ;
	else if (m)  sprintf(buffer,"%02d%02d%c",d,m,ns) ;
	else sprintf(buffer,"%02d%c",d,ns) ; 
	return(buffer) ;
}
char *deg2AE( double deg , char *buffer )
{
	int sec,s,m,d ;
	char ns ;
	if( deg > 0.0 ){
		sec = rint(deg*SRAD) ;
		ns = 'E' ;
	}else {	sec = rint(-deg*SRAD) ;
		ns = 'W' ;
	}
	d = sec/3600 ;
	m = (sec/60)%60 ;
	s = sec%60 ;
	if( s) sprintf(buffer,"%03d%02d%02d%c",d,m,s,ns) ;
	else if (m)  sprintf(buffer,"%03d%02d%c",d,m,ns) ;
	else sprintf(buffer,"%03d%c",d,ns) ; 
	return(buffer) ;
}
char *latLon2A( LatLon ll , char *buffer ) 
{
	char *cp,b1[30],b2[30] ;
	cp = deg2AN(ll.lat,b1) ;
	cp = deg2AE(ll.lon,b2) ;
	sprintf(buffer,"%s%s",b1,b2) ;
	return( buffer ) ;
}
Vec3 v3New(double lat, double lon)
{
        double sb,sl,cb,cl ;
	Vec3 v ;
        cb = cos(lat*M_PI/180.0) ;
        cl = cos(lon*M_PI/180.0) ;
        sb = sin(lat*M_PI/180.0) ;
        sl = sin(lon*M_PI/180.0) ;
        v.z = sb ;
        v.x = cb * cl ;
        v.y = cb * sl ;
	return(v) ;
}
Vec3 latLon2Vec3( LatLon ll)
{
	return(v3New(ll.lat,ll.lon)) ;
}
char *v32Ascii( Vec3 v, char *buffer )
{
	sprintf(buffer,"(%9.6f,%9.6f,%9.6f)",v.x,v.y,v.z) ;
	return(buffer) ;
}
LatLon v32LatLon( Vec3 v)
{
	double rad ;
	LatLon ll ;
	rad = 	sqrt(v.x*v.x + v.y*v.y) ;
	ll.lat = atan2(v.z,rad)*DEG ;
	ll.lon = atan2(v.y,v.x)*DEG ;
	return(ll) ;
}
char *v32A( Vec3 v, char *buffer) 
{
	return(latLon2A(v32LatLon(v),buffer))	 ;
}
double v3DotProd(Vec3 a, Vec3 b) 
{
	return(a.x*b.x + a.y*b.y + a.z*b.z) ;
}
Vec3 v3OuterProd(Vec3 a, Vec3 b)
{
	Vec3 v ;
	v.x = a.y*b.z  - a.z*b.y ;
	v.y = a.z*b.x  - a.x*b.z ;
	v.z = a.x*b.y  - a.y*b.x ;
	return(v) ;
}
double v3Length( Vec3 a )
{
	return( sqrt( a.x*a.x + a.y*a.y + a.z*a.z ) ) ;
}
int v3TestEQ( Vec3 a, Vec3 b)
/* test for equality, retun 1 if not equal */
{
	if( a.x != b.x ) return 0 ;
	if( a.y != b.y ) return 0 ;
	if( a.z != b.z ) return 0 ;
	return 1 ;
}
Vec3 v3Normalize( Vec3 a ) 
{
	Vec3 c ;
	double r ;
	r = sqrt(a.x*a.x + a.y*a.y + a.z*a.z) ;
	if( r == 0 ) return(a) ;
	c.x = a.x/r ; c.y = a.y/r ; c.z = a.z/r ;
	return(c) ;
}
double v3Angle2( Vec3 a, Vec3 b)
{
/* This computes angle between two vectors, the result is between 0 and PI.
   This is the shortest distance between two points on a sphere */
	double si,co ;
	co = v3DotProd(a,b) ;
	si = v3Length( v3OuterProd(a,b)) ;	
	return( atan2(si,co)) ;
}
Vec3 v3Scale( double scale, Vec3 a )
{
	Vec3 c ;
	c.x = scale*a.x ; c.y = scale*a.y ; c.z = scale*a.z ;
	return(c) ;
}
double v3Angle3A(Vec3 a, Vec3 b, Vec3 c)
/* return the angle abc             */
{
	Vec3 v1,v2 ;
	v1 = v3Normalize(v3OuterProd(b,a)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	return(v3Angle2(v1,v2) ) ;
}
double v3Turn(Vec3 a, Vec3 b, Vec3 c, double offset)
/* Return angle between ab and bc after adding offset
 * Result is in range between -pi and pi */
{
	Vec3 v1,v2,v3 ;
	double si,co,res ;
	v1 = v3Normalize(v3OuterProd(a,b)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	v3 = v3OuterProd(v1,v2) ;
	si = v3DotProd(b,v3) ;
	co = v3DotProd(v1,v2) ;
	res = atan2(si,co)  ;
	res += offset ;
	if( res < -M_PI) res += 2.0*M_PI ;
	if( res > M_PI) res -= 2.0*M_PI ;
	return(res) ;
}
double v3ArcTurn(Vec3 a, Vec3 b, Vec3 c, int mode)
/* Return angle between ba and bc, adding pi if mode
 * Result is in range between -pi and pi */
{
	Vec3 v1,v2,v3 ;
		char b1[30],b2[30],b3[30],b4[30] ;
	double si,co,res ;
	v1 = v3Normalize(v3OuterProd(b,a)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	v3 = v3OuterProd(v1,v2) ;
	si = v3DotProd(b,v3) ;
	co = v3DotProd(v1,v2) ;
	res = atan2(si,co)  ;
	if ( mode ) return(res) ;
	res += M_PI ;
	if( res > M_PI ) return(res - 2.0*M_PI ) ;
	return(res) ;
}

double v3Angle3(Vec3 a, Vec3 b, Vec3 c)
/* return the angle between ab and bc            */
/* finds angle between planes that include a,b and b,c.
   The results is signed ( range -PI to +PI )
*/
{
	Vec3 v1,v2,v3 ;
	double si,co,res ;
	v1 = v3Normalize(v3OuterProd(a,b)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	v3 = v3OuterProd(v1,v2) ;
	si = v3DotProd(b,v3) ;
	co = v3DotProd(v1,v2) ;
	res = atan2(si,co) ;
/*
	printf("%s  %s  %s  ",v32A(a,b1),v32A(b,b2),v32A(c,b3)) ;
	printf("res=%10.6f\n",res*DEG) ;
*/
	return(res) ;
}
double v3Vertex(Vec3 a, Vec3 b, Vec3 c)
/* return angle beetween ba og bc , range of result -PI to PI */
{
	Vec3 v1,v2,v3 ;
	double si,co,res ;
	v1 = v3Normalize(v3OuterProd(b,a)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	v3 = v3OuterProd(v1,v2) ;
	si = v3DotProd(b,v3) ;
	co = v3DotProd(v1,v2) ;
	return atan2(si,co) ;
}
double v3AngleV(Vec3 a, Vec3 b, Vec3 c, int antiClock)
{
	double an ;
	Vec3 v1,v2,v3 ;
	double si,co,res ;
	v1 = v3Normalize(v3OuterProd(b,a)) ;
	v2 = v3Normalize(v3OuterProd(b,c)) ;
	v3 = v3OuterProd(v1,v2) ;
	si = v3DotProd(b,v3) ;
	co = v3DotProd(v1,v2) ;
	an = atan2(si,co) ;
	if( antiClock && (an < 0.0 )) an += 2.0*M_PI ;
	if( (!antiClock) && (an > 0.0 )) an -= 2.0*M_PI ;
	return(an) ;
}
double v3AdjustAngle( double u )
{	/* Make sure    -pi < u <= pi       */
	while (u > M_PI) u -= 2.0*M_PI ;
	while (u <= -M_PI) u += 2.0*M_PI ;
	return(u) ;
}
char *v3DistA( Vec3 c, Vec3 n, char *buffer) 
{
	deg2AN(v3Angle2(c,n)*180.0/M_PI,buffer) ;
	buffer[strlen(buffer)-1] = 0 ;
	return(buffer) ;
}
double v3SumTurns(int n, Vec3 *node, Vec3 *center, ArcType *flag)
{
	int i,i0,i2 ;
	char line[200] ;
	double sum,an,an1,an2,turn,bias,deg ;
	sum = 0.0 ;
	for(i = 0 ; i < n ; i++) {
		i0 = (i+n-1) % n ;
		i2 = (i  +1) % n ;
		an2 = 0 ;
		if( flag[i0] != GC ) {
			if(flag[i] != GC) {
				an1 = v3ArcTurn(center[i0],node[i],center[i],
					flag[i]==flag[i0]) ;
				an2 = v3AngleV(node[i],center[i],node[i2],flag[i]==SCA) ;
			}else {
				if(flag[i0] == SCA) bias = -0.5*M_PI ;
				else bias = +0.5*M_PI ;
				an1 = v3Turn(center[i0],node[i],node[i2],bias) ;
			}
		} else {
			if(flag[i] != GC ) {
				if(flag[i]== SCA) bias = -0.5*M_PI ;
				else bias = 0.5*M_PI ;
				an1 = v3Turn(node[i0],node[i],center[i],bias) ;
				an2 = v3AngleV(node[i],center[i],node[i2],
						(flag[i]== SCA)) ;
			}else {
				an1 = v3Turn(node[i0],node[i],node[i2],0.0) ;
			}
		}
		an = an1 + an2 ;
		sum += an ;	
		if(flag[i]) {
			sprintf(line,"%-17s%-17s%2d%9.3f%9.3f%9.3f%9.3f %-6s %-6s",
				v32A(node[i],b1),
				v32A(center[i],b2),flag[i],
				an1*DEG,an2*DEG,an*DEG,sum*DEG,
				v3DistA(center[i],node[i],b3),
				v3DistA(center[i],node[i2],b4)) ;
		}else {
			sprintf(line,"%-17s%-17s%2d%9.3f%18s%9.3f",
				v32A(node[i],b1)," ",flag[i], 
				an1*DEG," ",sum*DEG) ;
		}
/*		CSDEBUGLV( AIRPLV,line );
		fprintf(stderr,"%s\n",line) ; */
	}
	return(sum) ;
}
Vec3 v3Sum( Vec3 a, Vec3 b)
{
	Vec3 c ;
	c.x = a.x + b.x ;
	c.y = a.y + b.y ;
	c.z = a.z + b.z ;
	return(c) ;
}
Vec3 v3Difference( Vec3 a, Vec3 b)
{
	Vec3 c ;
	c.x = a.x - b.x ;
	c.y = a.y - b.y ;
	c.z = a.z - b.z ;
	return(c) ;
}
double v3Angle(Vec3 aA, Vec3 bB, double a, double b)
/*
    Return the cosine of angle B for a triangle where locationvectors 
    for A and B and lengths of sides a and b are given. 
    Robust for small triangles.
*/
{
	double sum,diff,rr,sc2,sa2,sb2,cB;
	double sa22,sc22,sa,sc ;
	sum  = v3Length(v3Sum(aA,bB)) ;
	diff = v3Length(v3Difference(aA,bB)) ;
	rr = sqrt(sum*sum + diff*diff) ;
	sa2  = sin(a*0.5) ;
	sb2  = sin(b*0.5) ;
	sc2  = diff/rr ;
	sa22 = sa2*sa2 ;
	sc22 = sc2*sc2 ;
	sa = sin(a) ;
	sc = 2.0*sc2*sum/rr ;
	cB = (sa22+sc22-sb2*sb2-2.0*sa22*sc22)/(0.5*sa*sc) ;
	return(cB) ;
}
int v3Intersects( Vec3 aA, double r1,
	Vec3 bB, double r2, Vec3 *cC1, Vec3 *cC2)
{       
	Vec3 n1,n2,n3 ;
	double cosa,sina,sid,cod ;
	if( r1+r2 <= v3Angle2(aA,bB)) return(0) ;
	cosa = v3Angle(aA,bB,r1,r2) ;
	sina = sqrt(1.0 - cosa*cosa) ;
	sid = sin(r1) ;
	cod = cos(r1) ;
	n1 = v3Normalize(v3OuterProd(aA,bB)) ;
	n2 = v3OuterProd(n1,aA) ;
	n3 = v3Sum( v3Scale(cosa,n2) , v3Scale(-sina,n1) ) ;
	*cC1 = v3Sum( v3Scale(cod,aA)   , v3Scale(sid,n3) ) ;
	n3 = v3Sum( v3Scale(cosa,n2) , v3Scale(sina,n1) ) ;
	*cC2 = v3Sum( v3Scale(cod,aA)   , v3Scale(sid,n3) ) ;
	return(1) ;
}
double v3Area( Vec3 a, Vec3 b, Vec3 c)
{
/*	return area of trangle defined by ABC, the result is signed */
	double sum,s1,s2,s3 ;
	if ( v3TestEQ( a,b) ) return 0.0 ;
	if ( v3TestEQ( b,c) ) return 0.0 ;
	if ( v3TestEQ( c,a) ) return 0.0 ;
	s1 = v3Vertex(a,b,c); 
	s2 = v3Vertex(b,c,a); 
	s3 = v3Vertex(c,a,b);
	sum = s1+s2+s3 ;
/*	printf("%9.4f %9.4f %9.4f %9.4f\n",s1,s2,s3,sum) ; */
	if( sum > 0.0 ) sum -= M_PI ;
	if( sum < 0.0 ) sum += M_PI ;
	return sum ;
}
Vec3 v3DistAzi( Vec3 a, double dist, double azim )
/* return vector a moved dist in direction azim */
{
	Vec3 n1,n2,n3,north ;
	north.z = 1.0 ; north.x = 0.0 ; north.y = 0.0 ;
	n1 = v3Normalize(v3OuterProd(north,a)) ;
	n2 = v3OuterProd(a,n1) ;
	n3 = v3Sum( v3Scale(cos(azim),n2),v3Scale(sin(azim),n1)) ;
	return v3Sum( v3Scale(cos(dist),a),v3Scale(sin(dist),n3)) ;
}
Vec3 v3MMult( Mat33 *a , Vec3 b )
{
	Vec3 r ;
	r.x = a->a11*b.x + a->a12*b.y + a->a13*b.z ;
	r.y = a->a21*b.x + a->a22*b.y + a->a23*b.z ;
	r.z = a->a31*b.x + a->a32*b.y + a->a33*b.z ;
	return r ;
}
void transposeMat33( Mat33 *a, Mat33 *b, double scale)
	/* matrix transpose, with scale, would not work in place */
{	
	a->a11 = scale * b->a11 ;
	a->a12 = scale * b->a21 ;
	a->a13 = scale * b->a31 ;
	a->a21 = scale * b->a12 ;
	a->a22 = scale * b->a22 ;
	a->a23 = scale * b->a32 ;
	a->a31 = scale * b->a13 ;
	a->a32 = scale * b->a23 ;
	a->a33 = scale * b->a33 ;
}
void makeRotation( Mat33 *a, Vec3 b, Vec3 c, double scale )
{
/*	make rotation operator that maps b to c, with scale. */
	Vec3 n1,n2,n3 ;
	n1 = v3Normalize(v3OuterProd(c,b)) ;
	n2 = v3OuterProd(b,n1) ;
/*	v3Print(n1,"n1  ") ;
	v3Print(n2,"n2  ") ;
	v3Print(b,"b   ") ;
*/
	a->a11 = scale*n1.x ;
	a->a12 = scale*n1.y ;
	a->a13 = scale*n1.z ;
	a->a21 = scale*n2.x ;
	a->a22 = scale*n2.y ;
	a->a23 = scale*n2.z ;
	a->a31 = scale*b.x ;
	a->a32 = scale*b.y ;
	a->a33 = scale*b.z ;
}
void v3Print(FILE *fd, Vec3 a, char *text )
{
	LatLon ll ;
	ll = v32LatLon(a) ;
	fprintf(fd,"%s  %15.8f %15.8f (%15.11f %15.11f %15.11f)\n",
		text,ll.lat,ll.lon, a.x,a.y,a.z) ;
}
