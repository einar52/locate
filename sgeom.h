/* © Copyright  2008, 2021 Einar Kjartansson */
/*
tools for geometry on a sphere
*/

typedef struct { double lat, lon ; } LatLon ;
	/* lat and long should be in degrees */
typedef struct { double x,y,z ; } Vec3 ;
/* 3d vector, e.g. unit vector to indicate postion on sphere */
typedef enum { GC, SCC, SCA } ArcType ;
/*	GC	Great circle
	SCC	Small circle, clockwise
	SCA	Small circle, anticlockwise
*/
typedef struct { double a11,a21,a31,a12,a22,a32,a13,a23,a33 ; }  Mat33 ;
#define	DEG (180.0/M_PI) 
#define R2D (180.0/M_PI)
#define D2R (M_PI/180.0)

char *deg2SecN( double deg ,  char *buffer ) ;
char *deg2SecE( double deg ,  char *buffer ) ;
char *deg2AN( double deg ,  char *buffer ) ;
char *deg2AE( double deg ,  char *buffer ) ;
char *latLon2A( LatLon ll , char *buffer ) ;
/* convert LatLon to ascii according to ICAO conventions */

Vec3 v3New(double lat, double lon) ;
char *v32Ascii( Vec3 v, char *buffer ) ;
Vec3 latLon2Vec3( LatLon ll) ;
char *v32Ascii( Vec3 v, char *buffer ) ;
LatLon v32LatLon( Vec3 v) ;
char *v32A( Vec3 v, char *buffer) ;
void v3Print(FILE *fd, Vec3 a, char *text ) ;
double v3dotProd(Vec3 a, Vec3 b) ;
Vec3 v3OuterProd(Vec3 a, Vec3 b) ;
double v3Length( Vec3 a ) ;
int v3TestEQ( Vec3 a, Vec3 b) ;
/* test for equality, retun 1 if not equal */
Vec3 v3Normalize( Vec3 a ) ;
Vec3 v3Scale( double scale, Vec3 a ) ;
double v3Angle2( Vec3 a, Vec3 b) ;
double v3Turn(Vec3 a, Vec3 b, Vec3 c, double offset) ;
/* Return angle between ab and bc after adding offset
 * Result is in range between -pi and pi */
double v3Angle3A(Vec3 a, Vec3 b, Vec3 c) ;
double v3Angle3(Vec3 a, Vec3 b, Vec3 c) ;
	/* return the angle abc             */  
double v3Vertex(Vec3 a, Vec3 b, Vec3 c) ;
/* return angle beetween ba og bc , range of result -PI to PI */

double v3SumTurns(int n, Vec3 *node, Vec3 *center, ArcType *flag) ;

Vec3 v3Sum( Vec3 a, Vec3 b) ;
Vec3 v3Difference( Vec3 a, Vec3 b) ;
double v3Angle(Vec3 aA, Vec3 bB, double a, double b) ;
/*
    Return the cosine of angle B for a triangle where locationvectors 
    for A and B and lengths of sides a and b are given. 
    Robust for small triangles.
*/
int v3Intersects( Vec3 aA, double r1,
	Vec3 bB, double r2, Vec3 *cC1, Vec3 *cC2) ;
/* find possible third corners form triangle with corners at aA and bB,
	where lengths of sides are a and b
	0 is returned when there are no solutions.
*/
double v3Area( Vec3 a, Vec3 b, Vec3 c) ;
 /*      return area of trangle defined by ABC, the result is signed */

Vec3 v3DistAzi( Vec3 a, double dist, double azim ) ;
/* return vector a moved dist in direction azim */

Vec3 v3MMult( Mat33 *a , Vec3 b ) ;

void makeRotation( Mat33 *a, Vec3 b, Vec3 c, double scale ) ;
/*	make rotation operator that maps b to c, with scale. */

void transposeMat33( Mat33 *a, Mat33 *b, double scale) ;
	/* matrix transpose, with scale, would not work in place */

