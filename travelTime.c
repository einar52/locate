/*  © Copyright  2008, 2020 Einar Kjartansson  */

double travelTime( double dist )
{
	/* fitted to model Hengill - Hornafjordur, depth of 4 km */
	static double a[] = 
	 {4.78090,130.43203,-183.26574, 450.35999,-400.24155,128.81572} ;
	double sum,cen,xx ;
	int j ;
	sum = a[0] + a[1]*0.001*dist ;
	for( j = 2 ; j < 6 ; j++) {
		cen = 20.0*j - 32 ;
		xx = (dist-cen)*0.005 ;
		sum = sum + a[j]/(1.0 + xx*xx) ;
	}
	return sum ;
}
