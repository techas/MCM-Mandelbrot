/* gcc -fopenmp -lm mandel.c -o mandel */

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void color(int it, int maxit,FILE* fh);
void frame(int itermax, double x0, double y0, double magnify, int hyres, int hxres);

main()
{
	double x0 = 0.013438870532012129028364919004;
	double y0 = 0.655614218769465062251320027664;
	int itermax = 1000000;
	double magnify[200] = {
   0.000100000000000,   0.000106376485432,   0.000113159566528,
   0.000120375169802,   0.000128050874968,  0.000136216020355,
   0.000144901815049,   0.000154141458175,   0.000163970265800,
   0.000174425805911,   0.000185548042014,   0.000197379485882,
   0.000209965360044,   0.000223353770639,   0.000237595891284,
   0.000252746158678,   0.000268862480665,   0.000286006457576,
   0.000304243617677,   0.000323643667635,   0.000344280758952,
   0.000366233771390,   0.000389586614469,   0.000414428548184,
   0.000440854524183,   0.000468965548693,   0.000498869068584,
   0.000530679382066,   0.000564518075552,   0.000600514488398,
   0.000638806207266,   0.000679539592008,   0.000722870335095,
   0.000768964056702,   0.000817996937752,   0.000870156393319,
   0.000925641788971,   0.000984665202794,   0.001047452236001,
   0.001114242875233,   0.001185292409845,   0.001260872407681,
   0.001341271753068,   0.001426797751001,   0.001517777301732,
   0.001614558150262,   0.001717510215498,   0.001827027004177,
   0.001943527114930,   0.002067455838273,   0.002199286858606,
   0.002339524064745,   0.002488703475903,   0.002647395290480,
   0.002816206065495,   0.002995781034986,   0.003186806576246,
   0.003390012833314,   0.003606176507761,   0.003836123827418,
   0.004080733704412,   0.004340941094578,   0.004617740571069,
   0.004912190125854,   0.005225415213603,   0.005558613053441,
   0.005913057204994,   0.006290102436234,   0.006691189901716,
   0.007117852651001,   0.007571721488337,   0.008054531205965,
   0.008568127214899,   0.009114472598521,   0.009695655615936,
   0.010313897683787,   0.010971561867027,   0.011671161911101,
   0.012415371850064,   0.013207036227366,   0.014049180968355,
   0.014945024946065,   0.015897992284505,   0.016911725446448,
   0.017990099155778,   0.019137235207583,   0.020357518222611,
   0.021655612406318,   0.023036479376537,   0.024505397127943,
   0.026067980205769,   0.027730201165911,   0.029498413403417,
   0.031379375436648,   0.033380276739903,   0.035508765223261,
   0.037772976464675,   0.040181564806038,   0.042743736432097,
   0.045469284558624,   0.048368626864372,   0.051452845309859,
   0.054733728495201,   0.058223816718887,   0.061936449909706,
   0.065885818615068,   0.070087018240569,   0.074556106748144,
   0.079310166033331,   0.084367367216249,   0.089747040095843,
   0.095469747032875,   0.101557361544042,   0.108033151907647,
   0.114921870100370,   0.122249846405078,   0.130045090051290,
   0.138337396272962,   0.147158460192806,   0.156541997968413,
   0.166523875663253,   0.177142246335109,   0.188437695865931,
   0.200453398090524,   0.213235279816976,   0.226832196369604,
   0.241296118325361,   0.256682330157469,   0.273049641545533,
   0.290460612159805,   0.308981790778801,   0.328683969654205,
   0.349642455095314,   0.371937355307265,   0.395653886583227,
   0.420882699020889,   0.447720223008213,   0.476269037802799,
   0.506638263613665,   0.538943978684060,   0.573309662969503,
   0.609866670106889,   0.648754729478629,   0.690122480290852,
   0.734128039707011,   0.780939607208450,   0.830736107491937,
   0.883707874361460,   0.940057378228296,   1.000000000000000};

	int hxres = 500;
	int hyres = 500;
	int id, i;

	for(i=0; i<150; i++){
		magnify[i] = magnify[i]*1e4;
	}


#pragma omp parallel
	{
#pragma omp for
		for(i=0; i<150; i++) {
			id = omp_get_thread_num();
			fprintf(stdout,"core %i computing mag=%f\n",id,magnify[i]);	
			frame(itermax,x0,y0,magnify[i],hyres,hxres);
		}
	}
/* #pragma omp parallel num_threads(2)
	{
		id = omp_get_thread_num();
		ix = id + 8;
		frame(itermax,x0,y0,magnify[ix],hyres,hxres);
	}	
*/
}

void frame(int itermax, double x0, double y0, double magnify, int hyres, int hxres) 
{
	double x,xx,y,cx,cy;
	int iteration,hx,hy,ti;
	int tim = (int) log((double)itermax);
	FILE *fh;
	char fileName[100];
	sprintf(fileName,"frame-mag-%f.ppm",magnify);
	fh = fopen(fileName,"w");

	/* header for PPM output */
	fprintf(fh,"P6\n# mandel program\n");
	fprintf(fh,"%d %d\n255\n",hxres,hyres);

	for (hy=1;hy<=hyres;hy++)  {
		for (hx=1;hx<=hxres;hx++)  {
			cx = (((float)hx)/((float)hxres)-0.5)/magnify*3.0-x0;
			cy = (((float)hy)/((float)hyres)-0.5)/magnify*3.0-y0;
			x = 0.0; y = 0.0;
			for (iteration=1;iteration<itermax;iteration++)  {
				xx = x*x-y*y+cx;
				y = 2.0*x*y+cy;
				x = xx;
				if (x*x+y*y>100.0)
					break;
			}
			ti = (int)log((double)iteration);
			color(ti,tim,fh);
		}
	}
	fclose(fh);
}

void color(int it, int maxit,FILE* fh)  {
   int ix = (63*it)/maxit;
   int r,g,b;
   switch(ix) {
      case 0:
         r=128; g=0; b=0;
      break;
      case 1:
         r=139; g=0; b=0;
      break;
      case 2:
         r=151; g=0; b=0;
      break;
      case 3:
         r=162; g=0; b=0;
      break;
      case 4:
         r=174; g=0; b=0;
      break;
      case 5:
         r=185; g=0; b=0;
      break;
      case 6:
         r=197; g=0; b=0;
      break;
      case 7:
         r=209; g=0; b=0;
      break;
      case 8:
         r=220; g=0; b=0;
      break;
      case 9:
         r=232; g=0; b=0;
      break;
      case 10:
         r=241; g=2; b=0;
      break;
      case 11:
         r=249; g=7; b=0;
      break;
      case 12:
         r=253; g=15; b=0;
      break;
      case 13:
         r=255; g=28; b=0;
      break;
      case 14:
         r=255; g=43; b=0;
      break;
      case 15:
         r=255; g=58; b=0;
      break;
      case 16:
         r=255; g=74; b=0;
      break;
      case 17:
         r=255; g=90; b=0;
      break;
      case 18:
         r=255; g=105; b=0;
      break;
      case 19:
         r=255; g=119; b=0;
      break;
      case 20:
         r=255; g=131; b=0;
      break;
      case 21:
         r=255; g=142; b=0;
      break;
      case 22:
         r=255; g=154; b=0;
      break;
      case 23:
         r=255; g=166; b=0;
      break;
      case 24:
         r=255; g=179; b=0;
      break;
      case 25:
         r=255; g=192; b=0;
      break;
      case 26:
         r=255; g=206; b=0;
      break;
      case 27:
         r=255; g=219; b=0;
      break;
      case 28:
         r=254; g=231; b=1;
      break;
      case 29:
         r=252; g=240; b=3;
      break;
      case 30:
         r=248; g=247; b=7;
      break;
      case 31:
         r=243; g=251; b=12;
      break;
      case 32:
         r=236; g=253; b=19;
      break;
      case 33:
         r=228; g=254; b=27;
      break;
      case 34:
         r=220; g=255; b=35;
      break;
      case 35:
         r=210; g=255; b=45;
      break;
      case 36:
         r=200; g=255; b=55;
      break;
      case 37:
         r=189; g=255; b=66;
      break;
      case 38:
         r=177; g=255; b=78;
      break;
      case 39:
         r=165; g=255; b=90;
      break;
      case 40:
         r=153; g=255; b=102;
      break;
      case 41:
         r=142; g=255; b=113;
      break;
      case 42:
         r=131; g=255; b=124;
      break;
      case 43:
         r=121; g=255; b=134;
      break;
      case 44:
         r=111; g=255; b=144;
      break;
      case 45:
         r=101; g=255; b=154;
      break;
      case 46:
         r=91; g=255; b=164;
      break;
      case 47:
         r=80; g=255; b=175;
      break;
      case 48:
         r=70; g=255; b=185;
      break;
      case 49:
         r=60; g=255; b=195;
      break;
      case 50:
         r=49; g=254; b=206;
      break;
      case 51:
         r=38; g=249; b=217;
      break;
      case 52:
         r=26; g=235; b=229;
      break;
      case 53:
         r=13; g=199; b=242;
      break;
      case 54:
         r=5; g=157; b=250;
      break;
      case 55:
         r=1; g=120; b=254;
      break;
      case 56:
         r=0; g=91; b=254;
      break;
      case 57:
         r=0; g=66; b=252;
      break;
      case 58:
         r=0; g=43; b=247;
      break;
      case 59:
         r=0; g=25; b=236;
      break;
      case 60:
         r=0; g=12; b=220;
      break;
      case 61:
         r=0; g=5; b=198;
      break;
      case 62:
         r=0; g=1; b=172;
      break;
      case 63:
         r=0; g=0; b=143;
      break;
   }
   fputc((char)r,fh);
   fputc((char)g,fh);
   fputc((char)b,fh);
}


/*void color(int red, int green, int blue)  {
	fputc((char)red,stdout);
	fputc((char)green,stdout);
	fputc((char)blue,stdout);
}*/

