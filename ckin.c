#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <ctype.h>
#include <sys/stat.h>
#include <termios.h>

/*define the path to the mass-table file here before compiling*/
/*#define M_TAB   "./mass_table_2003.txt"*/
#define M_TAB 	"./mass_table_2020.txt"
#define M_TABYR 2020    /*mass table year*/
#define R_MAX 	10      /*max no. of isotopes in reaction*/
#define R0 	1.4	/*R_0 in fm*/
#define DMAX 	8	/*max number of decay particles stored*/
#define SMAX 	5	/*max number of special notation isotopes*/
#define MXEX    10      /*max. number of excitation and ang. mom. values*/
#define CHLEN   200     /*character length of char arrays*/

#ifdef HAVE_WCLBES
#define NUMOPT  7       /*number of options*/
#else
#define NUMOPT  6       /*number of options*/
#endif
/* Carl Wheldon April 2008 with physics input from Tzany Wheldon*/

/* To compile:
gcc ckin.c -Wall -pedantic -o ckin -lm -O2
    or, if you have libwclbes.a
gcc ckin.c -Wall -o ckin -lm -O2 -DHAVE_WCLBES ./libwclbes.a -lgfortran
    or, if you have the cern libraries:
gcc ckin.c -Wall -o ckin -lm -O2 -DHAVE_WCLBES -L./cernlib -lmathlib -lkernlib -lgfortran
*/

/*%%%% A program to work out some (relativistic) reaction kinematics %%%%*/

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* +++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int     ask_yn(char *questn, int yn);
int 	extract_react(char react[], int opt);
int 	kin(char react[], int opt);
int 	chk_react(char react[]);
void    coul_pen(int num);
void 	decay_part(void);
void 	format(char react[], int *fmt);
void 	get_ans(char ans[], int num, int md);
void    get_line(char ans[], int len);
int 	get_mode(int md);
int 	get_pars(char ans0[], double pars[], int num);
int 	qval(char react[], int md);
int	read_masses(int i, int md);
void 	skip_hash(FILE *file);
void	special();
void    store_colours();
void 	zero_par(void);

#ifdef HAVE_WCLBES
/*from library libwclbes.a*/
void    wclbes_(double complex *kr, double complex *eta, double complex *zlmin,
        int *nl, double complex f[], double complex g[], double complex fp[],
        double complex gp[], double complex *sig, int *kfn, int *mode,
        int *jfail, int *jpr);
#endif

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* +++++++++++++++++++++++++++ EXTERNAL VARIABLES +++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
double  amu, c, elec, eps0, fm, hbar, mev, pi;
int 	z[R_MAX], a[R_MAX], decz[DMAX], deca[DMAX], dec[DMAX], dk;
int 	sa[SMAX], sz[SMAX];
float	m[R_MAX], dm[R_MAX], ba[R_MAX], dba[R_MAX], q, dq, decm[DMAX];
float	bea[R_MAX], dbea[R_MAX];
char    el[R_MAX][2], decel[DMAX][2], sel[SMAX][2], ssel[SMAX][2], clr[10][12];
FILE	*fmtab;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int main(void)
{
    extern  FILE    *fmtab;
    int     flg = -1, md = 0;
    char    react[R_MAX*5] = "";        
    
    printf("\n\t\t**** Welcome to program ckin ****\n   This works out"
    	   " Q-values and some 2-body reaction kinematics.\n"
	    "\t  (Both relativistic and non-relativistic.)\n");

    /*store colours*/
    store_colours();
        
    /*open mass table file for reading*/
    if ((fmtab= fopen(M_TAB, "r" )) == NULL)
    {
    	printf("Cannot open mass-table file: %s \n", M_TAB);
    	return -1;			
    }
    else printf("Opened mass-table file: %s\n",M_TAB);
    
    /*skip inital lines starting with # in mass-table file*/
    skip_hash(fmtab);
          
    /*store decay particles*/
    decay_part();
    /*store special notation isotopes*/
    special();
    
    /*define constants*/
    amu = (double)1.66053886e-27;
    c = (double)299792458;
    elec = (double)1.60217653e-19;
    eps0 = (double)8.85418782e-12;
    fm = (double)1.0e-15;
    hbar = (double)1.05457148e-34;
    mev = (double)1.60217653e-13;
    pi = (double)3.14159265359;
    
    while (1)
    {
	if ((md = get_mode(md)) == 1 || md == 5)
	{
	    memset(react, ' ', R_MAX*5);
	    memset(react, '\0', 1);
    	    if (md == 1) printf("Enter reaction [e.g. c12(c12,o16,he4) or "
                "p(6he,n)] or decay [e.g. (55fe,e+)]\n");
            else if (md == 5) printf("Enter isotope(s) (e.g. 187Re,d,76ge)"
                                    " [max. %d]\n",R_MAX);
            
    	    get_line(react, R_MAX*5);
/*    	    strcpy(react,"12c(197au,197au)");*/
    	    printf("Reaction/decay/isotope: %s\n",react);
	    if (md != 5) flg = qval(react, 0);
	    else extract_react(react, md);
	}
    	else if (md == 2)
	{
	    if ( (flg = qval(react, 0)) < 0 )
		 printf("***You must enter a valid reaction!\n");
	}
    	else if (md == 3 || md == 4)
	{
	    if (dk != 0) printf("***You must enter a reaction not decay!\n");
	    else if (flg != 4)
		printf("***You must enter a reaction with 2 incoming and 2"
		" outgoing particles!\n");
	    else kin(react, md);
	}
    	else if (md == 6)
	{
	    if (dk != 0) printf("***You must enter a reaction not decay!\n");
	    else if (flg != 3)
		printf("***You must enter a RESONANT reaction with 2 incoming and 1"
		" outgoing particles!\n");
	    else kin(react, md);
	}
#ifdef HAVE_WCLBES
        else if (md == NUMOPT)
        {
            if (dk != 0) coul_pen(flg);
            else printf("***You must enter a valid decay!\n");
        }
#endif
 	else if (md == 0) break;
    }
     	    
    printf("\n\tGoodbye....\n\n");
    
    /*close mass-table file*/
    fclose(fmtab);
    return 0;
} /*END main()*/

/*==========================================================================*/
/* coul_pen: calculated penetrabilites and reduced widths for decays        */
/****************************************************************************/
#ifdef HAVE_WCLBES
void coul_pen(int num)
{
    /*decay threshold in MeV*/
    int lp, ifail = 0, ipr = 0, kfn, l[MXEX], mode, nl, nmex = 0, ernge = 0;
    int i, j, vb = 0, cntex = 0;
    double thresh = -q/(double)(1000.0), rwdth, wig, tmp[(MXEX*2)+2];
    double brch[MXEX], ecm, ex[MXEX], gam[MXEX], k[MXEX], mu, rf;
    double fd[22], gd[22], p[MXEX][22], rex = 0.0;
    double complex eta, f[22], g[22], fp[22], gp[22], kr, sig[22], zlmin;
    float tmpf = 0.0;
    char ans[CHLEN] = "" ;

/*wclbes function prototype
void wclbes_(double complex kr, double complex eta, double complex zlmin,
        int *nl, double complex f[], double complex g[], double complex fp[],
        double complex gp[], double complex *sig, int *kfn, int *mode,
        int *jfail, int *jpr);
*/
    
    /*zero arrays*/
    for (i = 0 ; i < MXEX; i++)
    {
        brch[i] = 0.0;
        gam[i] = 0.0;
        ex[i] = 0.0;
        l[i] = 0;
    }
    
    printf("Enter excitation energy (MeV) in the recoil decay product, %d%c%c.",
            a[2],el[2][0],isalpha(el[2][1]) ? el[2][1] : (char)0);
    printf(" [<Enter> for %.2f]\n",rex);
    get_line(ans, CHLEN);
    if (strlen(ans) == (int)0 && (ans[0] == ((char)0)));
    else tmpf = atof(ans);
    rex = (double)tmpf;
    thresh += rex;
    printf("Decay threshold: %7.3f MeV\n",thresh);
    
    ernge = 0;
    /*lp == 0: 1st calculation with Ex & l
      lp == 1: 2nd calculation with expt. width & branching ratios
                  for reduced widths*/
    for (lp = 0; lp < 2; lp++)
    {
        /*zero tmp array*/
        for (i = 0 ; i < (MXEX*2)+2; i++) tmp[i] = 0.0;
        while (1)
        {
            if (!lp) printf("Enter values of Ex(%d%c%c) (MeV) and l (max. %d pairs of values)\n"
               " (for all l values, 0 --> 8, enter -1 for l)\n"
               " (for a range of Ex from El-->Eh in Esteps, enter El Eh Estep  l -2 -2)\n"
               "                                               or El Eh Estep -1 -2 -2)\n",
                a[0],el[0][0],isalpha(el[0][1]) ? el[0][1] : (char)0,MXEX);
            else printf("Enter %d value(s) of Gamma(tot.) (keV) and branching ratio"
                " Gamma(%d%c%c)/Gamma(tot.)\n [for recoil(%d%c%c) Ex=%.3f MeV]\n",
                nmex,a[1],el[1][0],isalpha(el[1][1]) ? el[1][1] : (char)0,
                a[2],el[2][0],isalpha(el[2][1]) ? el[2][1] : (char)0,rex);
            
            get_line(ans, CHLEN);
            printf("0nmex = %d\n",nmex);
            if ( ( ((i = get_pars(ans, tmp, (2*MXEX))) % 2 != 0 || i == 0 ||
                     i > 2*MXEX) && !lp) || (lp == 1 && i != (2*nmex)) )
                printf("***You must enter values in pairs up to a max. of %d"
                    " pairs of values***\n",MXEX);
            else
            {
                /*set nmex equal to half the number of parameters entered*/
                if (!lp) nmex = i/2;
                break;
            }
            printf("1nmex = %d\n",nmex);
        }
    
        /*extract values from tmp to appropriate arrays*/
        for (i = 0; i < nmex; i++)
        {
            printf("2nmex = %d\n",nmex);
            /*all instances not involving a range of ex*/
            if (!lp && ( (nmex < 3 || nmex > 3) || (nmex == 3 && tmp[2*2] != -2 && tmp[(2*2)+1] != -2) ))
            {
                /*ZERO flag for Elow, Ehigh and Estep option*/
                ernge = 0;
                ex[i] = tmp[i*2];
                l[i] = (int)tmp[(i*2)+1];
            }
            /*process ex range option*/
            else if (!lp)
            {
                if (vb) printf("### %f,%f,%f,%f,%f,%f\n",tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]);
                /*set flag for Elow, Ehigh and Estep option*/
                ernge = 1;
                /*extract values*/
                ex[0] = tmp[0]; /*Elow*/
                ex[1] = tmp[1]; /*Ehigh*/
                ex[2] = tmp[2]; /*Estep*/
                l[0] = (int)tmp[3];
                nmex = (int)(((ex[1]-ex[0])/ex[2])+0.5);
                /*setting i to nmex terminates for loop*/
                i = nmex;
                if (vb) printf("nmex= %d, ex[0]=%f,ex[1]=%f,ex[2]=%f,l[0]=%d\n",
                        nmex,ex[0],ex[1],ex[2],l[0]);
            }
            else
            {
                gam[i] = tmp[i*2];
                brch[i] = tmp[(i*2)+1];
            }
        }
        printf("Entered:\n");
        
        if (!ernge)
        {
            for (i = 0; i < nmex; i++)
            {
                if (!lp) printf("%d) ex:  %.3f, l: %d\n",i,ex[i],l[i]);
                else printf("%d) Gam:  %.3f, \tBrch: %.3f\n",i,gam[i],brch[i]);
            }
        }
        if (vb) printf("0 got to here\n");
        /*ernge is flg for ex range*/
        if (ernge) printf("Erange: %f --> %f MeV (in %.4f MeV steps) for l = %d\n",
            ex[0],ex[1],ex[2],l[0]);
        if (vb) printf("1 got to here\n");
    
        /*prepare to call wclbes() from libwclbes or cernlibs*/
        mu = ba[1]*ba[2]*amu/(ba[1] + ba[2]);
        rf = R0*fm*(pow(ba[1],0.33333333) + pow(ba[2],0.33333333));
        
        if (!lp)
            printf("----------------------------------------------------------------------------\n");
        else
            printf("--------------------------------------------------------------------------------------------\n");
    
        if (!lp) printf("Using R0 = %6.3f fm\n",R0);
        if (!lp) printf("  Ex  \t  l   \t   k\t\t    kr\t\t   Barrier\t     Pl\n"
                    " (MeV)\t\t (m^-1)\t\t\t\t    penet.\t   (R-mat)\n");
        else     printf("  Ex  \t  l   \t   Gamma     Brch    Barrier\t  Red Width\t Wigner lim\t Red Width/ \n"
                    " (MeV)\t\t   (keV)\t      penet.\t    (MeV)\t    (MeV)\t  Wigner\n");
        if (!lp)
            printf("----------------------------------------------------------------------------\n");
        else
            printf("--------------------------------------------------------------------------------------------\n");
        for (i = 0; i < nmex; i++)
        {
            if (!lp)
            {
                if (ernge)
                {
                    ecm = (ex[0] - thresh)*mev;
                    i = 0;
                }
                else ecm = (ex[i] - thresh)*mev;
                if (vb) printf("Masses: ba[1]=%f ba[2] = %f, thesh=%f, q=%f ex[i]=%g\n",
                    ba[1],ba[2],thresh,q,ex[i]);
                /*centre-of-mass wavenumber*/ 
	        k[i] = sqrt(2.0*mu*ecm)/hbar;
	        kr = k[i]*rf;
                /*v = hbar*k/mu;*/
                eta = (z[1]*z[2]*elec*elec*mu/(4*pi*eps0*hbar*hbar*k[i]));
        
                /*set-up variables for wclbes function call*/
                kfn = 0; mode = 1; nl = 20; zlmin = 0;
                
                if (vb) printf("2 got to here\n");
        
/*              printf("kr %f, eta %f, k %g, rf %g\n",(double)kr,(double)eta,k,rf);*/
                /*calculate the Coulomb wavefunctions, requires libwclbes.a*/
                wclbes_(&kr,&eta,&zlmin,&nl,f,g,fp,gp,sig,&kfn,&mode,&ifail,&ipr);

                if (vb) printf("3 got to here\n");
            }
            /*if ( ( l[i] == -1 && ( i != 0 || ex[0] < (tmp[0]+tmp[2]) ) ) )*/
            if (l[i] == -1 && i != 0)
                printf("----------------------------------------------------------------------------\n");
            /*print results, here j is the l[] value,
                butonly if the l values match the user input*/
            for (j = 0; j < 9; j++)
            {
                if (ernge && l[0] == -1) j = cntex;
                if (vb) printf("4 got to here\n");
                if (!lp)
                {
                    fd[j] = (double)f[j];
                    gd[j] = (double)g[j];
                    p[i][j] = 1.0/(fd[j]*fd[j]+ gd[j]*gd[j]);
                    if (vb) printf("5 got to here\n");
                }
                if (vb) printf("6 got to here\n");
                /*check if j matches an entered l value snf ig so, print*/
                if ( (j == l[i] && ernge == 0) || (j == l[0] && ernge == 1) || l[i] == -1)
                {
                    if (!lp) printf("%7.3f\t%3d\t%11g\t%9g\t%11g\t%11g\n",
                    ex[i],j,k[i],k[i]*rf,p[i][j],(double)p[i][j]*creal(kr));
                    /*above check for printing an extra tab is because the %g
                    option prints no decimal places if they are all zero*/
                    else
                    {
                        /*calculate reduced width*/
                        rwdth = gam[i]*(double)1000.0*brch[i]*elec /
                        (p[i][j]*(double)2.0*k[i]*rf);
                        /*convert to MeV*/
                        rwdth /= mev;
                        /*calculate single-particle limit (Wigner)*/
                        wig = (double)3.0*hbar*hbar / ((double)2.0*mu*rf*rf);
                        /*convert to MeV*/
                        wig /= mev;
                        printf("%7.3f\t%3d\t%8g %8g %11g\t%11g\t%11g\t%11g\n",
                            ex[i],j,gam[i],brch[i],p[i][j],rwdth,wig,rwdth/wig);
                    }
                    if (ernge && l[0] != -1) ex[0] += ex[2];
                    else if (ernge && l[0] == -1)
                    {
                        /*increment ex to next value*/
                        ex[0] += ex[2];
                        /*set j to max (=9) to got to next ex calculation with current j*/
                        j = 8;
                        /*if reached final ex value of range
                            (i.e. more than half a step above end)*/
                        if (ex[0] > (ex[1] + ex[2]/2))
                        {
                            /*reset i to zero to calculate next set of ex, and reset ex[0]*/
                            i = 0;
                            ex[0] = tmp[0];
                            /*increment counter for j*/
                            cntex++;
                            /*set j to 9 when complete and i to nmex
                                to exit both the j and i for loops*/
                            if (cntex == 9)
                            {
                                j = 9;
                                i = nmex;
                            }
                        }
                    }
                }
            }
            
            /*printf("***After calcs0, ex[0]=%f, j= %d, l[i]=%d\n",ex[0],j,l[i]);*/
            if ( (!ernge && l[i] == -1 && i != (nmex - 1)) ||
                (ernge && l[0] == -1 && ( (ex[0] < (tmp[0]+tmp[2])) || j == 10)) )
            {
                /*printf("***After calcs1\n");*/
                if (!lp)
                    printf("----------------------------------------------------------------------------\n");
                else
                    printf("--------------------------------------------------------------------------------------------\n");
            }
            /*Legacy. Not sure what this is doing since ex[0] - ex[1] should be negative*/
            if (ernge && (ex[0] - ex[1] > (ex[2]/10)))
            {
                if (vb) printf("END1: ex[0] %f, ex[1] %f\n",ex[0],ex[1]);
                printf("----------------------------------------------------------------------------\n");
                return ;
            }   
        }
        /*printf("***After calcs2\n");*/
        if (!lp) printf("----------------------------------------------------------------------------\n");
        else printf("--------------------------------------------------------------------------------------------\n");
        if (ernge == 1 || (!lp && ask_yn("Calculate reduced widths from branching ratios? [y/n] (n)",0) != 1))
            return ;
    }
} /*END coul_pen()*/
#endif

/*==========================================================================*/
/* extract_react: extract masses from input. dk = 1 = decay, 0 = reaction   */
/****************************************************************************/
int extract_react(char react[], int opt)
{
    float   tmp = 0.0, semf = 0.0;
    int     num = 0, i = 0, iold = 0, j = 0, k = 0, fmt = -1, len = 1;
    char    ans[CHLEN] = "";
    
    /*zero paramters*/
    zero_par();
   
    /*preform some preliminary checks on the reaction/decay*/
    if (opt == 0 && chk_react(react) != 0) return -1;
    if (opt == 5 && strlen(react) < 1)
    {
	printf("    Invalid isotope\n");
	return -1;
    }
    
    /*determine if reaction or decay; decay should start with left bracket*/
    if ( ! strncmp(react, "(", 1) )
    {
	dk = 1;
	printf("Decay...");
    }
    else if (opt == 0) printf("Reaction...");
    else printf("Isotope (mass look-up)...");
    
    if (dk == 0) i = 0;
    else i = 1;
    iold = i;
    /*end of decay determination*/

    /*now determine format of input*/
    format(react, &fmt);
      
    /*read isotopes/particles*/
    while ( i < (strlen(react) - 1) || (opt == 5 && i < strlen(react)) )
    {
        /*skip over characters that aren't letters or digits*/
        while (isalnum(react[i]) == 0 && i < (5*R_MAX)) i++;
        /*store current position*/
        iold = i;
        
    	/*fmt = 1 is format c12 e.g.*/
    	if (fmt == 1)
    	{
            /*count number of letters (or signs for some decay particles)*/
    	    while ( (isalpha(react[i]) != 0 ||
    		    ! strncmp(react+i, "-", 1) ||
    		    ! strncmp(react+i, "+", 1)) && i < (5*R_MAX)) i++;
    	    
            /*copy letters to el[]*/
    	    strncpy(el[k], react+iold, (i - iold));
    	}
        /*store current position*/
    	iold = i;
        /*count number of digits*/
    	while ( isdigit((int)react[i]) != 0 && i < (5*R_MAX) ) i++;
    	
        /*copy digits to a[]*/
    	for (j = iold; j < i; j++)
    	{
    	    a[k] *=10;
    	    a[k] += (int)(react[j]-'0');
    	}
        /*store current position*/
    	iold = i;
        
    	/*fmt = 2 is format 12c e.g.*/
    	if (fmt == 2)
    	{
            /*count number of letters (or signs for some decay particles)*/
    	    while ( (isalpha(react[i]) != 0 ||
    		    ! strncmp(react+i, "-", 1) ||
    		    ! strncmp(react+i, "+", 1)) && i < (5*R_MAX) )
    		i++;
    	    
            /*copy letters to el[]*/
    	    strncpy(el[k], react+iold, (i - iold));
    	}
        /*increment i to next location of input and k to next isotope number*/
    	i++; k++;
        /*store current position*/
    	iold = i;
    }
    if (fmt == 2) fmt = 1;
    
    /*print what's been read in for debug purposes*/
/*    for (i = 0; i < k; i++)
 	printf("Debug: AEl(Z) M_xs: %3d%c%c\n", a[i],el[i][0],el[i][1]);*/
   
    /*check for special notation (shorthand) such as d, t, and p etc. provided A<5*/
    for (i = 0; i < k; i++)
    {
	len = 1;
	if (isalpha(el[i][1]) != 0) len++;
	for (j = 0; j < SMAX; j++)
	{
	    if (! strncmp(el[i], sel[j], len) && a[i] < 5)
	    {
	    	strncpy(el[i], ssel[j], 2);
	    	a[i] = sa[j];
	    	z[i] = sz[j];
		break;
	    }
	}
    }
/*    printf("Before read_masses:%d,%d,%d,%s,%f,:%f:%f,%f,%f\n",
    dec[0],a[0],z[0],el[0],m[0],ba[0],dba[0],bea[0],dbea[0]);*/
    
    /*perform checks on reactions/decay and read masses*/
    /*if more than one isotope*/
    if (k > 1)
    {
	for (i = 0; i < (k + 1); i++)
    	{
	    if (i == k && opt == 0)
	    {
		/*if not decay*/
	    	if (!dk)
	    	{
    	    	    z[k] = z[0] + z[1];
    	    	    a[k] = a[0] + a[1];
	    	}
		else
		{
    	    	    z[k] = z[0];
    	    	    a[k] = a[0];
		}
	    	for (j = (2 - dk); j < k ; j++)
    	    	{
    	    	    z[k] -= z[j];
    	    	    a[k] -= a[j];
    	    	}
	    	if (fmt == 1) fmt = 0;
	    }
	    /*check for additional isotopes in reaction*/
	    if (i == k && a[k] == 0 && z[k] == 0) continue;
	    num = read_masses(i, fmt);
/*            printf("opt0 fmt = %d\n",fmt);*/
	    if (num < 0) return num;
    	}
    }
    else num = read_masses(0, fmt);
/*    printf("opt5 fmt = %d\n",fmt);*/

    if (num < 0) return num;
    if (a[k] == 0 && z[k] == 0) k--;
    
/*    printf("After read_masses:%d,%d,%d,%s,%f,:%f:%f,%f,%f\n",
          dec[0],a[0],z[0],el[0],m[0],ba[0],dba[0],bea[0],dbea[0]);*/
    
    /*reorder no.s 0 and 1, so that 0 is the beam*/
    if (!dk && !opt)
    {
    	tmp = dec[0]; dec[0] = dec[1]; dec[1] = tmp;
    	tmp = a[0]; a[0] = a[1]; a[1] = tmp;
    	tmp = z[0]; z[0] = z[1]; z[1] = tmp;
    	tmp = el[0][0]; el[0][0] = el[1][0]; el[1][0] = tmp;
    	tmp = el[0][1]; el[0][1] = el[1][1]; el[1][1] = tmp;
    	tmp = m[0]; m[0] = m[1]; m[1] = tmp;
    	tmp = dm[0]; dm[0] = dm[1]; dm[1] = tmp;
    	tmp = ba[0]; ba[0] = ba[1]; ba[1] = tmp;
    	tmp = dba[0]; dba[0] = dba[1]; dba[1] = tmp;
    	tmp = bea[0]; bea[0] = bea[1]; bea[1] = tmp;
    	tmp = dbea[0]; dbea[0] = dbea[1]; dbea[1] = tmp;
    }
    /*printf("k+1=%d\n",k+1);*/
    /*print isotope properties*/
    printf("\n Masses and binding energies in keV\n");
    for (i = 0; i < (k + 1); i++)
    {
	/*if necessary, convert first element letter to capital*/
	if (el[i][0] > 96 && dec[i] != 1) el[i][0] -= 32;
	/*if necessary, uncapitalise second element letter*/
	if (el[i][1] < 91 && el[i][1] > 64) el[i][1] += 32;
	/*for neutrons don't capitalise*/
	if (a[i] == 1 && el[i][0] == 78) el[i][0] += 32;      
        /*set any trailing spaces ((char)32) to null characters (char)0
            in element names*/
        if (el[i][1] == (char)32) el[i][1] = (char)0;
	
        if (dec[i] != 1) printf("AEl(Z) M_xs: %3d%c%c(%3d) %13.5f +/- %.5f %11.6f +/- %.6f\n",
		a[i],el[i][0],(isalpha(el[i][1]) ? el[i][1] : ' '),
	    	z[i],m[i],dm[i],ba[i],dba[i]);
	else printf("Decay particle: %c%c\n",el[i][0],el[i][1]);
        /*if 2-body reaction (i.e. k+1 == 4 && a[2]>a[3] swap 3 and 4
            for kinematics calculation, so first kinematic solution is as expected
            with lightest particle taking most energy at forward angles.
            Ask this only after processing all nuclei/particles*/  
        if (i == k && !dk && !opt && (k+1) == 4 && a[2] > a[3])
        {
            /*ask user if they want to swap ejectile and recoil order
                to follow usual convention of lightest particle first*/
            sprintf(ans,"%sSwap order of %d%c%c and %d%c%c so lightest particle is 1st product [y/n] (y)?%s\n"
                "%s[For kinematics, the angle entered corresponds to that of the 1st reaction product.]%s",
                clr[3],a[2],el[2][0],(isalpha(el[2][1]) ? el[2][1] : (char)(3)),
                a[3],el[3][0],(isalpha(el[3][1]) ? el[3][1] : (char)(3)),clr[0],clr[3],clr[0]);
            if (ask_yn(ans,1))
            {
    	        tmp = dec[2]; dec[2] = dec[3]; dec[3] = tmp;
    	        tmp = a[2]; a[2] = a[3]; a[3] = tmp;
    	        tmp = z[2]; z[2] = z[3]; z[3] = tmp;
    	        tmp = el[2][0]; el[2][0] = el[3][0]; el[3][0] = tmp;
    	        tmp = el[2][1]; el[2][1] = el[3][1]; el[3][1] = tmp;
    	        tmp = m[2]; m[2] = m[3]; m[3] = tmp;
    	        tmp = dm[2]; dm[2] = dm[3]; dm[3] = tmp;
    	        tmp = ba[2]; ba[2] = ba[3]; ba[3] = tmp;
    	        tmp = dba[2]; dba[2] = dba[3]; dba[3] = tmp;
    	        tmp = bea[2]; bea[2] = bea[3]; bea[3] = tmp;
    	        tmp = dbea[2]; dbea[2] = dbea[3]; dbea[3] = tmp;
            }
        }
        
        /*for mass-look-up print SEMF and measured B.E.s*/
	if (opt == 5 && dec[i] != 1)
	{
	    semf = (15.8*a[i]) - (18.3*pow(a[i],0.6666666)) -
		(0.714*z[i]*(z[i]-1))/pow(a[i],0.33333) -
		    (23.2*(a[i]-(2*z[i]))*(a[i]-(2*z[i]))/a[i]);
	    if (z[i] % 2 == 0 && (a[i] - z[i]) % 2 == 0) semf += 12/sqrt(a[i]);
	    else if (z[i] % 2 != 0 && (a[i] - z[i]) % 2 != 0)
		semf -= 12/sqrt(a[i]);	
	    /*convert to keV*/
	    semf *= 1000;
	    
	    printf("   SEMF: %f keV, Meas: %12.3f +/- %.3f keV, Meas-SEMF: %10.3f keV\n",
		    semf,bea[i]*a[i],dbea[i]*a[i],((bea[i]*a[i])-semf));
            printf("   BE/A, Meas: %12.3f +/- %.3f keV\n",bea[i],dbea[i]);
	}
    }
    return (k + 1);
} /*END extract_react()*/

/*==========================================================================*/
/* format: perform preliminary checks on reaction/decay entered   	    */
/****************************************************************************/
void format(char react[], int *fmt)
{
    int     i = 0;
    
    /*find position of first digit*/
    while (isdigit((int)react[i]) == 0 && i < (int)strlen(react)) i++;
    
    if (i < (int)strlen(react))
    {
    	/*fmt = 2 is format 12c e.g.*/
	if (i == 0 && (isalpha(react[1]) != 0 || isalpha(react[2]) != 0 ||
		isalpha(react[3]) != 0) ) *fmt = 2;
	else if (i > 0 && isalnum(react[i-1]) == 0 && isalnum(react[i+1]) != 0)
                    *fmt = 2;
        /*fmt = 1 is format c12 e.g.*/
	else if ( i > 0 && isalpha(react[i-1]) != 0 &&
                  ((isalnum(react[i+1]) == 0) || (isdigit(react[i+1]) != 0)) )
                    *fmt = 1;
        else
        {
            printf("****ERROR...input format not detected\n");
            return ;
        }
    }
    else *fmt = 1;
    printf("input has format %d\n",*fmt);
} /*END format()*/

/*==========================================================================*/
/* kin: calculate reaction kinematics	    	    	    	    	    */
/****************************************************************************/
int kin(char react[], int opt)
{
    static float ejf = 0.0, ext = 0.0, e0 = 0.0;
    static double tf[3];
    double  e[4], ec[4], p[4], etot = 0.0, ec_in = 0.0, ec_out = 0.0;
    double  ex = 0.0, rt = 0.0, tlo = 0.0, thi = 0.0, tstep = 0.0;
    double  ecr[4], tcr[4], v[4], vc[4], vc_in = 0.0, vc_out = 0.0, qv = 0.0;
    double  d = 0.0, con = 0.0;
    double  vcr_in = 0.0, vcr_out = 0.0, gam_vcr_in = 0.0, gam_vcr_out = 0.0;
    double  tcr_in = 0.0, tcr_out = 0.0, vcr[4], ri = 0.0, rf = 0.0;
    double  bar[4], er[4], pr[4], etotr = 0.0, tr[4];
    double  vr[4], phir = 0.0, phird = 0.0, e2m;
    double  theta = 0.0, thetac = 0.0, thetad = 0.0, thetacr = 0.0;
    double  phi = 0.0, phid = 0.0, rsxs = 0.0, thcutoff = 0.0;
    double  thcutoffr = 0.0, mui, muf;
    double  coul = 0.0, coulc = 0.0, thgraz = 0.0, thgrazc = 0.0;
    double  pi = 3.1415926535, bim = 0.0, rmin = 0.0;
    double  arr = 0.0, brr = 0.0, crr = 0.0, ej = 0.0, lin, lout;
    double  pc_in = 0.0, pc_out = 0.0;
/*    double  ecr_in = 0.0, ecr_out = 0.0, dca, gam[4], pcr[4];
    double  ar = 0.0, br = 0.0, cr = 0.0, phic = 0.0, phicr = 0.0;*/
    int     flg = 0, i = 0, j = 0, nnr = 0, n2s = 0, nrs = 0;
    char    ans[CHLEN] = "";
    
    e2m = c*c;
    e2m *= amu;
    e2m /= mev;
    
    for (i = 0; i < 4; i++)
    {
	bar[i] = 0.0;
/*	gam[i] = 0.0;*/
	e[i] = 0.0;
	er[i] = 0.0;
	ec[i] = 0.0;
	ecr[i] = 0.0;
	p[i] = 0.0;
	pr[i] = 0.0;
/*	pcr[i] = 0.0;*/
	tcr[i] = 0.0;
	tr[i] = 0.0;
	v[i] = 0.0;
	vr[i] = 0.0;
	vc[i] = 0.0;
	vcr[i] = 0.0;
    }
    if (e0 == 0.0 ) for (i = 0; i < 3; i++) tf[i] = 0.0;
    
    if (opt == 3 || opt == 4)
    {
        printf("Enter the beam energy in the lab. [MeV]");
        printf(" [<Enter> for %.2f]\n",e0);
        get_line(ans, CHLEN);
        if (strlen(ans) == (int)0 && (ans[0] == ((char)0)));
        else e0 = atof(ans);
        e[0] = (double)e0;

        if (opt == 3)
        {
	    printf("Enter scattering angle (lab), or range (low high step)"
	        " [deg.]\n");
    	    printf(" [<Enter> for %.3f %.3f %.3f]\n",tf[0],tf[1],tf[2]);
        }
        else if (opt == 4)
        {
	    printf("Enter scattering angle (lab), [deg.]\n");
    	    printf(" [<Enter> for %.3f]\n",tf[0]);
        }
        /*extract 3 numbers from string ans0*/
        get_line(ans, CHLEN);
        if (strlen(ans) == (int)0 && (ans[0] == ((char)0)));
        else if ( (i = get_pars(ans, tf, 3)) < 3)
        {
	    if (i == 2) tf[2] = 1.0;
	    else if (i == 1)
	    {
	        tf[2] = 1.0;
	        tf[1] = tf[0];
	    }
        }
        if (opt == 3) printf(" -->Angles entered: %.3f %.3f %.3f\n",
	    tf[0],tf[1],tf[2]);
        else if (opt == 4) printf(" -->Angle entered: %.3f\n",tf[0]);
    
        tlo = (double)tf[0];
        thi = (double)tf[1];
        tstep = (double)tf[2];
        if (tstep == 0.0) tstep = 5.0;
    }
    
    if (flg < 0) return flg;
    /*convert the reaction Q-value, q, to MeV*/
    qv = (double)q/1000;
    
    /*usual two-body kinematics (opt == 6 is resonant reaction)*/
    if (opt == 3 || opt == 6)
    {
	printf("Enter the excitation energy [MeV]");
    	printf(" [<Enter> for %.3f]\n",ext);
    	get_line(ans, CHLEN);
    	if (strlen(ans) == (int)0 && (ans[0] == ((char)0)));
    	else ext = atof(ans);
    	ex = (double)ext;
        /*resonant reaction beam*/
        if (opt == 6)
        {
            /*calculate the beam energy*/
            e[0] = (ex - qv)*(ba[1]+ba[0])/ba[1];
	    printf("\n The %d%c%c beam energy is = %.3f MeV "
                "to reach %.3f MeV in %d%c%c*\n",
                a[0],el[0][0],el[0][1],e[0],ex,a[2],el[2][0],el[2][1]);
	    return 0;
        }
    }
    else if (opt == 4) /*solve for ex knowing ejectile energy*/
    {
	printf("Enter the ejectile energy [MeV]");
    	printf(" [<Enter> for %.3f]\n",ejf);
    	get_line(ans, CHLEN);
    	if (strlen(ans) == (int)0 && (ans[0] == ((char)0)));
    	else ejf = atof(ans);
    	ej = (double)ejf;
     	p[0] = sqrt((double)2.0*ba[0]*e[0]);
   	p[2] = sqrt((double)2.0*ba[2]*ej);
	p[3] = sqrt( (p[0]*p[0]) + (p[2]*p[2]) -
		((double)2.0*p[0]*p[2]*cos( tlo*pi/(double)180.0 )) );
    	e[3] = (p[3]*p[3])/((double)2.0*ba[3]);
	ex = e[0] + qv - ej - e[3];
	printf("\n Excitation energy = %.3f MeV\n",ex);
	return 0;
    }
    
    for (i = 0; i < 4; i++) bar[i] = ba[i]*e2m;
    
    etot = e[0] + qv - ex;
        
    tr[0] = e[0];
    er[0] = tr[0] + bar[0];
    er[1] = bar[1];
    etotr = er[0] + bar[1] - ex;
/*    printf("etotr - m = %f\n",(etotr - bar[2] - bar[3]));*/
    
    /*beam energy in COM frame*/
    ec[0] = e[0]*ba[1]/(ba[0] + ba[1]);
    /*Coulomb barrier*/
/*    ri = (double)R0*(pow(a[0],0.33333333)+pow(a[1],0.33333333)+(double)2.0);*/
    ri = (double)R0*(pow(ba[0],0.33333333)+pow(ba[1],0.33333333));
/*    rf = (double)R0*(pow(a[2],0.33333333)+pow(a[3],0.33333333)+(double)2.0);*/
    rf = (double)R0*(pow(ba[2],0.33333333)+pow(ba[3],0.33333333));
    coulc = ((double)1.44*z[0]*z[1])/ri;
    coul = ( (double)1.0 + (ba[0]/ba[1]) )*coulc;
    /*grazing angle*/
    thgrazc = ( ((double)2.0*ec[0]/coulc) - (double)1.0 );
    if (thgrazc > 1.0)
    {
    	thgrazc = (double)2.0*asin( (double)1.0 /
	    	    ( ((double)2.0*ec[0]/coulc) - (double)1.0 ) );
    	thgraz = atan( sin(thgrazc)/( (ba[0]/ba[1]) + cos(thgrazc) ) );
    }
    else thgrazc = -1.0;

    /*ec_in: initial energy of centre-of-mass frame*/
    ec_in = e[0]*ba[1]/(ba[0] + ba[1]);
    /*ec_out: final energy of centre-of-mass frame*/
    ec_out = ec_in + qv - ex;
    /*print a summary of the initial information*/
    printf("\n Kinetic energy\t\t\t: %g MeV\n",e[0]);
    printf(" Coulomb barrier \t\t: (CM) %g : (LAB) %g MeV\n",coulc,coul);
    printf(" Centre-of-mass energy (in)\t: %g MeV\n",ec_in);
    printf(" Q-value\t\t\t: %g +/- %7.5f MeV\n",qv,(double)(dq/1000));
    printf(" Kinetic energy (out)\t\t: %g MeV\n",etot);
    printf(" Excitation energy\t\t: %g MeV\n",ex);
    if (thgrazc != -1.0)
    	printf(" Grazing angle\t\t\t: (CM) %g : (LAB) %g [deg.]\n",
	    thgrazc*((double)180.0/pi),thgraz*((double)180.0/pi));
    else printf(" Grazing angle\t\t\t: below the barrer!\n");
        
    /*reduced masses before and after reaction*/
    mui = ba[0]*ba[1]/(ba[0] + ba[1]);
    muf = ba[2]*ba[3]/(ba[2] + ba[3]);
    
    /*incoming and outgoing momenta and angular momenta*/
    pc_in = sqrt(2.0*ec_in*mev*mui*amu);
    pc_out = sqrt(2.0*ec_out*mev*muf*amu);
    lin = pc_in*ri*fm/hbar;
    lout = pc_out*rf*fm/hbar;

    thetad = tlo;    
    while ( thetad <= thi )
    {
	/*calculate some coefficients for quadratic equations*/
	theta = (double)(pi/180)*thetad;

	p[0] = sqrt((double)2.0*ba[0]*e[0]);
	
	pr[0] = sqrt((er[0]*er[0]) - (ba[0]*ba[0]*e2m*e2m));
	
	/*constants for solving the quadratic equation for pr[2]*/
	con = (etotr*etotr) - (pr[0]*pr[0]) - (bar[3]*bar[3]) + (bar[2]*bar[2]);
	arr = ((double)4*pr[0]*pr[0]*cos(theta)*cos(theta)) -
		((double)4*etotr*etotr);	
	brr = (double)4*con*pr[0]*cos(theta);	
	crr = (con*con) - ((double)4*etotr*etotr*bar[2]*bar[2]);
	
		
	/*constants for solving the quadratic equation for er[2]*/
/*	ar = ((double)4*pr[0]*pr[0]*cos(theta)*cos(theta)) -
		((double)4*etotr*etotr);	
	br = ( ((double)4*pow(etotr,3)) - ((double)4*pr[0]*pr[0]*etotr)
		+ ((double)4*bar[2]*bar[2]*etotr)
		- ((double)4*bar[3]*bar[3]*etotr) );	
	cr = (double)-4*bar[2]*bar[2]*pr[0]*pr[0]*cos(theta)*cos(theta);
	cr -= ( pow(bar[3],4) + pow(bar[2],4) + pow(pr[0],4) + pow(etotr,4) );
	cr += (double)2*bar[3]*bar[3]*bar[2]*bar[2]
		- (double)2*bar[3]*bar[3]*pr[0]*pr[0]
		+ (double)2*bar[3]*bar[3]*etotr*etotr
		+ (double)2*bar[2]*bar[2]*pr[0]*pr[0]
		- (double)2*bar[2]*bar[2]*etotr*etotr
		+ (double)2*pr[0]*pr[0]*etotr*etotr;*/
	
	/*rt is the sqrt part of the quadratic solution for er[2]*/
	rt = pow( (p[0]*cos(theta)), (double)2.0 );
	rt -= ((double)1.0 + (ba[3]/ba[2]))*( (p[0]*p[0])
		- ((double)2.0*ba[3]*etot) );
	rt = sqrt(rt);
	
	/*work out maximum scattering angle*/
	if (thetad == tlo)
	{
	    thcutoff = ((double)1.0 + (ba[3]/ba[2]))*( (p[0]*p[0])
		- ((double)2.0*ba[3]*etot) );
	    thcutoff /= p[0]*p[0];
	    
	    if (thcutoff >= -1.0 && thcutoff <= 1.0)
		thcutoff = acos(sqrt(thcutoff));
	    else thcutoff = pi;
	    printf(" Theta cut off (non-Rel.)\t:%7.3f [deg.]\n",
		    thcutoff*((double)180.0/pi));
	
	    /*get 3 coefficients and solve for thetacutoffr*/
	    /*thcutoffr here is cos^2(thcutoffr)*/
	    thcutoffr = sqrt( crr/(crr -(con*con)) );
	    thcutoffr *= (etotr/pr[0]);
	    if (thcutoffr >= -1.0 && thcutoffr <= 1.0)
		thcutoffr = acos(thcutoffr);
	    else thcutoffr = pi;
	    printf(" Theta cut off (Rel.)\t\t:%7.3f [deg.]\n",
		    thcutoffr*((double)180.0/pi));
    	}
        
	if (thetad == tlo)
	{
	    nnr = ask_yn("Suppress non-relativistic solutions? [y/n] (n)",0);
            if (nnr == 0)
	        nrs = ask_yn("Suppress relativistic solutions? [y/n] (y)",1);
            
	    n2s = ask_yn("Suppress second kinematic solution? [y/n] (y)",1);
	}
        
	if ( (theta > thcutoff && nnr == 0) || (theta > thcutoffr && nrs == 0) )
	{
	    printf("\t**** Maximum scattering angle exceeded ****\n");
	    break;
	}        
        
	/*There can be 2 solutions*/
	for (j = 0; j < 2; j++)
	{    	    
	    p[2] = p[0]*cos(theta);
	    if (j == 0)
	    {
		p[2] += rt;
	      	pr[2] = (double)-1*brr - sqrt(brr*brr - ((double)4*arr*crr));
	    }
	    else 
	    {
		p[2] -= rt;
	      	pr[2] = (double)-1*brr + sqrt(brr*brr - ((double)4*arr*crr));
	    }
	    p[2] /= (double)1.0 + (ba[3]/ba[2]);
	    pr[2] /= (double)2*arr;
	    
	    /*get signs correct after sqrt*/
	    if (p[2] < 0.0) pr[2] *= -1.0;
	    
	    if ( p[2] < 0.0 && n2s == 0 && j == 1 )
	    {
		printf("\tNo second kinematic solution.\n");
		n2s = 1;
	    }
	    
	    if ( pr[2] < 0.0 && nnr == 0 && j == 1)
	    {
		printf("\tNo second kinematic solution.\n");
		nnr = 1;
	    }
	    	    
	    /*for kinetic energy (tr) subract rest mass energy*/
	    er[2] = sqrt( (pr[2]*pr[2]) + (bar[2]*bar[2]) );
	    tr[2] = er[2] - bar[2];
	    er[3] = etotr - er[2];
	    tr[3] = er[3] - bar[3];
	    
	    e[2] = (p[2]*p[2])/((double)2.0*ba[2]);
	    e[3] = etot - e[2];
	
	    p[3] = sqrt((double)2.0*ba[3]*e[3]);
	    pr[3] = sqrt(er[3]*er[3] - bar[3]*bar[3]);	    
	    
	    if (p[3] < 0.0) pr[3] *= -1.0;
	    
	    phi = asin( (p[2]/p[3])*sin(theta) );
	    phir = asin( (pr[2]/pr[3])*sin(theta) );
            /*check for sign of phi*/
            if ( (p[0] - p[2]*cos(theta)) < 0.0 && phi < 90.0 && phi > 0.0 )
            {
                phi = pi - phi;
                phir = pi - phir;
            }
            
	    phid = ((double)180.0/pi)*phi;
	    phird = ((double)180.0/pi)*phir;
	    	    
	    /*vr are the relativistic speeds/c, i.e. beta*/
	    for (i = 0; (i < 4) && (er[i] != (double)0.0); i++)
	    {
		vr[i] = pr[i]/er[i];
/*		gam[i] = (double)1/sqrt((double)1 - (vr[i]*vr[i]));*/
	    }
	    
	    vcr_in = pr[0]/etotr;
	    gam_vcr_in = (double)1/sqrt((double)1 - (vcr_in*vcr_in));
	    	    
	    d = (bar[0] + bar[1])/(bar[2] + bar[3]) ;
	    d *= (bar[0] + bar[1])/(bar[2] + bar[3]) ;
	    d *= gam_vcr_in*gam_vcr_in*vcr_in*vcr_in;	    
	    vcr_out = sqrt(d/( (double)1 + d) );
	    gam_vcr_out = (double)1/sqrt((double)1 - (vcr_out*vcr_out));
	    	    
	    ecr[0] = gam_vcr_in*(er[0] - (pr[0]*vcr_in));
	    tcr[0] = ecr[0] - bar[0];
	    ecr[1] = gam_vcr_in*(er[1] - (pr[1]*vcr_in));
	    tcr[1] = ecr[1] - bar[1];
	    ecr[2] = gam_vcr_out*(er[2] - (pr[2]*vcr_out*cos(theta)));
	    tcr[2] = ecr[2] - bar[2];
	    ecr[3] = gam_vcr_out*(er[3] - (pr[3]*vcr_out*cos(phir)));
	    tcr[3] = ecr[3] - bar[3];
	    
	    tcr_in = tcr[0] + tcr[1];
/*	    tcr_out = tcr_in + qv - ex;*/
	    tcr_out = tcr_in + bar[0] + bar[1] - bar[2] - bar[3] - ex;
/*	    ecr_in = tcr_in + bar[0] + bar[1];
	    ecr_out = tcr_out + bar[2] + bar[3];*/
	    	    	    
	    vcr[0] = (vr[0] - vcr_in)/((double)1 - (vr[0]*vcr_in));
	    vcr[1] = (vr[1] - vcr_in)/((double)1 - (vr[1]*vcr_in));	    
	    vcr[2] = sqrt( (double)1 - ((bar[2]*bar[2])/(ecr[2]*ecr[2])) );
	    vcr[3] = sqrt( (double)1 - ((bar[3]*bar[3])/(ecr[3]*ecr[3])) );
	    vcr[3] *= (double)-1;
	    
/*	    pcr[0] = ecr[0]*vcr[0];
	    pcr[1] = ecr[1]*vcr[1];
	    pcr[2] = ecr[2]*vcr[2];
	    pcr[3] = ecr[3]*vcr[3];*/
            
/*	    pcr[3] = (double)-1*pcr[2];*/
	    
	    
	    thetacr = (pr[2]*sin(theta));
	    thetacr /= gam_vcr_out*( (pr[2]*cos(theta)) - (vcr_out*er[2]) );
	    thetacr = atan(thetacr);
/*	    thetacr = (vr[2]*cos(theta) - vcr_out)/((double)1 - (vr[2]*cos(theta)*vcr_out));
	    thetacr = acos(thetacr/vcr[2]);
	    thetacr = (vr[3]*cos(phir) - vcr_out)/((double)1 - (vr[3]*cos(phir)*vcr_out));
	    thetacr = acos(thetacr/vcr[3]);*/
	    if (thetacr < (double)0.0) thetacr += pi;
	    
/*	    phicr = pi - thetacr;*/
		    
	    for (i = 0; i < 4; i++) v[i] = p[i]/ba[i];
	    
	    /*calculate com frame values*/
	    ec[0] = ba[1]/(ba[0] + ba[1]);
	    ec[0] *= ba[1]/(ba[0] + ba[1]);
	    ec[0] *= e[0];
	    
	    ec[1] = ba[0]/(ba[0] + ba[1]);
	    ec[1] *= ba[1]/(ba[0] + ba[1]);
	    ec[1] *= e[0];
	    
	    ec[2] = ba[3]/(ba[2] + ba[3]);
	    ec[2] *= ec_out;
	    
	    ec[3] = ba[2]/(ba[2] + ba[3]);
	    ec[3] *= ec_out;
	    	    
	    vc_in = (v[0]*ba[0])/(ba[0] + ba[1]);
/*	    vc_out = vc_in*(ba[0] + ba[1])/(ba[2] + ba[3]);*/
	    vc_out = (v[0]*ba[0])/(ba[2] + ba[3]);
	    
	    vc[0] = v[0] - vc_in;
	    vc[1] = v[1] - vc_in;
	    vc[2] = sqrt((double)2.0*ec[2]/ba[2]);
	    vc[3] = (double)-1.0*(ba[2]/ba[3])*vc[2];
	    	
/*	    thetac = asin( (sin(theta))*sqrt(e[2]/ec[2]) );*/
/*	    thetac = asin( sin(theta)*(v[2]/vc[2]) );*/
/*	    thetac = atan( (v[2]*sin(theta))/(v[2]*cos(theta) - vc_out) );*/    
	    thetac = acos( (v[2]*cos(theta) - vc_out)/vc[2] );	    
/*    	    phic = pi - thetac;*/
	    
	    /*twokin uses the following equations for vc[2], but is the
	    same as that calculated here. But for calculating
	    rsxs, twokin uses 0.5*mu*vc[2]^2 as the energy term. This is
	    presumably to correct for finite recoil mass. Here, the classical
	    Rutherford xs is calculated assuming infinite target mass. The
	    two converge as target mass increases. Note the twokin formula
	    needs transforming to the lab frame via rsxslab
	    However, I think the twokin calculations are wrong*/
/*    	    vc[2] = sqrt( (ba[3]/(ba[2]*(ba[0]+ba[1])))
		    *2.0*(( ( ba[1]/(ba[0]+ba[1]) )*e0)+qv) );*/
/*    	    rsxs = z[0]*z[1];
	    rsxs /= vc[2]*vc[2]*(mu/2.0)*sin(thetac/(double)2.0)*sin(thetac/(double)2.0);*/
	    rsxs = z[0]*z[1];
	    rsxs /= (e0)*sin(theta/(double)2.0)*sin(theta/(double)2.0);
	    rsxs *= rsxs;
	    rsxs *= (double)1.3;
/*	    rsxslab *= ( sin(thetac)*sin(thetac) )/
    	    	    	    ( sin(theta)*sin(theta)*cos(thetac-theta) );*/
	    /*already in mb, not fm^2, so don't multiply by 10*/
	    
	    /*distance of closest approach (Classical) p397 Krane*/
	    /*bim is impact parameter*/
/*	    dca = 1.44*z[0]*z[1]/e0;*/
	    bim = (((1.44/2.0)*z[0]*z[1])/e0)*sqrt( (1+cos(theta))/(1-cos(theta)) );
	    rmin = bim*cos(theta/2.0)/(1-sin(theta/2.0));
	    
	    /*if first time round loop*/
	    if (thetad == tlo && j == 0)
	    {
	    	printf(" v_cm(in)/c (rel)\t: %.4f (non-rel: %.4f)\n"
                       " v_cm(out)/c (rel)\t: %.4f (non-rel: %.4f)\n",vcr_in,
                        (vc_in*sqrt(mev/amu)/c),vcr_out,(vc_out*sqrt(mev/amu)/c));
	    	printf(" ec_cm(in)   \t\t: %g (non-rel: %g) MeV\n"
                       " ec_cm(out)  \t\t: %g (non-rel: %g) MeV\n",tcr_in,ec_in,tcr_out,ec_out);
	    	printf(" ang. mom. L(in) \t: %g hbar\n"
                       " ang. mom. L(out)\t: %g (deltaL: %g) hbar\n\n",
                        lin,lout,lin-lout);
                
                /*print reaction to screen*/
                qval(react, 3);

	    	printf("\t                                                   Ruth. Scatt.  Dist. clos.\n");
	    	printf("\ttheta2(L)  theta2(CM)   E2(L)    E3(L)    theta3(L) dsig/dOmeg   approach, d\n");
	      	printf("\t [deg.]      [deg.]     [MeV]    [MeV]     [deg.]    [mb/sr]        (fm)\n");
	    }
	    /*if want non-relativistic solutions, j=0*/
	    if (j == 0 && nnr == 0)
	    	printf("Sol.%d: %8.3f   %8.3f   %8.3f  %8.3f  %8.3f   %11g  %8.3f\n",
		    j+1,thetad,(thetac*(180/pi)),e[2],e[3],phid,rsxs,rmin);
	    /*if want second non-relativistic solutions, j=1*/
	    else if (nnr == 0 && j == 1 && n2s == 0)
	    	printf("Sol.%d:            %8.3f   %8.3f  %8.3f  %8.3f\n",
		    j+1,(thetac*(180/pi)),e[2],e[3],phid);
		
	    /*if want relativistic solutions, j=0*/
	    if (j == 0 && nrs == 0)
	    	printf("Rel.%d: %8.3f   %8.3f   %8.3f  %8.3f  %8.3f   %11g  %8.3f\n",
		    j+1,thetad,(thetacr*(180/pi)),tr[2],tr[3],phird,rsxs,rmin);
	    /*if want second relativistic solutions j=1*/
	    else if (j == 1 && n2s == 0 && nrs == 0)
		printf("Rel.%d:            %8.3f   %8.3f  %8.3f  %8.3f\n",
		j+1,(thetacr*(180/pi)),tr[2],tr[3],phird);
	    if (j == 1 && (n2s == 0 || (nnr == 0 && nrs == 0.0)) )
                    printf("-------\n");
	}
	thetad += tstep;
    }
    return 0;
} /*END kin()*/

/*==========================================================================*/
/* ask_yn: 0 means no, 1 means yes. yn == -1 for no default, else yn default*/
/****************************************************************************/
int ask_yn(char *questn, int yn)
{
    char ans[3] = "";
    
    while (1)
    {
    	printf("%s\n", questn);
        /*0 option to get_ans doesn't allow carriage returns*/
        if (yn == -1) get_ans(ans,1,0);
        else get_ans(ans,1,1);
        if (yn != -1 && (ans[0] == '\n')) return yn;
        else if (ans[0] == 'y' || ans[0] == 'Y') return 1;
	else if (ans[0] == 'n' || ans[0] == 'N') return 0;
    }
} /*END ask_yn()*/

/*==========================================================================*/
/* chk_react: perform preliminary checks on reaction/decay entered   	    */
/****************************************************************************/
int chk_react(char react[])
{
    int     bd = 0;
    /*perform some preliminary checks on the reaction entered*/
    /*min. decay length 5, e.g. (a,g) & min. reaction length 6, e.g. a(d,g)*/
    if ( (strlen(react) < 6 && react[0] != '(') ||
         (strlen(react) < 5 && react[0] == '(') )
    {
    	bd = 1;
    	printf(" Reaction length < 6 or decay length < 5 length characters...\n");	
    }
    /*check for left and right brackets with comma between,
        first check if all are present*/
    else if ( ! strchr(react,'(') || ! strchr(react,',') || ! strchr(react,')') )
    {
        bd = 1;
        printf(" One left bracket, one right bracket and one comma not found...\n");
    }
    /*check that comma is at least two array elements from each bracket*/
    else if ( ((int)(strchr(react,',') - strchr(react,'(')) <= 1) ||
              ((int)(strchr(react,')') - strchr(react,',')) <= 1) )
    {
        bd = 1;
        printf(" Spacing between brackets and comma is not correct\n");
    }
    /*check for letter/number number in first or second elements of input*/
    else if (isalnum((int)react[0]) == 0 && isalnum((int)react[1]) == 0)
    {
    	bd = 1;
    	printf(" First and second characters are not alphanumeric...\n"); 	
    }
    
    if (bd == 1) printf("***Bad reaction or decay. Check input***\n");

    return bd;
} /*END chk_react()*/

/*==========================================================================*/
/* decay_part: load decay particle into arrays	    	    	    	    */
/****************************************************************************/
void decay_part(void)
{
    extern  float  decm[DMAX];
    extern  int  decz[DMAX], deca[DMAX];
    extern  char decel[DMAX][2];
    int     i;
    
    /*store decay particles in arrays*/
    /*gamma-ray*/
    memcpy(decel[0], "g", 1);
    decz[0] = 0;
    deca[0] = 0;
    decm[0] = 0;
    /*beta-*/
    memcpy(decel[1], "e-", 2);
    decz[1] = -1;
    deca[1] = 0;
    decm[1] = 0.0;
    memcpy(decel[2], "b-", 2);
    decz[2] = -1;
    deca[2] = 0;
    decm[2] = 0.0;
    /*beta+*/
    memcpy(decel[3], "e+", 2);
    decz[3] = 1;
    deca[3] = 0;
    decm[3] = 2*510.998918;
/*    decm[3] = 0.0;*/
    memcpy(decel[4], "b+", 2);
    decz[4] = 1;
    deca[4] = 0;
    decm[4] = 2*510.998918;
/*    decm[4] = 0.0;*/
    /*electron capture*/
    memcpy(decel[5], "ec", 2);
    decz[5] = 1;
    deca[5] = 0;
    decm[5] = 0.0;
     /*2beta+*/
    memcpy(decel[6], "ee", 2);
    decz[6] = -2;
    deca[6] = 0;
    decm[6] = 0.0;
    memcpy(decel[7], "bb", 2);
    decz[7] = -2;
    deca[7] = 0;
    decm[7] = 0.0;
    
    printf("Allowed particles/decay modes are: ");
    for (i = 0; i < DMAX; i++) printf(" %c%c%c",decel[i][0],decel[i][1],
	 (DMAX-i-1 ? ',' : ' '));
    printf("\n");
} /*END decay_part()*/

/*==========================================================================*/
/* get_ans: get answer which can be carriage return if md = 1. 	    	    */
/****************************************************************************/
void get_ans(char ans[], int num, int md)
{
    int     i;
    struct  termios newt, oldt;
    
    /*This code is a modified version of a function in the RadWare
        software suite, by D.C. Radford*/
    while (1)
    {
    	tcgetattr(0, &oldt);
	newt = oldt;
    	newt.c_lflag &= ~ICANON;
    	newt.c_cc[VMIN] = 1;
    	newt.c_cc[VTIME] = 0;
	/*handle sigs*/
/*    	newt.c_lflag |= ISIG;*/
    	tcsetattr(0, TCSANOW, &newt);
    	i = 0;
    	while( (ans[i++] = (char)getchar()) != '\n' && i < num) ;
    	
	tcsetattr(0, TCSANOW, &oldt);
    	if (ans[i-1] != '\n') printf("\n");
	else if (md == 0 && ans[0] == '\n') continue;
	
    	ans[i] = '\0';
    	return ;
    }
} /*END get_ans()*/

/*==========================================================================*/
/* get_line: read a line into ans, of max length len   	    	    	    */
/****************************************************************************/
void get_line(char ans[], int len)
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(char)*len);
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = getchar()) != '\n' && c != EOF && i < len) ans[i++] = c;
    
    ans[i] = '\0';
} /*END get_line()*/

/*==========================================================================*/
/* get_mode: get mode from user input	    	    	    	    	    */
/****************************************************************************/
int get_mode(int md)
{
    char    ans[10] = "";
    
    while(1)
    {
    	printf("\n 1) Enter/change a reaction or decay\n");
    	printf(" 2) Calculate the reaction Q value\n");
    	printf(" 3) Work out some kinematics\n");
    	printf(" 4) Input ejectile energy and solve for excitation energy\n");
    	printf(" 5) Mass look-up and SEMF/B.E.s\n");
    	printf(" 6) Calculate resonant-reaction beam energy for a given excitation energy\n");
#ifdef HAVE_WCLBES
    	printf(" %1d) Calculate penetrabilities and reduced widths\n",NUMOPT);
#endif
    	printf(" 0) Quit\n");
	
	get_ans(ans,1,0);
 	if (ans[0] - '0' >= 0 && ans[0] - '0' <= NUMOPT)
	{
	    md = ans[0] - '0';
	    return md;
	}
    }
} /*END get_mode()*/

/*==========================================================================*/
/* get_pars: extract comma or space separated numbers from string ans0	    */
/****************************************************************************/
int get_pars(char ans0[], double pars[], int num)
{
    int     i, j = 0, k, minus, p = 0;
    char    ans1[200] = "";

    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
	memset(ans1,'\0',sizeof(ans1));
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	if (k > 0) p++;
	j++;
	pars[i] = (double)atof(ans1);
    }
    return p;
} /*END get_pars()*/
	
/*==========================================================================*/
/* qval: calculate Q-value. If md == 1, qvalue and reaction will be printed */
/****************************************************************************/
int qval(char react[], int md)
{
    extern float   q, dq;
    int     i = 0, num = 0;
    
    if (md != 3) num = extract_react(react, 0);
    else num = 4;
    if (num < 0)
    {
	printf("***ERROR. Invalid reaction.\n");
	return num;
    }
    if (!dk)
    {
	q = m[0] + m[1];
    	dq = dm[0]*dm[0] + dm[1]*dm[1];
    }
    else
    {
	q = m[0];
    	dq = dm[0]*dm[0];
    }
    for (i = 2 - dk; i < num; i++)
    {
	q -= m[i];
	dq += dm[i]*dm[i];
    }
    dq = sqrt(dq);
    
    if (md == 0 || md == 3)
    {
    	if (md == 3) printf("\tfor the reaction: ");
    	else printf("Q-value = %f +/- %f MeV\n\tfor the %s: ",
		q/1000,dq/1000,(dk ? "decay" : "reaction"));
        
    	if (!dk) 
	{
	    if (dec[1] != 1) printf("%d%c%c(",a[1],el[1][0],el[1][1]);
	    else printf(" %c%c(",el[1][0],el[1][1]);
	    if (dec[0] != 1) printf("%d%c%c,", a[0],el[0][0],el[0][1]);
	    else printf(" %c%c,", el[0][0],el[0][1]);
	}
	else printf("(%d%c%c,",
	    a[0],el[0][0],el[0][1]);
        
    	for (i = 2 - dk; i < (num - 1); i++)
	{
	    if (dec[i] != 1) printf("%c%d%c%c",i > (2 -dk) ? ' ' : (char)0,
                                a[i],el[i][0],el[i][1]);
	    else printf(" %c%c",el[i][0],el[i][1]);
	}
    
    	printf(")%d%c%c\n",a[num-1],el[num-1][0],el[num-1][1]);
    }
    return num;
} /*END qval()*/

/*==========================================================================*/
/* read_masses: read masses, z, a etc. from mass-table file	    	    */
/****************************************************************************/
int read_masses(int p, int fmt)
{
    int     c = 0, cnt = 0, i = 0, ln = 0, len = 1, ta = 0, tn = 0, tz = 0;
    int     rdchrs, dig = 0;
    char    tel[R_MAX] = "", tans1[3] = "";
    fpos_t  pos1;
	    
    /*get length of decay particle name*/
    if (isalpha((int)el[p][1]) != 0 ||
    	    ! strncmp(el[p]+1, "-", 1) ||
    	    ! strncmp(el[p]+1, "+", 1)) len++;
    
    /*check list of decay particles for a match*/
    for (i = 0; (i < DMAX && fmt != 0 && a[p] == 0); i++)
    {
/*	printf("el[p] = %c%c, decel[%d] = %c%c strlen = %d\n",
		el[p][0],el[p][1],i,decel[i][0],decel[i][1],len);*/
	if (! strncmp(el[p],decel[i],len) )
	{
	    dec[p] = 1;
	    z[p]=decz[i];
	    a[p]=deca[i];
	    m[p]=decm[i];
/*	    printf("Found decay particle match: %c%c %d\n",
		    decel[i][0],decel[i][1],p);*/
	    return 0;
	}
    }
    /*reset length*/
    len = 1;
    
    /*first few format statements for each line of the mass table*/
    /*a1,i3,i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,*/
    /* ,N-Z, N, Z, A,   el,       m_xs,dm_xs,*/
        
    fgetpos(fmtab, &pos1);
    
    /*get length of element name*/
    if (fmt == 1) if (isalpha((int)el[p][1]) != 0) len++;
    
    /*read m and el*/
    while (1)
    {
	cnt++;
    	/*read the first 4 characters*/
    	for (i = 0; i < 4; i++)
	{
	    c = fgetc(fmtab);
	    memset(&tel[i], ' ', 1);
	}

    	/*read N, Z, A, el*/
    	fscanf(fmtab, "%d %d %d %s", &tn, &tz, &ta, tel);
	if (strlen(tel) == 2) ln = 2;
	else ln = len;
	
/*      if (cnt < 6) printf("In read_masses: tel=%s\n",tel);
	if (cnt < 6) printf("Debug1: fmt=%d, tn=%d tz=%d, ta=%d tel=%c%c len=%d\n",
    	    	fmt,tn,tz,ta,tel[0],tel[1], (int)strlen(tel));
	if (cnt < 6) printf("Debug1: a[p]=%d, z[p]=%d, tel=%s len=%d el[p] = %c%c\n",
    	    	a[p],z[p],tel, (int)strlen(tel), el[p][0],el[p][1]);*/
		
    	if ( (fmt == 0 && tz == z[p] && ta == a[p]) || (fmt == 1 && ta == a[p]
		&& ( ! strncasecmp(el[p],tel,ln)) ) )
    	{
    	    if (fmt == 0) memcpy(el[p], tel, strlen(tel));
	    if (fmt == 1) z[p] = tz;
    	    /*read the next 6 characters*/
    	    for (i = 0; i < 6; i++) c = fgetc(fmtab);
    	    
    	    /*read m and dm*/
    	    fscanf(fmtab, "%f", &m[p]);
	    /*read possible hash (#)*/
	    c = fgetc(fmtab);
    	    fscanf(fmtab, "%f", &dm[p]);
	    /*read possible hash (#)*/
	    c = fgetc(fmtab);
	    
    	    /*read next few parameters*/
    	    fscanf(fmtab, "%f", &bea[p]);
	    /*read possible hash (#)*/
	    c = fgetc(fmtab);
    	    fscanf(fmtab, "%f", &dbea[p]);
	    /*read possible hash (#)*/
	    c = fgetc(fmtab);
    	    fscanf(fmtab, "%s",tans1);
            
/*            printf("In read_masses: dbea=%f, tans1:%s:\n",dbea[p],tans1);*/

    	    /*read the next rdchrs characters: for M_TABYR <=2016 read 21
                                               for M_TABYR >=2020 read 25*/
    	    if (M_TABYR < 2020) rdchrs = 21;
            else rdchrs = 25;
            for (i = 0; i < rdchrs; i++) c = fgetc(fmtab);
	    
    	    /*read amu and damu*/	    
    	    fscanf(fmtab, "%d %f %f", &dig, &ba[p], &dba[p]);
	    ba[p] /= 1000000;
	    ba[p] += dig;
            dba[p] /= 1000000;
    	    break;
    	}
    	else
    	{
    	    while( (c = fgetc(fmtab)) != '\n' && (c != EOF) ) ;
	    
	    if (c == EOF)
	    {
		printf("***No match found in mass table for: "
			"z=%d, el=%c%c a=%d\n",z[p],el[p][0],el[p][1],a[p]);
    	    	fsetpos(fmtab, &pos1);
		return -1;
	    }
    	}
    }
    /*reset file position to pos1*/
    fsetpos(fmtab, &pos1);
    return 0;
} /*END read_masses()*/

/*==========================================================================*/
/* skip_hash: skip comment lines starting with hash (#) at start of file    */
/****************************************************************************/
void skip_hash(FILE *file)
{
    int     i = 0, hash = 0;
    fpos_t  pos;
     	   
    for (i = 0; i < 1000; i++)
    {
    	/*store position corresponding to start of line*/
    	fgetpos(file, &pos);
    	if ( (hash = fgetc(file)) != '#' )
    	{
    	    /*reset position to start of line*/
    	    fsetpos(file, &pos);
    	    break;
    	}
    	/*read rest of line*/
    	while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
    	    ;
    }
} /*END skip_hash()*/

/*==========================================================================*/
/* special: load special notation isotopes into arrays	    	    	    */
/****************************************************************************/
void special(void)
{
    extern  int  sz[SMAX], sa[SMAX];
    extern  char sel[SMAX][2], ssel[SMAX][2];
    int     i;
    
    /*neutron*/
    memcpy(sel[0], "n", 1);
    memcpy(ssel[0], "n", 1);
    sa[0] = 1;
    sz[0] = 0;
    /*proton*/
    memcpy(sel[1], "p", 1);
    memcpy(ssel[1], "H", 1);
    sa[1] = 1;
    sz[1] = 1;
    /*deuteron*/
    memcpy(sel[2], "d", 1);
    memcpy(ssel[2], "H", 1);
    sa[2] = 2;
    sz[2] = 1;
    /*triton*/
    memcpy(sel[3], "t", 1);
    memcpy(ssel[3], "H", 1);
    sa[3] = 3;
    sz[3] = 1;
    /*alpha*/
    memcpy(sel[4], "a", 1);
    memcpy(ssel[4], "He", 2);
    sa[4] = 4;
    sz[4] = 2;

    printf("Allowed special notation for isotopes is: ");
    for (i = 0; i < SMAX; i++)
	printf(" %c%c (%d%c%c)%c",sel[i][0],sel[i][1],sa[i],
	  ssel[i][0],ssel[i][1],(SMAX-i-1 ? ',' : ' '));
    printf("\n");
} /*END special()*/

/*==========================================================================*/
/* store_colours: store colours in clr[][] array                            */
/****************************************************************************/
void store_colours()
{
    /*none*/
    strcpy(clr[0], "\033[0;0m");
    /*white on red background*/
    strcpy(clr[1], "\033[0;41m");
    /*green no background*/
    strcpy(clr[2], "\033[0;32m");
    /*white on blue background*/
    strcpy(clr[3], "\033[0;44m");
    /*bold*/
    strcpy(clr[4], "\033[1m");
    /*underline*/
    strcpy(clr[5], "\033[0;4m");
} /*END store_colours()*/  

/*==========================================================================*/
/* zero_par: zero external parameters	    	    	    	    	    */
/****************************************************************************/
void zero_par(void)
{
    extern int     z[R_MAX],a[R_MAX], dec[DMAX];
    extern float   m[R_MAX], dm[R_MAX];
    extern char    el[R_MAX][2], decel[DMAX][2];
    int     i = 0, j = 0;
    
    for (i = 0; i < DMAX; i++)
    {
	dec[i] = 0;
    }
        
    /*zero external variables*/
    for (i = 0; i < R_MAX; i++)
    {
	dk = 0;
	z[i] = 0;
	a[i] = 0;
	m[i] = 0.0;
	dm[i] = 0.0;
	ba[i] = 0.0;
	dba[i] = 0.0;
	bea[i] = 0.0;
	dbea[i] = 0.0;
    	for (j = 0; j < 2; j++) memset(&el[i][j], ' ', 1);
    }
    q = 0.0; dq = 0.0;
    return ;
} /*END zero_par()*/
