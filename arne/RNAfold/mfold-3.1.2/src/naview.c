/*
*   NAVIEW -- A program to make a modified radial drawing of an RNA
*   secondary structure. 
*   
*   Copyright (c) 1988 Robert E. Bruccoleri
*   Copying of this software, in whole or in part, is permitted
*   provided that the copies are not made for commercial purposes,
*   appropriate credit for the use of the software is given, this
*   copyright notice appears, and notice is given that the copying
*   is by permission of Robert E. Bruccoleri. Any other copying
*   requires specific permission. 
*
*   Usage: RADIAL ct-file PLT2-file
*
*          ct-file is produced by Zuker's RNA folding program.
*          
*   See R. Bruccoleri and G. Heinrich, Computer Applications in the
*   Biosciences, 4, 167-173 (1988) for a full description.
*   
*   Michael Zuker has made the following changes (1989):
*
*   bases[i].hist_num : This is the historical numbering of the ith
*                       base in the original sequence. This number is
*                       picked up as the last item of each record in
*                       the .ct file.
*
*   istart : This is the first position to label on the plot. The
*            historical number for this base is the smallest multiple
*            of label_rate that is >= the historical number of the
*            first base.
*
*   In November 1997, Michael Zuker made a number of changes to bring
*   naview up to modern standards. All functions defined in naview are
*   now declared before main() with arguments and argument types.
*   When functions are defined, their argument types are declared 
*   with the function and these definitions are removed after the '{'.
*   The 'void' declaration was used as necessary for functions.
*
*   The troublesome na_scanf function was deleted and replaced by
*   scanf. Finally, there is now no default for the minimum separation
*   of bases. A floating point number must be entered. However, as
*   before an entry < 0 will be moved up to 0 and an entry > 0.5
*   will be reduced to 0.5.
*/

/* Feb 3, 1998 by Darrin Stewart  The line much below:
   csz=csz/scale; makes the character size depend on the internal
   scaling of naview.  It works best when character size is entered as .3


   Feb 6, 1998 by D.S the .5 above was changed to 1.0
   This was done to keep letters from being written on top
   of each other in loops.
*/

#include <math.h>
#include <stdio.h>

typedef int LOGICAL;
#define logical LOGICAL

#define f77_true -1
#define f77_false 0
#define null 0
#define true 1
#define false 0
#ifdef vms
#   define FATAL_ERROR 4
#   define SUCCESS 1
#else
#   define FATAL_ERROR 1
#   define SUCCESS 0
#endif

#define type_alloc(type) (type *) malloc(sizeof(type))

#define struct_alloc(structure_name) type_alloc(struct structure_name)

#define add_double_list(head,tail,newp) {\
	(newp)->next = (newp)->prev = null; \
        if ((head) == null) (head) = (tail) = (newp); \
	else { \
	     (tail)->next = (newp); \
	     (newp)->prev = (tail); \
	     (tail) = (newp); \
	     } \
	}

static float pi = 3.141592653589793;
static double pid = 3.141592653589793;
static float rad120 = (120.0/180.0*3.141592653589793);
static float dtor = (3.141592653589793/180.0);
static float bignum = 1.7e38;
static float anum = 9999.0;
static float eps1 = 1.0e-6;
static int bigint = 2147483647;


/*
*   Function data type definitions
*/

extern float minf2(),maxf2();
	     
FILE *inf,*outf;

struct base {
    char name;
    int mate,hist_num;
    float x,y;
    logical extracted;
    struct region *region;
    } *bases;
    
struct region {
    int start1,end1,start2,end2;
    };
	
struct loop {
    int nconnection;
    struct connection **connections;
    int number;
    int depth;
    logical mark;
    float x,y,radius;
    };
    
struct connection {
    struct loop *loop;
    struct region *region;
    int start,end;	    /* Start and end form the 1st base pair of the region. */
    float xrad,yrad,angle;
    logical extruded;	    /* True if segment between this connection and
			       the next must be extruded out of the circle */
    logical broken;	    /* True if the extruded segment must be drawn long. */
    };

int nbase,nregion,loop_count;

struct loop *anchor,*root,*loops;

struct region *regions;

extern struct loop *construct_loop();

struct radloop
    {
	float radius;
	int loopnumber;
	struct radloop *next,*prev;
    };

struct radloop *rlphead,*rlptail;

char *dev;

float xsize,ysize;

float lencut;

logical mark_loops,draw_bases;

float csz;

int label_rate;

logical dot_pairs;

int mosaicx,mosaicy;

float glob_rot;

#define LINEMX 500
char line[LINEMX],title[LINEMX];

logical debug = false;

void read_in_bases();
void find_regions();
void dump_loops();
void find_central_loop();
void determine_depths();
void traverse_loop(struct loop *lp,struct connection *anchor_connection);
void determine_radius(struct loop *lp,float lencut);
void generate_region(struct connection *cp);
void construct_extruded_segment(struct connection *cp,struct connection *cpnext);
void find_center_for_arc(int n,float b,float *hp,float *thetap);
void write_plt2_output();
void trim(char *st,int *stlenp);
void chcnbl(char *st,int stlen);


main(argc,argv)
int argc;
char *argv[];

{   char *cp;
    int fd;
    logical change_rad;
    int ilp,linel;
    float r;
    struct radloop *rlp;

    if (argc != 3) {
	printf("Usage is RADIAL ct-file PLT2-file, please try again.\n");
	exit(FATAL_ERROR);
	}
    inf = fopen(argv[1],"r");
    if (inf == null) {
	printf("Unable to open input file %s.\n",argv[1]);
	exit(FATAL_ERROR);
	}
    if ((cp = fgets(line,LINEMX,inf)) == NULL) {
	printf("Unable to read header record.\n");
	exit(FATAL_ERROR);
	}
#   ifdef vms
	if (strlen(line) <= 1) {
	    if ((cp = fgets(line,LINEMX,inf)) == NULL) {
		printf("Unable to read header record.\n");
		exit(FATAL_ERROR);
		}
	    }
#   endif
    line[76] = '\0';
    linel = strlen(line);
    chcnbl(line,linel);
    trim(line,&linel);
    line[linel] = '\0';
    sscanf(line,"%5d",&nbase);
    printf("%d bases specified in file. Title is \n%s\n",nbase,&line[5]);
    strcpy(title,&line[5]);
    outf = fopen(argv[2],"w");
    if (outf == NULL) {
	printf("Unable to open output file %s.\n",argv[2]);
	exit(FATAL_ERROR);
	}
    bases = (struct base *) malloc(sizeof(struct base)*(nbase+1));
    regions = (struct region *) malloc(sizeof(struct region)*(nbase+1));
    read_in_bases();
    printf("Please specify a PLT2 device string:\n");
    fgets(line,LINEMX,stdin);
    linel = strlen(line);
    if (line[linel-1] == '\n') {
	line[--linel] = '\0';
	}
    dev = (char *) malloc(linel+1);
    strcpy(dev,line);
    printf("Please specify the plotter device dimensions in cm.\n");
    printf("Horizontal then vertical sizes separated by a space:\n");
    scanf("%f%f",&xsize,&ysize);
    draw_bases = ask("Should bases be displayed?");
    mark_loops = ask("Should loops numbers be displayed?");
    lencut = 1.0;
    printf("Please specify the minimum permissible separation between bases.\n");
    printf("A size > %f will be reduced to %f.\n",lencut,lencut);
    scanf("%f",&lencut);
    lencut = maxf2(0.0,minf2(1.0,lencut));
    rlphead = rlptail = null;
    change_rad = ask("Do you want to change loop radii?");
    if (change_rad) {
	do {
	    printf("Please specify a loop number and radius. Type zeroes to quit:\n");
	    r = 0.0;
	    scanf("%d%f",&ilp,&r);
	    if (ilp > 0) {
		if (r <= 0.0)
		    printf("Bad radius specified. It must be positive.\n");
		else {
		    rlp = struct_alloc(radloop);
		    rlp->radius = r;
		    rlp->loopnumber = ilp;
		    add_double_list(rlphead,rlptail,rlp);
		    }
		}
	    }
	while    (ilp > 0);
	}
    printf("Please specify the character size:\n");
    scanf("%f",&csz);
    dot_pairs = ask("Do you want base pairs drawn as dots?");
    printf("What frequency for labeling bases (specify 0 for none)?\n");
    scanf("%d",&label_rate);
    if (label_rate < 0) label_rate = 0;
    if (ask("Do you want a mosaic?")) {
	printf("How many display pages/screens for the horizontal dimension?\n");
	scanf("%d",&mosaicx);
	printf("How many display pages/screens for the vertical dimension?\n");
	scanf("%d",&mosaicy);
	if (mosaicx < 1) mosaicx = 1;
	if (mosaicy < 1) mosaicy = 1;
	}
    else {
	mosaicx = 1;
	mosaicy = 1;
	}
    printf("Please specify a rotation for the drawing (positive is clockwise):\n");
    glob_rot = 0.0;
    scanf("%f",&glob_rot);
    printf("Now constructing tree of loops...\n");
    find_regions();
    loop_count = 0;
    loops = (struct loop *) malloc(sizeof(struct loop)*(nbase+1));
    anchor = construct_loop(0);
    find_central_loop();
    if (debug) dump_loops();
    printf("Now constructing the drawing...\n");
    traverse_loop(root,null);
    printf("Now writing PLT2 commands...\n");
    write_plt2_output();
    printf("All done.\n");
    exit(SUCCESS);
    }

void read_in_bases()
/*
*   Read in the bases from the CT-file.
*/

{   char *cp;
    int i;

/*
*   Set up an origin.
*/
    bases[0].name = 'o';
    bases[0].mate = 0;
    bases[0].hist_num = 0;
    bases[0].extracted = false;
    bases[0].x = anum;
    bases[0].y = anum;
    i = 1;
    for (;;) {
	cp = fgets(line,LINEMX,inf);
	if (cp == null) break;
	sscanf(line,"%*d %1s %*d %*d %d %d",&bases[i].name,&bases[i].mate
                                           ,&bases[i].hist_num);
/*
*	The code in this program depends on mate being false (namely zero)
*	where there is no mate.
*/
	bases[i].extracted = false;
	bases[i].x = anum;
	bases[i].y = anum;
	i++;
	}
    if (--i != nbase) {
	printf("Number of bases read from file, %d, doesn't match header record.\n",i);
	nbase = i;
	}
    }

void find_regions()
/*
*   Identifies the regions in the structure.
*/


{   int i,mate,j,nb1;
    logical *mark,found;
    
    nb1 = nbase + 1;
    mark = (int *) malloc(sizeof(int)*nb1);
    for (i = 0; i < nb1; i++) mark[i] = false;
    nregion = 0;
    for (i=0; i<=nbase; i++) {
	if ( (mate = bases[i].mate) && !mark[i]) {
	    regions[nregion].start1 = i;
	    regions[nregion].end2 = mate;
	    mark[i] = true;
	    mark[mate] = true;
	    bases[i].region = bases[mate].region = &regions[nregion];
            for (i++,mate--; 
		 i<mate && bases[i].mate == mate; 
		 i++,mate--) {
		mark[i] = mark[mate] = true;
		bases[i].region = bases[mate].region = &regions[nregion];
		}
	    regions[nregion].end1 = --i;
	    regions[nregion].start2 = mate+1;
	    if (debug) {
		if (nregion == 0) printf("\nRegions are:\n");
		printf("Region %d is %d-%d and %d-%d with gap of %d.\n",
		       nregion+1,regions[nregion].start1,regions[nregion].end1,
		       regions[nregion].start2,regions[nregion].end2,
		       regions[nregion].start2-regions[nregion].end1+1);
		}
	    nregion++;
	    }
	}
    free(mark);
    }

struct loop *construct_loop(ibase)
/*
*   Starting at residue ibase, recursively constructs the loop containing
*   said base and all deeper bases. 
*/

int ibase;

{   int i,mate;
    struct loop *retloop,*lp;
    struct connection *cp,**cpp;
    struct region *rp;
    struct radloop *rlp;
    
    retloop = &loops[loop_count++];
    retloop->nconnection = 0;
    retloop->connections = (struct connection **) malloc(sizeof(struct connection *));
    retloop->depth = 0;
    retloop->number = loop_count;
    retloop->radius = 0.0;
    for (rlp = rlphead;  rlp;  rlp = rlp->next)
	if (rlp->loopnumber == loop_count) retloop->radius = rlp->radius;
    i = ibase;
    do {
	if ((mate = bases[i].mate) != 0) {
	    rp = bases[i].region;
	    if (!bases[rp->start1].extracted) {
                if (i == rp->start1) {
		    bases[rp->start1].extracted = true;
		    bases[rp->end1].extracted = true;
		    bases[rp->start2].extracted = true;
		    bases[rp->end2].extracted = true;
		    lp = construct_loop(rp->end1 < nbase ? rp->end1+1 : 0);	
		    }
		else if (i == rp->start2){
		    bases[rp->start2].extracted = true;
		    bases[rp->end2].extracted = true;
		    bases[rp->start1].extracted = true;
		    bases[rp->end1].extracted = true;
		    lp = construct_loop(rp->end2 < nbase ? rp->end2+1 : 0);
		    }
		else {
		    printf("\n\nError detected in construct_loop. i = %d not found in region table.\n",i);
		    exit(FATAL_ERROR);
		    }
		retloop->connections = (struct connection **) 
				       realloc(retloop->connections,
					       (++retloop->nconnection+1) * 
					       sizeof(struct connection *));
		retloop->connections[retloop->nconnection-1] = cp = 
		    struct_alloc(connection);
		retloop->connections[retloop->nconnection] = null;
		cp->loop = lp;
		cp->region = rp;
		if (i == rp->start1) {
		    cp->start = rp->start1;
		    cp->end = rp->end2;
		    }
		else {
		    cp->start = rp->start2;
		    cp->end = rp->end1;
		    }
		cp->extruded = false;
		cp->broken = false;
		lp->connections = (struct connection **) 
				  realloc(lp->connections,
					  (++lp->nconnection+1) * 
					     sizeof(struct connection *));
		lp->connections[lp->nconnection-1] = cp = 
		    struct_alloc(connection);
		lp->connections[lp->nconnection] = null;
		cp->loop = retloop;
		cp->region = rp;
		if (i == rp->start1) {
		    cp->start = rp->start2;
		    cp->end = rp->end1;
		    }
		else {
		    cp->start = rp->start1;
		    cp->end = rp->end2;
		    }
		cp->extruded = false;
		cp->broken = false;
		}
	    i = mate;
	    }
	if (++i > nbase) i = 0;
        }
    while (i != ibase);
    return retloop;
    }

void dump_loops()
/*
*   Displays all the loops.
*/

{   int il,ilp,irp;
    struct loop *lp;
    struct connection *cp,**cpp;
    
    printf("\nRoot loop is #%d\n",(root-loops)+1);
    for (il=0; il < loop_count; il++) {
	lp = &loops[il];
	printf("Loop %d has %d connections:\n",il+1,lp->nconnection);
	for (cpp = lp->connections; cp = *cpp; cpp++) {
	    ilp = (cp->loop - loops) + 1;
	    irp = (cp->region - regions) + 1;
	    printf("  Loop %d Region %d (%d-%d)\n",
		       ilp,irp,cp->start,cp->end);
	    }
	}
    }

void find_central_loop()
/*
*   Find node of greatest branching that is deepest.
*/

{   struct loop *lp;
    int maxconn,maxdepth,i;

    determine_depths();
    maxconn = 0;
    maxdepth = -1;
    
    for (i=0; i<loop_count; i++) {
	lp = &loops[i];
	if (lp->nconnection > maxconn) {
	    maxdepth = lp->depth;
	    maxconn = lp->nconnection;
	    root = lp;
	    }
	else if (lp->depth > maxdepth && lp->nconnection == maxconn) {
	    maxdepth = lp->depth;
	    root = lp;
	    }
	}
    }

void determine_depths()
/*
*   Determine the depth of all loops.
*/

{   struct loop *lp;
    int i,j;
    
    for (i=0; i<loop_count; i++) {
	lp = &loops[i];
	for (j=0; j<loop_count; j++) loops[j].mark = false;
	lp->depth = depth(lp);
	}
    }

depth(lp)
/*
*   Determines the depth of loop, lp. Depth is defined as the minimum
*   distance to a leaf loop where a leaf loop is one that has only one
*   or no connections.
*/

struct loop *lp;

{   struct connection *cp,**cpp;
    int count,ret,d;
    
    if (lp->nconnection <= 1) return 0;
    if (lp->mark) return -1;
    lp->mark = true;
    count = 0;
    ret = 0;
    for (cpp=lp->connections; cp = *cpp; cpp++) {
	d = depth(cp->loop);
	if (d >= 0) {
	    if (++count == 1) ret = d;
	    else if (ret > d) ret = d;
	    }
	}
    lp->mark = false;
    return ret+1;
    }

void traverse_loop(lp,anchor_connection)
/*
*   This is the workhorse of the display program. The algorithm is
*   recursive based on processing individual loops. Each base pairing
*   region is displayed using the direction given by the circle diagram,
*   and the connections between the regions is drawn by equally spaced
*   points. The radius of the loop is set to minimize the square error
*   for lengths between sequential bases in the loops. The "correct"
*   length for base links is 1. If the least squares fitting of the
*   radius results in loops being less than 1/2 unit apart, then that
*   segment is extruded. 
*   
*   The variable, anchor_connection, gives the connection to the loop
*   processed in an previous level of recursion.
*/

struct loop *lp;
struct connection *anchor_connection;

{   float sumn,sumd,xs,ys,xe,ye,xn,yn,angle,angleinc,r,dt;
    float radius,xc,yc,xo,yo,astart,aend,a,mindit,dit;
    struct connection *cp,*cpnext,**cpp,*acp,*cpprev;
    struct region *rp;
    int end,start,i,j,l,n,mate,ic,iprev,i2,nl,ns2,nstem,imindit;
    float astartnext,da,xlr,ylr,xla,yla,angle_lp,radius1,ci,maxang;
    int count,icstart,icend,icmiddle,icroot;
    logical done,done_all_connections,rooted,found;
    int iend,sign;
    float midx,midy,nrx,nry,mx,my,vx,vy,dotmv,nmidx,nmidy,dotcn;
    int icstart1,icup,icdown,icnext,il1,il2,il3,il4,direction;
    float xlp,ylp,rlp,dan,da1,dx,dy,rr,rr1,rr2,xcross,ycross,rcross,xout,yout;
    float cpx,cpy,cpnextx,cpnexty,cnx,cny,rcn,rc,lnx,lny,rl,ac,acn,sx,sy,dcp;
    int maxloop,imaxloop;
    
    angleinc = 2 * pi / (nbase+1);
    acp = null;
    icroot = -1;
    for (cpp=lp->connections, ic = 0; 
         cp = *cpp; 
         cpp++, ic++) {
/*	xs = cos(angleinc*cp->start);
	ys = sin(angleinc*cp->start);
	xe = cos(angleinc*cp->end);
	ye = sin(angleinc*cp->end); */
	xs = -sin(angleinc*cp->start);
	ys = cos(angleinc*cp->start);
	xe = -sin(angleinc*cp->end);
	ye = cos(angleinc*cp->end);
	xn = ye-ys;
	yn = xs-xe;
	r = sqrt(xn*xn + yn*yn);
	cp->xrad = xn/r;
	cp->yrad = yn/r;
	cp->angle = atan2(yn,xn);
	if (cp->angle < 0.0) cp->angle += 2*pi;
	if (anchor_connection != null &&
	    anchor_connection->region == cp->region) {
	    acp = cp;
	    icroot = ic;
	    }
        }
	
set_radius:
    determine_radius(lp,lencut);
    radius = lp->radius;
    if (anchor_connection == null) xc = yc = 0.0;
    else {
	xo = (bases[acp->start].x+bases[acp->end].x) / 2.0;
	yo = (bases[acp->start].y+bases[acp->end].y) / 2.0;
	xc = xo - radius * acp->xrad;
	yc = yo - radius * acp->yrad;
	}
	
/*
*   The construction of the connectors will proceed in blocks of
*   connected connectors, where a connected connector pairs means
*   two connectors that are forced out of the drawn circle because they
*   are too close together in angle.
*/

/*
*   First, find the start of a block of connected connectors
*/

    if (icroot == -1) 
	icstart = 0;
    else icstart = icroot;
    cp = lp->connections[icstart];
    count = 0;
    if (debug) printf("Now processing loop %d\n",lp->number);
    done = false;
    do {
	j = icstart - 1;
	if (j < 0) j = lp->nconnection - 1;
	cpprev = lp->connections[j];
	if (!connected_connection(cpprev,cp)) {
	    done = true;
	    }
	else {
	    icstart = j;
	    cp = cpprev;
	    }
	if (++count > lp->nconnection) {
/*
*	    Here everything is connected. Break on maximum angular separation
*	    between connections. 
*/
	    maxang = -1.0;
	    for (ic = 0;  ic < lp->nconnection;  ic++) {
		j = ic + 1;
		if (j >= lp->nconnection) j = 0;
		cp = lp->connections[ic];
		cpnext = lp->connections[j];
		ac = cpnext->angle - cp->angle;
		if (ac < 0.0) ac += 2*pi;
		if (ac > maxang) {
		    maxang = ac;
		    imaxloop = ic;
		    }
		}
	    icend = imaxloop;
	    icstart = imaxloop + 1;
	    if (icstart >= lp->nconnection) icstart = 0;
	    cp = lp->connections[icend];
	    cp->broken = true;
	    done = true;
	    }
	} while    (!done);
    done_all_connections = false;
    icstart1 = icstart;
    if (debug) printf("Icstart1 = %d\n",icstart1);
    while (!done_all_connections) {
	count = 0;
	done = false;
	icend = icstart;
	rooted = false;
	while (!done) {
	    cp = lp->connections[icend];
	    if (icend == icroot) rooted = true;
	    j = icend + 1;
	    if (j >= lp->nconnection) {
		j = 0;
		}
	    cpnext = lp->connections[j];
	    if (connected_connection(cp,cpnext)) {
		if (++count >= lp->nconnection) break;
		icend = j;
		}
	    else {
		done = true;
		}
	    }
	icmiddle = find_ic_middle(icstart,icend,anchor_connection,acp,lp);
	ic = icup = icdown = icmiddle;
	if (debug)
	    printf("IC start = %d  middle = %d  end = %d\n",
		   icstart,icmiddle,icend);
        done = false;
	direction = 0;
	while (!done) {
	    if (direction < 0) {
		ic = icup;
		}
	    else if (direction == 0) {
		ic = icmiddle;
		}
	    else {
		ic = icdown;
		}
	    if (ic >= 0) {
		cp = lp->connections[ic];
		if (anchor_connection == null || acp != cp) {
		    if (direction == 0) {
			astart = cp->angle - asin(1.0/2.0/radius);
			aend = cp->angle + asin(1.0/2.0/radius);
			bases[cp->start].x = xc + radius*cos(astart);
			bases[cp->start].y = yc + radius*sin(astart);
			bases[cp->end].x = xc + radius*cos(aend);
			bases[cp->end].y = yc + radius*sin(aend);
			}
		    else if (direction < 0) {
			j = ic + 1;
			if (j >= lp->nconnection) j = 0;
			cp = lp->connections[ic];
			cpnext = lp->connections[j];
			cpx = cp->xrad;
			cpy = cp->yrad;
			ac = (cp->angle + cpnext->angle) / 2.0;
			if (cp->angle > cpnext->angle) ac -= pi;
			cnx = cos(ac);
			cny = sin(ac);
			lnx = cny;
			lny = -cnx;
			da = cpnext->angle - cp->angle;
			if (da < 0.0) da += 2*pi;
			if (cp->extruded) {
			    if (da <= pi/2) rl = 2.0;
			    else rl = 1.5;
			    }
			else {
			    rl = 1.0;
			    }
			bases[cp->end].x = bases[cpnext->start].x + rl*lnx;
			bases[cp->end].y = bases[cpnext->start].y + rl*lny;
			bases[cp->start].x = bases[cp->end].x + cpy;
			bases[cp->start].y = bases[cp->end].y - cpx;
			}
		    else {
			j = ic - 1;
			if (j < 0) j = lp->nconnection - 1;
			cp = lp->connections[j];
			cpnext = lp->connections[ic];
			cpnextx = cpnext->xrad;
			cpnexty = cpnext->yrad;
			ac = (cp->angle + cpnext->angle) / 2.0;
			if (cp->angle > cpnext->angle) ac -= pi;
			cnx = cos(ac);
			cny = sin(ac);
			lnx = -cny;
			lny = cnx;
			da = cpnext->angle - cp->angle;
			if (da < 0.0) da += 2*pi;
			if (cp->extruded) {
			    if (da <= pi/2) rl = 2.0;
			    else rl = 1.5;
			    }
			else {
			    rl = 1.0;
			    }
			bases[cpnext->start].x = bases[cp->end].x + rl*lnx;
			bases[cpnext->start].y = bases[cp->end].y + rl*lny;
			bases[cpnext->end].x = bases[cpnext->start].x - cpnexty;
			bases[cpnext->end].y = bases[cpnext->start].y + cpnextx;
			}
		    }
		}
	    if (direction < 0) {
		if (icdown == icend) {
		    icdown = -1;
		    }
		else if (icdown >= 0) {
		    if (++icdown >= lp->nconnection) {
			icdown = 0;
			}
		    }
		direction = 1;
		}
	    else {
		if (icup == icstart) icup = -1;
		else if (icup >= 0) {
		    if (--icup < 0) {
			icup = lp->nconnection - 1;
			}
		    }
		direction = -1;
		}
	    done = icup == -1 && icdown == -1;
            }
	icnext = icend + 1;
	if (icnext >= lp->nconnection) icnext = 0;
	if (icend != icstart && (! (icstart == icstart1 && icnext == icstart1))) {
/*
*	    Move the bases just constructed (or the radius) so
*	    that the bisector of the end points is radius distance
*	    away from the loop center.
*/
	    cp = lp->connections[icstart];
	    cpnext = lp->connections[icend];
	    dx = bases[cpnext->end].x - bases[cp->start].x;
	    dy = bases[cpnext->end].y - bases[cp->start].y;
	    midx = bases[cp->start].x + dx/2.0;
	    midy = bases[cp->start].y + dy/2.0;
	    rr = sqrt(dx*dx + dy*dy);
	    mx = dx / rr;
	    my = dy / rr;
	    vx = xc - midx;
	    vy = yc - midy;
	    rr = sqrt(dx*dx + dy*dy);
	    vx /= rr;
	    vy /= rr;
	    dotmv = vx*mx + vy*my;
	    nrx = dotmv*mx - vx;
	    nry = dotmv*my - vy;
	    rr = sqrt(nrx*nrx + nry*nry);
	    nrx /= rr;
	    nry /= rr;
/*
*	    Determine which side of the bisector the center should be.
*/
	    dx = bases[cp->start].x - xc;
	    dy = bases[cp->start].y - yc;
	    ac = atan2(dy,dx);
	    if (ac < 0.0) ac += 2*pi;
	    dx = bases[cpnext->end].x - xc;
	    dy = bases[cpnext->end].y - yc;
	    acn = atan2(dy,dx);
	    if (acn < 0.0) acn += 2*pi;
	    if (acn < ac) acn += 2*pi;
	    if (acn - ac > pi) sign = -1;
	    else sign = 1;
	    nmidx = xc + sign*radius*nrx;
	    nmidy = yc + sign*radius*nry;
	    if (rooted) {
		xc -= nmidx - midx;
		yc -= nmidy - midy;
		}
	    else {
		for (ic=icstart; ; ++ic >= lp->nconnection ? (ic = 0) : 0) {
		    cp = lp->connections[ic];
		    i = cp->start;
		    bases[i].x += nmidx - midx;
		    bases[i].y += nmidy - midy;
		    i = cp->end;
		    bases[i].x += nmidx - midx;
		    bases[i].y += nmidy - midy;
		    if (ic == icend) break;
		    }
		}
	    }
	icstart = icnext;
	done_all_connections = icstart == icstart1;
	}
    for (ic=0; ic < lp->nconnection; ic++) {
	cp = lp->connections[ic];
	j = ic + 1;
	if (j >= lp->nconnection) j = 0;
	cpnext = lp->connections[j];
	dx = bases[cp->end].x - xc;
	dy = bases[cp->end].y - yc;
	rc = sqrt(dx*dx + dy*dy);
	ac = atan2(dy,dx);
	if (ac < 0.0) ac += 2*pi;
	dx = bases[cpnext->start].x - xc;
	dy = bases[cpnext->start].y - yc;
	rcn = sqrt(dx*dx + dy*dy);
	acn = atan2(dy,dx);
	if (acn < 0.0) acn += 2*pi;
	if (acn < ac) acn += 2*pi;
	dan = acn - ac;
	dcp = cpnext->angle - cp->angle;
	if (dcp <= 0.0) dcp += 2*pi;
	if (fabs(dan-dcp) > pi) {
	    if (cp->extruded) {
		printf("Warning from traverse_loop. Loop %d has crossed regions\n",
		       lp->number);
		}
	    else {
		cp->extruded = true;
		goto set_radius;	    /* Forever shamed */
		}
	    }
	if (cp->extruded) {
	    construct_extruded_segment(cp,cpnext);
	    }
	else {
	    n = cpnext->start - cp->end;
	    if (n < 0) n += nbase + 1;
	    angleinc = dan / n;
	    for (j = 1;  j < n;  j++) {
		i = cp->end + j;
		if (i > nbase) i -= nbase + 1;
		a = ac + j*angleinc;
		rr = rc + (rcn-rc)*(a-ac)/dan;
		bases[i].x = xc + rr*cos(a);
		bases[i].y = yc + rr*sin(a);
		}
	    }
	}
    for (ic=0; ic < lp->nconnection; ic++) {
	if (icroot != ic) {
	    cp = lp->connections[ic];
	    generate_region(cp);
	    traverse_loop(cp->loop,cp);
	    }
	}
    n = 0;
    sx = 0.0;
    sy = 0.0;
    for (ic = 0;  ic < lp->nconnection;  ic++) {
	j = ic + 1;
	if (j >= lp->nconnection) j = 0;
	cp = lp->connections[ic];
	cpnext = lp->connections[j];
	n += 2;
	sx += bases[cp->start].x + bases[cp->end].x;
	sy += bases[cp->start].y + bases[cp->end].y;
	if (!cp->extruded) {
	    for (j = cp->end + 1;  j != cpnext->start;  j++) {
		if (j > nbase) j -= nbase + 1;
		n++;
		sx += bases[j].x;
		sy += bases[j].y;
		}
	    }
	}
    lp->x = sx / n;
    lp->y = sy / n;
    }

void determine_radius(struct loop *lp,float lencut)
/*
*   For the loop pointed to by lp, determine the radius of
*   the loop that will ensure that each base around the loop will
*   have a separation of at least lencut around the circle.
*   If a segment joining two connectors will not support this separation,
*   then the flag, extruded, will be set in the first of these
*   two indicators. The radius is set in lp.
*   
*   The radius is selected by a least squares procedure where the sum of the
*   squares of the deviations of length from the ideal value of 1 is used
*   as the error function.
*/

{   float mindit,ci,dt,sumn,sumd,radius,dit,da;
    int count,i,j,end,start,imindit;
    struct connection *cp,*cpnext;
    static float rt2_2 = 0.7071068;

    count = 0;
    do {
	mindit = 1.0e10;
	for (sumd=0.0, sumn=0.0, i=0;
	     i < lp->nconnection;
	     i++) {
	    cp = lp->connections[i];
	    j = i + 1;
	    if (j >= lp->nconnection) j = 0;
	    cpnext = lp->connections[j];
	    end =  cp->end;
	    start = cpnext->start;
	    if (start < end) start += nbase + 1;
	    dt = cpnext->angle - cp->angle;
	    if (dt <= 0.0) dt += 2*pi;
	    if (!cp->extruded) 
		ci = start - end;
	    else {
		if (dt <= pi/2) ci = 2.0;
		else ci = 1.5;
		}
	    sumn += dt * (1.0/ci + 1.0);
	    sumd += dt * dt / ci;
	    dit = dt/ci;
	    if (dit < mindit && !cp->extruded && ci > 1.0) {
		mindit = dit;
		imindit = i;
		}
	    }
	radius = sumn/sumd;
	if (radius < rt2_2) radius = rt2_2;
	if (mindit*radius < lencut) {
	    lp->connections[imindit]->extruded = true;
	    }
	} while (mindit*radius < lencut);
    if (lp->radius > 0.0)
	radius = lp->radius;
    else lp->radius = radius;
    }

logical    connected_connection(cp,cpnext)
/*
*   Determines if the connections cp and cpnext are connected
*/

struct connection *cp,*cpnext;

{

    if (cp->extruded) {
	return true;
	}
    else if (cp->end+1 == cpnext->start) {
	return true;
	}
    else {
	return false;
	}
    }

int    find_ic_middle(icstart,icend,anchor_connection,acp,lp)
/*
*   Finds the middle of a set of connected connectors. This is normally
*   the middle connection in the sequence except if one of the connections
*   is the anchor, in which case that connection will be used.
*/

int icstart,icend;
struct connection *anchor_connection,*acp;
struct loop *lp;

{   int count,ret,ic,i;
    logical done;

    count = 0;
    ret = -1;
    ic = icstart;
    done = false;
    while (!done) {
	if (count++ > lp->nconnection * 2) {
	    printf("Infinite loop detected in find_ic_middle\n");
	    exit(FATAL_ERROR);
	    }
	if (anchor_connection != null && lp->connections[ic] == acp) {
	    ret = ic;
	    }
	done = ic == icend;
	if (++ic >= lp->nconnection) {
	    ic = 0;
	    }
	}
    if (ret == -1) {
	for (i=1, ic=icstart; i<(count+1)/2; i++) {
	    if (++ic >= lp->nconnection) ic = 0;
	    }
	ret = ic;
	}
    return ret;
    }

void generate_region(cp)
/*
*   Generates the coordinates for the base pairing region of a connection
*   given the position of the starting base pair.
*/

struct connection *cp;

{   int l,start,end,i,mate;
    struct region *rp;

    rp = cp->region;
    l = 0;
    if (cp->start == rp->start1) {
	start = rp->start1;
	end = rp->end1;
	}
    else {
	start = rp->start2;
	end = rp->end2;
	}
    if (bases[cp->start].x > anum - 100.0 ||
        bases[cp->end].x > anum - 100.0) {
	printf("Bad region passed to generate_region. Coordinates not defined.\n");
	exit(FATAL_ERROR);
	}
    for (i=start+1; i<=end; i++) {
	l++;
	bases[i].x = bases[cp->start].x + l*cp->xrad;
	bases[i].y = bases[cp->start].y + l*cp->yrad;
	mate = bases[i].mate;
	bases[mate].x = bases[cp->end].x + l*cp->xrad;
	bases[mate].y = bases[cp->end].y + l*cp->yrad;
	}
    }

void construct_circle_segment(start,end)
/*
*   Draws the segment of residue between the bases numbered start
*   through end, where start and end are presumed to be part of a base
*   pairing region. They are drawn as a circle which has a chord given
*   by the ends of two base pairing regions defined by the connections.
*/

int start,end;

{   float dx,dy,rr,h,angleinc,midx,midy,xn,yn,nrx,nry,mx,my,a;
    int l,j,i;

    dx = bases[end].x - bases[start].x;
    dy = bases[end].y - bases[start].y;
    rr = sqrt(dx*dx + dy*dy);
    l = end - start;
    if (l < 0) l += nbase + 1;
    if (rr >= l) {
	dx /= rr;
	dy /= rr;
	for (j = 1;  j < l;  j++) {
	    i = start + j;
	    if (i > nbase) i -= nbase + 1;
	    bases[i].x = bases[start].x + dx*(float)j/(float)l;
	    bases[i].y = bases[start].y + dy*(float)j/(float)l;
	    }
	}
    else {
	find_center_for_arc(l-1,rr,&h,&angleinc);
	dx /= rr;
	dy /= rr;
	midx = bases[start].x + dx*rr/2.0;
	midy = bases[start].y + dy*rr/2.0;
	xn = dy;
	yn = -dx;
	nrx = midx + h*xn;
	nry = midy + h*yn;
	mx = bases[start].x - nrx;
	my = bases[start].y - nry;
	rr = sqrt(mx*mx + my*my);
	a = atan2(my,mx);
	for (j = 1;  j < l;  j++) {
	    i = start + j;
	    if (i > nbase) i -= nbase + 1;
	    bases[i].x = nrx + rr*cos(a+j*angleinc);
	    bases[i].y = nry + rr*sin(a+j*angleinc);
	    }
	}
    }

void construct_extruded_segment(cp,cpnext)
/*
*   Constructs the segment between cp and cpnext as a circle if possible.
*   However, if the segment is too large, the lines are drawn between
*   the two connecting regions, and bases are placed there until the
*   connecting circle will fit.
*/

struct connection *cp,*cpnext;

{   float astart,aend1,aend2,aave,dx,dy,a1,a2,ac,rr,da,dac;
    int start,end,n,nstart,nend;
    logical collision;

    astart = cp->angle;
    aend2 = aend1 = cpnext->angle;
    if (aend2 < astart) aend2 += 2*pi;
    aave = (astart + aend2) / 2.0;
    start = cp->end;
    end = cpnext->start;
    n = end - start;
    if (n < 0) n += nbase + 1;
    da = cpnext->angle - cp->angle;
    if (da < 0.0) {
	da += 2*pi;
	}
    if (n == 2) construct_circle_segment(start,end);
    else {
	dx = bases[end].x - bases[start].x;
	dy = bases[end].y - bases[start].y;
	rr = sqrt(dx*dx + dy*dy);
	dx /= rr;
	dy /= rr;
	if (rr >= 1.5 && da <= pi/2) {
	    nstart = start + 1;
	    if (nstart > nbase) nstart -= nbase + 1;
	    nend = end - 1;
	    if (nend < 0) nend += nbase + 1;
	    bases[nstart].x = bases[start].x + 0.5*dx;
	    bases[nstart].y = bases[start].y + 0.5*dy;
	    bases[nend].x = bases[end].x - 0.5*dx;
	    bases[nend].y = bases[end].y - 0.5*dy;
	    start = nstart;
	    end = nend;
	    }
	do {
	    collision = false;
	    construct_circle_segment(start,end);
	    nstart = start + 1;
	    if (nstart > nbase) nstart -= nbase + 1;
	    dx = bases[nstart].x - bases[start].x;
	    dy = bases[nstart].y - bases[start].y;
	    a1 = atan2(dy,dx);
	    if (a1 < 0.0) a1 += 2*pi;
	    dac = a1 - astart;
	    if (dac < 0.0) dac += 2*pi;
	    if (dac > pi) collision = true;
	    nend = end - 1;
	    if (nend < 0) nend += nbase + 1;
	    dx = bases[nend].x - bases[end].x;
	    dy = bases[nend].y - bases[end].y;
	    a2 = atan2(dy,dx);
	    if (a2 < 0.0) a2 += 2*pi;
	    dac = aend1 - a2;
	    if (dac < 0.0) dac += 2*pi;
	    if (dac > pi) collision = true;
	    if (collision) {
		ac = minf2(aave,astart + 0.5);
		bases[nstart].x = bases[start].x + cos(ac);
		bases[nstart].y = bases[start].y + sin(ac);
		start = nstart;
		ac = maxf2(aave,aend2 - 0.5);
		bases[nend].x = bases[end].x + cos(ac);
		bases[nend].y = bases[end].y + sin(ac);
		end = nend;
		n -= 2;
		}
	    } while    (collision && n > 1);
	}
    }

void find_center_for_arc(int n,float b,float *hp,float *thetap)
/*
*   Given n points to be placed equidistantly and equiangularly on a
*   polygon which has a chord of length, b, find the distance, h, from the
*   midpoint of the chord for the center of polygon. Positive values
*   mean the center is within the polygon and the chord, whereas
*   negative values mean the center is outside the chord. Also, the
*   radial angle for each polygon side is returned in theta.
*    
*   The procedure uses a bisection algorithm to find the correct
*   value for the center. Two equations are solved, the angles
*   around the center must add to 2*pi, and the sides of the polygon
*   excluding the chord must have a length of 1.
*/

{   float h,hhi,hlow,r,disc,theta,e,phi;
    int iter;
#define maxiter 500

    hhi = (n+1) / pi;
    hlow = - hhi - b/(n+1-b);
    iter = 0;
    do {
	h = (hhi + hlow) / 2.0;
	r = sqrt(h*h + b*b/4.0);
	disc = 1.0 - 1.0/2.0/(r*r);
	if (fabs(disc) > 1.0) {
	    printf("Unexpected large magnitude discriminant = %g\n",disc);
	    exit(FATAL_ERROR);
	    }
	theta = acos(disc);
	phi = acos(h/r);
	e = theta * (n+1) + 2*phi - 2*pi;
	if (e > 0.0) {
	    hlow = h;
	    }
	else {
	    hhi = h;
	    }
	} while    (fabs(e) > 0.0001 && ++iter < maxiter);
    if (iter >= maxiter) {
	printf("Iteration failed in find_center_for_arc\n");
	h = 0.0;
	theta = 0.0;
	}
    *hp = h;
    *thetap = theta;
    }
#undef maxiter

void write_plt2_output()
/*
*   Writes the coordinates as a PLT2 command stream.
*/

{   int i,mate,imx,imy,pagecnt,istart;
    float xmin,xmax,ymin,ymax,scalex,scaley,scale,scalecsz;
    float xmn,xmx,ymn,ymx;
    float xn,yn,r,x1,y1,x2,y2,xs,ys,xoff,yoff,ct,st,xr,yr;
    struct loop *lp;

#define okx(x) (xmn-scale <= (x) && (x) <= xmx+scale)
#define oky(y) (ymn-scale <= (y) && (y) <= ymx+scale)

    fprintf(outf,"DEV %s\n",dev);
    fprintf(outf,"CSZ %10.4f\n",csz);
    ct = cos(glob_rot*pi/180.0);
    st = sin(glob_rot*pi/180.0);
    for (i=0; i<=nbase; i++) {
	xr = bases[i].x*ct + bases[i].y*st;
	yr = bases[i].y*ct - bases[i].x*st;
	bases[i].x = xr;
	bases[i].y = yr;
	}
    for (i=0;  i<loop_count;  i++)
	{
	    lp = loops + i;
	    xr = lp->x*ct + lp->y*st;
	    yr = lp->y*ct - lp->x*st;
	    lp->x = xr;
	    lp->y = yr;
	}
    xmax = ymax = -bignum;
    xmin = ymin = bignum;
    for (i=0; i<=nbase; i++) {
	if (bases[i].x > anum - 100.0) {
	    printf("\nError in write_plt2_output -- base %d position is undefined.\n",i);
	    }
	else {
	    xmax = maxf2(xmax,bases[i].x);
	    ymax = maxf2(ymax,bases[i].y);
	    xmin = minf2(xmin,bases[i].x);
	    ymin = minf2(ymin,bases[i].y);
	    }
	}
    scalex = (xmax-xmin) / xsize / (float) mosaicx;
    scaley = (ymax-ymin) / ysize / (float) mosaicy;
    scale = maxf2(scalex,scaley) * 1.02;
    /* The line below was added by Darrin Stewart $$$$$$$$*/
    csz=csz/scale;
    /* It works best when csz is .3, but numbers from .4 to .1 are reasonable*/
    /* The character size is now depended on the internal scaling of the image*/
    /* This seems to be quite accurate. */
    xs = xsize * scale;
    ys = ysize * scale;
    xoff = (xs*mosaicx - xmax + xmin) / 2.0;
    yoff = (ys*mosaicy - ymax + ymin) / 2.0;
    fprintf(outf,"SA %g\n",scale);
    fprintf(outf,"ORI 0.0 0.0\n");
    pagecnt = 0;
    for (imx=1; imx <= mosaicx; imx++) {
	for (imy=1; imy <= mosaicy; imy++) {
	    if (++pagecnt > 1) {
		fprintf(outf,"DUMP\n");
		}
	    if (title[0] != '\0') {
		if (imx == (mosaicx+1)/2 && imy == 1) {
		    fprintf(outf,"CSZ 0.4\n");
		    fprintf(outf,"CTA %10.3f 1.0 \"%s\"\n",xsize/2.0,title);
		    fprintf(outf,"CSZ %10.4f\n",csz);
		    }
		}
	    xmn = xmin + (imx-1)*xs - xoff;
	    ymn = ymin + (imy-1)*ys - yoff;
	    xmx = xmin + imx*xs - xoff;
	    ymx = ymin + imy*ys - yoff;
	    fprintf(outf,"OD %10.3f %10.3f\n",xmn,ymn);
	    fprintf(outf,"COLOR WHITE\n");
	    if (bases[0].x != anum)
		fprintf(outf,"BRI 2\nMOV %10.3f %10.3f\nCIA 0.05\nCIA 0.1\nCIA 0.15\n",
			     bases[0].x,-bases[0].y);
	    fprintf(outf,"BRI 3\n");
	    if (draw_bases) {
		fprintf(outf,"CM BASES %d\n",nbase);   /* Zuker adds nbase */
		for (i=1; i <= nbase; i++) {
		    if (bases[i].x != anum)
			if (okx(bases[i].x) && oky(bases[i].y))
			    fprintf(outf,"CTX %10.3f %10.3f \"%c\"\n",bases[i].x,-bases[i].y,bases[i].name);
		    }
		}
	    fprintf(outf,"CM SEQUENCE LINES\n");
	    if (draw_bases) scalecsz = scale * csz * 1.8;
	    else scalecsz = 0.0;
	    for (i=0; i <= nbase-1; i++) {
		if (bases[i].x != anum && bases[i+1].x != anum) {
		    if (okx(bases[i].x) && oky(bases[i].y) ||
		        okx(bases[i+1].x) && oky(bases[i+1].y)) {
			xn = bases[i+1].x - bases[i].x;
			yn = bases[i+1].y - bases[i].y;
			r = sqrt(xn*xn + yn*yn);
			if (r > scalecsz) {
			    xn /= r;
			    yn /= r;
			    x1 = bases[i].x + xn*scalecsz/2.0;
			    y1 = bases[i].y + yn*scalecsz/2.0;
			    x2 = bases[i+1].x - xn*scalecsz/2.0;
			    y2 = bases[i+1].y - yn*scalecsz/2.0;
			    fprintf(outf,"LI %10.3f %10.3f 0.0 %10.3f %10.3f 0.0\n",
				    x1,-y1,x2,-y2);
			    }
			}
		    }
		}
	   fprintf(outf,"CM BASE PAIRING LINES WITH BASE PAIRS\n");
	   fprintf(outf,"COLOR RED\nBRI 5\n");

#define draw_1bp \
	    if (dot_pairs)\
		{\
		    x1 = (bases[i].x + bases[mate].x) / 2.0;\
		    y1 = (bases[i].y + bases[mate].y) / 2.0;\
		    if (okx(x1) && oky(y1)) \
			fprintf(outf,"MOV %10.3f %10.3f %5d %5d\nCIA 0.0\n",x1,-y1,i,mate);\
		}\
	    else\
		{\
		    xn = bases[mate].x - bases[i].x;\
		    yn = bases[mate].y - bases[i].y;\
		    r = sqrt(xn*xn + yn*yn);\
		    if (r > scalecsz) {\
			xn /= r;\
			yn /= r;\
			x1 = bases[i].x + xn*scalecsz/2.0;\
			y1 = bases[i].y + yn*scalecsz/2.0;\
			x2 = bases[mate].x - xn*scalecsz/2.0;\
			y2 = bases[mate].y - yn*scalecsz/2.0;\
			if (okx(x1) && oky(y1) || \
			    okx(x2) && oky(y2)) \
			    fprintf(outf,"LI %10.3f %10.3f 0.0 %10.3f %10.3f 0.0 %5d %5d\n",\
				    x1,-y1,x2,-y2,i,mate);\
			}\
		}

	   for (i=0; i <= nbase; i++) {
		if ((mate = bases[i].mate) && i < mate ) {
		    if (bases[i].x != anum && bases[mate].x != anum) {
			if (bases[i].name == 'G' && bases[mate].name == 'C' ||
			    bases[i].name == 'C' && bases[mate].name == 'G')
			    {
				draw_1bp;
			    }
			}
		    }
		}

	   fprintf(outf,"COLOR MAGENTA\nBRI 1\n");
	   for (i=0; i <= nbase; i++) {
		if ((mate = bases[i].mate) && i < mate ) {
		    if (bases[i].x != anum && bases[mate].x != anum) {
			if (!(bases[i].name == 'G' && bases[mate].name == 'C' ||
			      bases[i].name == 'C' && bases[mate].name == 'G'))
			    {
				draw_1bp;
			    }
			}
		    }
		}

	    fprintf(outf,"COLOR WHITE\nBRI 3\n");
	    if (label_rate > 0) {
		fprintf(outf,"CM LABELS %d\nCSZ %10.4f\n",nbase,csz);
                istart = label_rate*(1 + bases[1].hist_num/label_rate) 
                            - bases[1].hist_num + 1;
		for (i=istart; i <= nbase; i += label_rate) {
		    if (bases[i].x != anum) {
			if (okx(bases[i].x) && oky(bases[i].y)) {
			    float dx,dy,angle;
			    int ia;
			    if (i == nbase) {
				dx = bases[i].x - bases[i-1].x;
				dy = bases[i].y - bases[i-1].y;
				}
			    else {
				dx = bases[i+1].x - bases[i].x;
				dy = bases[i+1].y - bases[i].y;
				}
			    angle = atan2(dy,dx)/dtor - 90.0;
			    ia = angle;
			    fprintf(outf,"TEX %d %10.3f %10.3f \"  %d\"\n"
                             ,ia,bases[i].x,-bases[i].y,bases[i].hist_num);
			    }
			}
		    }
		}
	    if (mark_loops) {
		fprintf(outf,"CM LOOP LABELS\n");
		for (i=0;  i<loop_count;  i++)
		    {
			lp = loops + i;
			if (okx(lp->x) && oky(lp->y)) 
			    fprintf(outf,"CTX %10.3f %10.3f \"%d\"\n",lp->x,-lp->y,lp->number);
		    }
		}
	    }
	}
    }


float minf2(x1,x2)
/* 
*   Computes minimum of two floating point numbers.
*/

float x1,x2;

{   return x1 < x2 ? x1 : x2;
    }

/****************************************************/

float maxf2(x1,x2)
/* 
*   Computes maximum of two floating point numbers.
*/

float x1,x2;

{   return x1 > x2 ? x1 : x2;
    }


int ask(question)
/*
*   Ask a yes/no question and return the answer as Boolean variable.
*/

char *question;

{
    char answer[2];
    char help[10];
    int done,ret;

    strcpy(help,"");
    done = false;
    while (!done) {
	printf("%s %s",question,help);
	scanf("%1s",answer);
	if (answer[0] == 'Y' || answer[0] == 'y') {
	    ret = true;
	    done = true;
	    }
	else if (answer[0] == 'N' || answer[0] == 'n') {
	    ret = false;
	    done = true;
	    }
	else {
	    strcpy(help,"(Y or N) ");
	    }
    }
    return ret;
}

void trim(st,stlenp)
/*
*   Trims blanks off the end of st. *stlenp gives the current length
*   of st.
*/

char *st;
int *stlenp;

{
    int stlen = *stlenp;

    while (stlen > 0) {
	if (st[stlen-1] != ' ') break;
	stlen -= 1;
	}
    *stlenp = stlen;
}

void chcnbl(st,stlen)
/*
*   Converts all non-printable control characters into blanks.
*/

char *st;
int stlen;

{
    int i;

    for (i = 1; i <= stlen; i++) {
	if (*st < 32 || *st > 126) *st = ' ';
	st += 1;
	}
}
