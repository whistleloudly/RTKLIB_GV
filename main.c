#include <stdio.h>
#include "rtklib.h"

#define MAXFILE     16                  /* max number of input files */

static int nepoch = 0;            /* number of observation epochs */
static char proc_rov[64] = "";   /* rover for current processing */
static char proc_base[64] = "";   /* base station for current processing */

static nav_t navs = { 0 };          /* navigation data */


/* show message and check break ----------------------------------------------*/
static int checkbrk(const char* format, ...)
{
    va_list arg;
    char buff[1024], * p = buff;
    if (!*format) return showmsg("");
    va_start(arg, format);
    p += vsprintf(p, format, arg);
    va_end(arg);
    if (*proc_rov && *proc_base) sprintf(p, " (%s-%s)", proc_rov, proc_base);
    else if (*proc_rov) sprintf(p, " (%s)", proc_rov);
    else if (*proc_base) sprintf(p, " (%s)", proc_base);
    return showmsg(buff);
}

/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char** infile,
    const int* index, int n, const prcopt_t* prcopt,
    nav_t* nav)
{
    int i, ind = 0, nobs = 0, rcv = 1;

    trace(3, "readobsnav: ts=%s n=%d\n", time_str(ts, 0), n);

    nav->eph = NULL; nav->n = nav->nmax = 0;
    nav->geph = NULL; nav->ng = nav->ngmax = 0;
    nav->seph = NULL; nav->ns = nav->nsmax = 0;
    nepoch = 0;

    for (i = 0; i < n; i++) {
        /* read rinex obs and nav file */
        if (readrnxt(infile[i], rcv, ts, te, ti, prcopt->rnxopt[rcv <= 1 ? 0 : 1], nav) < 0) {
            checkbrk("error : insufficient memory");
            trace(1, "insufficient memory\n");
            return 0;
        }
    }

    if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
        checkbrk("error : no nav data");
        trace(1, "\n");
        return 0;
    }

    /* delete duplicated ephemeris */
    uniqnav(nav);

    return 1;
}

/* get exe folder ------------------------------------------------------------*/
static void getexefolder(char* pathfile, const char* exedir)
{
    char tempstr[_MAX_DIR] = "";
    char drive[_MAX_DRIVE] = "";
    char dir[_MAX_DIR] = "";
    char fname[_MAX_FNAME] = "";
    char ext[_MAX_EXT] = "";
    strcpy(tempstr, exedir);
    _splitpath(tempstr, drive, dir, fname, ext);
    strcpy(pathfile, drive);
    strcat(pathfile, dir);
}

int main(int argc, char** argv)
{
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = { "" };
    gtime_t ts = { 0 }, te = { 0 };
    double tint = 0.0, es[] = { 2000,1,1,0,0,0 }, ee[] = { 2000,12,31,23,59,59 };

    int svh[MAXSAT], n = MAXSAT;
    double* rs, * dts, * var, *azel, r, dop[4];
    rs = mat(6, n); dts = mat(2, n); var = mat(1, n);
    char* timestr = "2022 01 14 08 00 00";
    gtime_t time;
    str2time(timestr, 0, 19, &time);

    /* wmk20140220 */
    char* infile[MAXFILE] = { "" }, outfile[_MAX_DIR] = { "/0" };
    char cfgfile[_MAX_DIR] = "";

    FILE* fp1, * fp2;
    char nsatfile[_MAX_DIR] = "", dopfile[_MAX_DIR] = "";

    /* wmk20140220   obs and nav file path*/
    char FilePath[_MAX_DIR] = "";
    int ifile = 0;
    /* wmk20140220 allocate memory */
    for (ifile = 0; ifile < MAXFILE; ifile++)
    {
        infile[ifile] = malloc(_MAX_DIR * sizeof(char));
    }
    infile[0] = "E:\\地大研究生文件\\RTKLIB_Code_Read\\2.43\\rtklib2.4.3_b34\\算例\\rinex\\tohy0213\\22_0113\\BRDM00DLR_S_20220130000_01D_MN.rnx";
    int* index = {0};


    /* load options from configuration file */
    getexefolder(cfgfile, argv[0]);
    strcat(cfgfile, "opt.conf");
    resetsysopts();
    if (!loadopts(cfgfile, sysopts)) return -1;
    getsysopts(&prcopt, &solopt, &filopt);


    getexefolder(nsatfile, argv[0]);
    strcat(nsatfile, "nsat.txt");
    if (fp1 = fopen(nsatfile, "w+"))
        printf("nsfile open");
    else
        printf("nsfile open error!");

    getexefolder(dopfile, argv[0]);
    strcat(dopfile, "dop.txt");
    if (fp2 = fopen(dopfile, "w+"))
        printf("dopfile open");
    else
        printf("dopfile open error!");

    /* read obs and nav data */
    if (!readobsnav(ts, te, tint, infile, index, 1, &prcopt_default, &navs)) return 0;

    //rtkinit(&rtk, popt);


    //char satid[4];
    //for (int i = 0; i < MAXSAT; i++) {
    //    satno2id(i+1, satid);
    //    int prn;
    //    int sys = satsys(i + 1, &prn);
    //    printf("sat number:%d sat id:%s sat prn:%d satsys:%d \n", i + 1, satid, prn, sys);
    //}

    /* satellite positons, velocities and clocks */
    satposs(time, n, &navs, prcopt.sateph, rs, dts, var, svh);

    /*生成经纬度格网点*/
    int ns = 0;
    double lat = 91.0,lon = -182.0;
    while (lat > -90) {
        lat = lat - 1; 
        while (lon < 180) {
            lon = lon + 2;
            ns = 0;
            double refpos[3] = { lat ,lon ,10 }, refpos_rad[3] = { lat * PI / 180,lon * PI / 180, 80 }, refecef[3] = { 0,0,0 }, e[3] = { 0,0,0 };
            azel = zeros(2, n);
            pos2ecef(refpos_rad, refecef);

            /* excluded satellite? */
            for (int i = 0; i < MAXSAT; i++) {

                //azel[i * 2] = azel[1 + i * 2] = 0.0;

                if (satexclude(i + 1, var[i], svh[i], &prcopt)) continue;

                /* geometric distance */
                if ((r = geodist(rs + i * 6, refecef, e)) <= 0.0) continue;

                /* test elevation mask */
                if (satazel(refpos_rad, e, azel + i * 2) < prcopt.elmin) continue;
                ns++;

                char satid[4];
                satno2id(i + 1, satid);
                //printf("the retain sat: %d %s azel: %8.4f %8.4f opt.elmin: %8.4f \n", i + 1, satid, azel[i * 2] * 180 / 3.14, azel[1 + i * 2] * 180 / 3.14, prcopt.elmin * 180 / 3.14);
            }



            dops(ns, azel, prcopt.elmin, dop);
            if (dop[1] <= 0.0 || dop[1] > prcopt.maxgdop) {
                dop[1] = 0;
            }
            printf("refpos:%f %f %f  ns:%d dop:%f\n", refpos[0],refpos[1],refpos[2], ns, dop[1]);
            //fprintf(fp, "%8.4f %8.4f %8.4f %3d %4.3f\n", refpos[0], refpos[1], refpos[2], ns, dop[1]);
            //if (lon == 180) {
            //    fprintf(fp, "%8.4f %8.4f\n", refpos[1], refpos[0]);
            //}
            //else {
            //    fprintf(fp,"%8.4f %8.4f ", refpos[1], refpos[0]);
            //}
            if (lon == 180) {
                if (fp1 && fp2) {
                    fprintf(fp1, "%3d\n", ns);
                    fprintf(fp2, "%5.3f\n", dop[1]);
                }
            }
            else {
                if (fp1 && fp2) {
                    fprintf(fp1, "%3d\t", ns);
                    fprintf(fp2, "%5.3f\t", dop[1]);
                }
            }
            //sol->dop[0] = dop[0];   /* GDOP */
            //sol->dop[1] = dop[1];   /* PDOP */
            //sol->dop[2] = dop[2];   /* HDOP */
            //sol->dop[3] = dop[3];   /* VDOP */
        }
        lon = -182.0;
    }


	return 0;
}