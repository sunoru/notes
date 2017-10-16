#ifndef ITS_H
#define ITS_H

struct its_struct
{
    int nb, mcycle;
    char *fb_filename;
    char *norml_filename;
    double gfsum, bgfsum, highT, lowT, beta0;
    double *fb;
    double *norml;
    double *normlold;
    double *mybeta;
    double *gf;
    double *ratio;
    double *pratio;
    double *rb;
    double *rbfb;
    double *bgf;
};

/*
    its_error: to output the error.
*/
int *its_error(int err_id);

/*
    its_init: to do the initialization.
    default:
    n=200, highT=400, lowT=280,
    fb_filename="fb.dat", norml_filename="norml.dat"
*/
struct its_struct *its_init(
    double objtemp,
    int n, double highT, double lowT,
    char *fb_filename, char *norml_filename
);

/*
    its_free: to free the obj.
*/
int its_free(struct its_struct *obj);

/*
    its_updatefb: to update the sampling.
*/
int its_updatefb(struct its_struct * obj);

/*
    its_force: to calculate the bias energy and forces.
*/
int its_force(struct its_struct *obj, double v, double vshift, int step, double **f, int natom);

#endif
