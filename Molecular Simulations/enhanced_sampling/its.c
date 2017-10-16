#include "its.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* may use specific lib */

int *its_error(int err_id)
{
    switch (err_id)
    {
        case 1:
            fprintf(stderr, "Allocation aborted.n"); break;
        case 2:
            fprintf(stderr, "Deallocation aborted.n"); break;
        default:
            break;
    }
    return NULL;
}

struct its_struct *its_init(
    double objtemp,
    int nb, double highT, double lowT,
    char *fb_filename, char *norml_filename
)
{
    int i;
    double dtemp;
    struct its_struct *obj;
    double p;
    FILE *fb_file, *norml_file;

    if ((obj = (struct its_struct *)malloc(sizeof(struct its_struct))) == NULL)
        return (struct its_struct *)its_error(1);
    obj->mcycle = 0;
    obj->nb = nb;
    obj->highT = highT;
    obj->lowT = lowT;
    obj->fb_filename = fb_filename;
    obj->norml_filename = norml_filename;
    obj->beta0 = exp(300.0 / 273.15 / objtemp);
    /* may be wrong... */

    if ((obj->fb = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->norml = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->normlold = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->mybeta = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->gf = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->ratio = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->pratio = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->rb = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->rbfb = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    if ((obj->bgf = (double *)calloc(obj->nb, sizeof(double))) == NULL)
        return (struct its_struct *)its_error(1);
    

    dtemp = (obj->highT-obj->lowT) / (obj->nb-1);
    for (i = 0; i < obj->nb; i++)
        obj->mybeta[i] = 300.0 / (obj->lowT + dtemp*i) / objtemp;
        /* why such a expression..? */
    p = exp(-0.005*obj->nb);
    for (i = 0; i < obj->nb; i++)
        obj->fb[i] = p;
    for (i = 0; i < obj->nb-1; i++)
     /* why obj->nb-1? ...solved */
        obj->norml[i] = 0;

    fb_file = fopen(obj->fb_filename, "a");
    norml_file = fopen(obj->norml_filename, "a");
    for (i = 0; i < obj->nb-1; i++)
    {
        fprintf(fb_file, "%dt%fn", obj->mcycle, obj->fb[i]);
        fprintf(norml_file, "%dt%fn", obj->mcycle, obj->norml[i]);
    }
    fprintf(fb_file, "%dt%fn", obj->mcycle, obj->fb[i]);
    fclose(fb_file);
    fclose(norml_file);

    for (i = 0; i < obj->nb; i++)
    {
        obj->gf[i] = 0;
        obj->bgf[i] = 0;
    }
    return obj;
}

int its_free(struct its_struct *obj)
{
    free(obj->fb);
    free(obj->norml);
    free(obj->normlold);
    free(obj->mybeta);
    free(obj->gf);
    free(obj->ratio);
    free(obj->pratio);
    free(obj->rb);
    free(obj->rbfb);
    free(obj->bgf);
    free(obj);

    return 0;
}

int its_updatefb(struct its_struct * obj)
{
    int i;
    FILE *fb_file, *norml_file;

    obj->mcycle++;
    for (i = 0; i < obj->nb-1; i++)
    {
        obj->rb[i] = (obj->rbfb[i]+obj->rbfb[i+1]) * 0.5;
        obj->ratio[i] = obj->fb[i] - obj->fb[i+1];
        obj->normlold[i] = obj->norml[i];
    }
    if (obj->mcycle == 1)
        for (i = 0; i < obj->nb-1; i++)
        {
            obj->norml[i] = obj->rb[i];
            obj->normlold[i] = -10000000000;
        }
    else
        for(i = 0; i < obj->nb-1; i++)
            if (obj->norml[i] > obj->rb[i])
                obj->norml[i] = obj->norml[i] + log(1.0 + exp(-obj->rb[i]+obj->norml[i]));
            else
                obj->norml[i] = obj->rb[i] + log(1.0 + exp(-obj->rb[i]+obj->norml[i]));

    for (i = 0; i < obj->nb-1; i++)
        if (obj->normlold[i] > obj->rbfb[i+1] - obj->rbfb[i] + obj->rb[i])
            obj->ratio[i] += obj->normlold[i] - obj->norml[i] + log(1.0 +
                exp(-obj->normlold[i] + obj->rbfb[i+1] - obj->rbfb[i] + obj->rb[i]));
        else
            obj->ratio[i] = obj->rbfb[i+1] - obj->rbfb[i] + obj->rb[i] + obj->ratio[i] - obj->norml[i] +
                log(1.0 + exp(obj->normlold[i] - obj->rbfb[i+1] + obj->rbfb[i] - obj->rb[i]));

    obj->pratio[0] = 0.0;
    for (i = 0; i < obj->nb-1; i++)
        obj->pratio[i+1] = obj->pratio[i] + obj->ratio[i];
    for (i = 0; i < obj->nb; i++)
        obj->fb[i] = -obj->pratio[i];

    fb_file = fopen(obj->fb_filename, "a");
    norml_file = fopen(obj->norml_filename, "a");
    for (i = 0; i < obj->nb-1; i++)
    {
        fprintf(fb_file, "%dt%fn", obj->mcycle, obj->fb[i]);
        fprintf(norml_file, "%dt%fn", obj->mcycle, obj->norml[i]);
    }
    fprintf(fb_file, "%dt%fn", obj->mcycle, obj->fb[obj->nb-1]);
    fclose(fb_file);
    fclose(norml_file);

    for (i = 0; i < obj->nb-1; i++)
        obj->rbfb[i] = 0.0;

    return 0;
}

int its_force(struct its_struct *obj, double v, double vshift, int step, double **f, int natom)
{
    int i;
    double vb;

    vb = v - vshift;
    for (i = 0; i < obj->nb; i++)
    {
        obj->gf[i] = -obj->mybeta[i]*vb + obj->fb[i];
        obj->bgf[i] = obj->gf[i] + log(obj->mybeta[i]);
    }

    obj->gfsum = obj->gf[0];
    obj->bgfsum = obj->bgf[0];
    for (i = 1; i < obj->nb; i++)
    {
        if (obj->bgfsum > obj->bgf[i])
            obj->bgfsum = obj->bgfsum + log(1.0 + exp(obj->bgf[i] - obj->bgfsum));
        else
            obj->bgfsum = obj->bgf[i] + log(1.0 + exp(obj->bgfsum - obj->bgf[i]));

        if (obj->gfsum > obj->gf[i])
            obj->gfsum = obj->gfsum + log(1.0 + exp(obj->gf[i] - obj->gfsum));
        else
            obj->gfsum = obj->gf[i] + log(1.0 + exp(obj->gfsum - obj->gf[i]));
    }

    if ((step == 1)||(step%100 == 0))
        for (i = 0; i < obj->nb; i++)
            obj->rbfb[i] = obj->gf[i];
    else
        for (i = 0; i < obj->nb; i++)
            if (obj->rbfb[i] > (obj->gf[i] - obj->gfsum))
                obj->rbfb[i] += log(1.0 + exp(obj->gf[i] - obj->gfsum - obj->rbfb[i]));
            else
                obj->rbfb[i] = obj->gf[i] - obj->gfsum + log(1.0 + exp(obj->rbfb[i] - obj->gf[i] + obj->gfsum));

    for (i = 0; i < 3*natom; i++) *f[i] = *f[i] * exp(obj->bgfsum - obj->gfsum) / obj->beta0;

    return 0;
}
