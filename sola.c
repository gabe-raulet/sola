#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <sndfile.h>

#ifndef TWOPI
#define TWOPI (6.283185307179586)
#endif

int read_input_wave_file(const char *fname, float **sig, int *chans, int *samps, int *rate);
int write_output_wave_file(const char *fname, const float *sig, int chans, int samps, int rate);
int sola(float *insig, int insamps, float **outsig, int *outsamps, int N, int Sa, double alpha, int cor);

int main(int argc, char *argv[])
{
    int Sa = 256; /* analysis hop size */
    int N = 2048; /* block length */
    double alpha = 1.; /* time scaling factor */
    int cor = 0;

    const char *infile = NULL;
    const char *outfile = NULL;

    int c;
    while ((c = getopt(argc, argv, "i:o:S:N:a:Ch")) >= 0)
    {
        if      (c == 'i') infile = optarg;
        else if (c == 'o') outfile = optarg;
        else if (c == 'S') Sa = atoi(optarg);
        else if (c == 'N') N = atoi(optarg);
        else if (c == 'a') alpha = atof(optarg);
        else if (c == 'C') cor = 1;
        else if (c == 'h')
        {
            fprintf(stderr, "Usage: %s [options] -i <input.wav> -o <output.wav>\n", argv[0]);
            fprintf(stderr, "Options: -S INT   analysis hop size [%d]\n", Sa);
            fprintf(stderr, "         -N INT   block length [%d]\n", N);
            fprintf(stderr, "         -a FLOAT time scaling factor [%.2f]\n", alpha);
            fprintf(stderr, "         -C       use auto-correlation\n");
            fprintf(stderr, "         -h       help message\n");
            return -1;
        }
    }

    if (!infile) { fprintf(stderr, "error: missing -i <input.wav>\n"); return 1; }
    if (!outfile) { fprintf(stderr, "error: missing -o <output.wav>\n"); return 1; }

    assert((alpha >= 0.25 && alpha <= 2.));

    int channels, insamps, outsamps, inrate, outrate;
    float *insig, *outsig;


    int err = read_input_wave_file(infile, &insig, &channels, &insamps, &inrate);

    if (err < 0)
    {
        fprintf(stderr, "error: unable to read input file '%s'\n", infile);
        return -1;
    }

    assert((channels == 1));

    sola(insig, insamps, &outsig, &outsamps, N, Sa, alpha, cor);

    if (err < 0)
    {
        fprintf(stderr, "error: something went wrong during resampling\n");
        return -1;
    }

    err = write_output_wave_file(outfile, outsig, channels, outsamps, inrate);

    if (err < 0)
    {
        fprintf(stderr, "error: unable to write output file '%s'\n", outfile);
        return -1;
    }

    return 0;
}

int read_input_wave_file(const char *fname, float **sig, int *chans, int *samps, int *rate)
{
    if (!fname || !sig || !chans || !samps || !rate)
        return -1;

    SNDFILE *sf;
    SF_INFO info;

    sf = sf_open(fname, SFM_READ, &info);
    if (!sf) return -1;

    int frames = info.frames;
    int channels = info.channels;
    int samplerate = info.samplerate;

    float *buf = malloc(channels*frames*sizeof(float));

    if (!buf)
    {
        sf_close(sf);
        return -1;
    }

    sf_readf_float(sf, buf, frames);
    sf_close(sf);

    *sig = buf;
    *chans = channels;
    *samps = frames;
    *rate = samplerate;

    return 0;
}

int write_output_wave_file(const char *fname, const float *sig, int chans, int samps, int rate)
{
    if (!fname || !sig)
        return -1;

    SNDFILE *sf;
    SF_INFO info;

    info.channels = chans;
    info.samplerate = rate;
    info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

    sf = sf_open(fname, SFM_WRITE, &info);
    if (!sf) return -1;

    sf_writef_float(sf, sig, samps);
    sf_close(sf);

    return 0;
}

void correlate(float *r, float *x1, float *x2, int L)
{
    for (int k = 0; k < L; ++k)
    {
        double val = 0.;

        for (int n = 0; n < L-k; ++n)
        {
            val += x1[n]*x2[n+k];
        }

        r[k] = val/L;
    }
}

int argmax(float *r, int L)
{
    int i;
    int rmax=0;
    double v = DBL_MIN;

    for (i = 0; i < L; ++i)
    {
        if (r[i] > v)
        {
            v = r[i];
            rmax = i;
        }
    }

    return rmax;
}

void overlap(float *yp, float *xp, int N, int L)
{
    double cur = 0;
    double delta = 1./L;

    int i = 0;

    while (i < L)
    {
        yp[i] *= (1. - cur);
        yp[i] += cur*xp[i];
        cur += delta;
        i++;
    }

    while (i < N)
    {
        yp[i] = xp[i];
        i++;
    }
}

int sola(float *insig, int insamps, float **outsig, int *outsamps, int N, int Sa, double alpha, int cor)
{
    if (!insig || !outsig || !outsamps)
        return -1;

    int Ss = (int)(Sa*alpha);
    int L = (int)(256*(alpha/2.));
    int M = (insamps + Sa - 1) / Sa;

    int n = N + (M-1)*Sa;
    int m = N + (M-1)*Ss;

    float *x = calloc(n, sizeof(float));
    float *y = calloc(m, sizeof(float));
    float *r = malloc(L*sizeof(float));

    memcpy(x, insig, insamps*sizeof(float));
    memcpy(y, x, N*sizeof(float));

    float *xp = x + Sa;
    float *yp = y + Ss;

    int km = 0;
    float *xl1, *xl2;

    for (int ni = 1; ni < M; ++ni)
    {
        if (cor)
        {
            correlate(r, yp, xp, L);
            km = argmax(r, L);
        }

        overlap(yp+km, xp, N, L);

        yp += Ss;
        xp += Sa;
    }

    *outsig = y;
    *outsamps = m;

    return 0;
}
