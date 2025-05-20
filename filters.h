#ifndef FILTERS_H
#define FILTERS_H

void median_filter (unsigned char *img, int width, int height, int channels, int kernel);
void gauss_blur (unsigned char *img, int width, int height, int channels, float sigma);
void toGrayscale (unsigned char *img, int width, int height, int channels);
void edge_detecting_Canny (unsigned char *img, int width, int height, int channels, float down, float up);
void sharpen (unsigned char *img, int width, int height, int channels);
void ridge (unsigned char *img, int width, int height, int channels);
void histogramEqualization (unsigned char *img, int width, int height, int channels);

#endif