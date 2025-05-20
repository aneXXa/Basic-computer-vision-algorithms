#include "filters.h"

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// --- Helper functions ---

// Comparison function for qsort
int compare(const void* a, const void* b) {
    return (*(unsigned char*)a - *(unsigned char*)b);
}

// --- Median Filter ---
/**
 * @brief Applies a median filter to an image for noise reduction.
 * 
 * @details This function processes an image by replacing each pixel with the median value
 *          in a neighborhood of the specified kernel size. It effectively removes "salt-and-pepper" noise.
 *          Edge pixels are handled by clamping coordinates (replicating border pixels).
 * 
 * @param[in,out] img   Pointer to the image pixel array in [R,G,B,...] or other format.
 *                      The data is overwritten with the filtered result.
 * @param[in] width     Image width in pixels.
 * @param[in] height    Image height in pixels.
 * @param[in] channels  Number of channels (e.g., 3 for RGB).
 * @param[in] kernel    Filter kernel size.
 * 
 * @note The kernel size must be an odd number. If even, behavior is undefined.
 * @warning The function allocates temporary memory of size `width * height * channels` bytes.
 *          On allocation failure, it prints an error to stderr and returns early.
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image (640x480)
 * median_filter(image, 640, 480, 3, 5); // 5x5 median filter
 * @endcode
 * 
 * @see The function relies on `compare()` for sorting pixel values (must be defined separately).
 */
void median_filter(unsigned char* img, int width, int height, int channels, int kernel) {
    int pad = kernel / 2;
    unsigned char* temp = malloc(width * height * channels * sizeof(unsigned char));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed in median_filter\n");
        return;
    }
    
    memcpy(temp, img, width * height * channels * sizeof(unsigned char));

    // Process all pixels including edges
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < channels; c++) {
                unsigned char window[kernel * kernel];
                int index = 0;
                int actual_pixels = 0;

                // Sample pixels from the kernel window
                for (int wy = -pad; wy <= pad; wy++) {
                    for (int wx = -pad; wx <= pad; wx++) {
                        int px = x + wx;
                        int py = y + wy;
                        
                        // Handle edge cases by clamping coordinates
                        px = (px < 0) ? 0 : ((px >= width) ? width - 1 : px);
                        py = (py < 0) ? 0 : ((py >= height) ? height - 1 : py);
                        
                        int pos = (py * width + px) * channels + c;
                        window[index] = temp[pos];
                        index++;
                        actual_pixels++;
                    }
                }

                // Sort only the actual pixels we sampled
                qsort(window, actual_pixels, sizeof(unsigned char), compare);
                
                // Store median value
                img[(y * width + x) * channels + c] = window[actual_pixels / 2];
            }
        }
    }
    free(temp);
}

// --- Gaussian Blur ---

/**
 * @brief Applies a 1D Gaussian blur to an image in either horizontal or vertical direction.
 * 
 * @details This is a helper function for separable Gaussian blur. It processes the image by convolving
 *          each pixel with a 1D Gaussian kernel in the specified direction (horizontal or vertical).
 *          Edge pixels are handled by clamping (only valid kernel positions contribute to the result).
 * 
 * @param[in] src           Input image pixel array (format: [R,G,B,...] or other format).
 * @param[out] dst          Output image pixel array (must be pre-allocated).
 * @param[in] width         Image width in pixels.
 * @param[in] height        Image height in pixels.
 * @param[in] channels      Number of channels (e.g., 3 for RGB).
 * @param[in] kernel        1D Gaussian kernel array (normalized).
 * @param[in] kernel_size   Size of the kernel (must be odd).
 * @param[in] is_horizontal Direction flag: 1 for horizontal pass, 0 for vertical.
 * 
 * @note The kernel must be normalized (sum of weights = 1.0). 
 * @warning `src` and `dst` must not overlap. For in-place filtering, use a temporary buffer.
 */
static void gauss_1d(unsigned char* src, unsigned char* dst, int width, int height, int channels, const float* kernel, int kernel_size, int is_horizontal) {
    int radius = kernel_size / 2;

    // Apply 1D Gaussian filter in specified direction
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < channels; c++) {
                float sum = 0.0f;
                float weight_sum = 0.0f;

                for (int k = -radius; k <= radius; k++) {
                    if (is_horizontal) {
                        int px = x + k;
                        if (px >= 0 && px < width) {
                            float weight = kernel[k + radius];
                            sum += src[(y * width + px) * channels + c] * weight;
                            weight_sum += weight;
                        }
                    } else {
                        int py = y + k;
                        if (py >= 0 && py < height) {
                            float weight = kernel[k + radius];
                            sum += src[(py * width + x) * channels + c] * weight;
                            weight_sum += weight;
                        }
                    }
                }
                dst[(y * width + x) * channels + c] = (unsigned char)(sum / weight_sum + 0.5f);
            }
        }
    }
}

/**
 * @brief Applies a Gaussian blur to an image using a separable kernel.
 * 
 * @details This function smooths an image by convolving it with a 2D Gaussian kernel, implemented
 *          efficiently as two 1D passes (horizontal + vertical). The kernel size is determined by
 *          the `sigma` parameter (3*sigma rule). Memory is allocated for a temporary buffer.
 * 
 * @param[in,out] img      Pointer to the image pixel array in [R,G,B,...] or other format.
 *                         The data is overwritten with the filtered result.
 * @param[in] width        Image width in pixels.
 * @param[in] height       Image height in pixels.
 * @param[in] channels     Number of channels (e.g., 3 for RGB).
 * @param[in] sigma        Standard deviation of the Gaussian kernel (controls blur strength).
 * 
 * @note Larger `sigma` values produce stronger blur. Kernel size is auto-calculated as `6*sigma + 1`.
 * @warning Allocates a temporary buffer of size `width * height * channels`. On allocation failure,
 *          prints an error to stderr and returns early.
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image (640x480)
 * gauss_blur(image, 640, 480, 3, 2.0f); // Blur with sigma=2.0
 * @endcode
 */
void gauss_blur(unsigned char* img, int width, int height, int channels, float sigma) {
    // Calculate kernel size (3*sigma rule)
    int kernel_size = (int)(sigma * 3.0f) * 2 + 1;
    if (kernel_size < 3) kernel_size = 3;

    unsigned char* temp = malloc(width * height * channels * sizeof(unsigned char));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed in gauss_blur\n");
        return;
    }

    // Generate Gaussian kernel
    int radius = kernel_size / 2;
    float kernel[kernel_size];
    float sum = 0.0f;

    for (int i = -radius; i <= radius; i++) {
        float val = (1/sqrt(2* M_PI * sigma * sigma)) * exp(-(i * i) / (2.0f * sigma * sigma));
        kernel[i + radius] = val;
        sum += val;
    }

    // Normalize kernel
    for (int i = 0; i < kernel_size; i++) {
        kernel[i] /= sum;
    }

    // Apply separable Gaussian filter (horizontal then vertical)
    gauss_1d(img, temp, width, height, channels, kernel, kernel_size, 1);
    gauss_1d(temp, img, width, height, channels, kernel, kernel_size, 0);

    free(temp);
}

// --- Grayscale Conversion (YUV) ---

/**
 * @brief Converts an RGB image to grayscale using luminance (Y) component.
 * 
 * @details This function transforms an RGB image to grayscale by calculating the luminance (Y) value
 *          for each pixel using the standard CCIR 601 weights (0.299R + 0.587G + 0.114B) and storing
 *          the result in all three channels (R=G=B=Y). The result is a monochromatic image that
 *          preserves perceived brightness.
 * 
 * @param[in,out] img      Pointer to the image pixel array in [R,G,B,...] or other format.
 *                         The data is overwritten with the filtered result.
 * @param[in] width        Image width in pixels.
 * @param[in] height       Image height in pixels.
 * @param[in] channels     Number of channels (e.g., 3 for RGB).
 * 
 * @note 
 * - Uses ITU-R BT.601 luminance coefficients for RGB to grayscale conversion.
 * - If `channels > 3`, extra channels are left unchanged (e.g., alpha channel).
 * 
 * @warning 
 * - The function assumes RGB pixel order (channels 0=R, 1=G, 2=B).
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image
 * toYUV(image, 640, 480, 3); // Convert to grayscale
 * @endcode
 */
void toGrayscale(unsigned char* img, int width, int height, int channels) {
    for (int i = 0; i < width * height; i++) {
        unsigned char r = img[i * channels + 0];
        unsigned char g = img[i * channels + 1];
        unsigned char b = img[i * channels + 2];

        // Calculate luminance (Y component)
        unsigned char y = (unsigned char)(0.299f * r + 0.587f * g + 0.114f * b);

        // Set all channels to luminance value
        img[i * channels + 0] = y;
        img[i * channels + 1] = y;
        img[i * channels + 2] = y;
    }
}

// --- Canny Edge Detector ---

// Sobel operators for edge detection
const int SOBEL_X[3][3] = {
    {-3, 0, 3},
    {-10, 0, 10},
    {-3, 0, 3}
};

const int SOBEL_Y[3][3] = {
    {-3, -10, -3},
    {0, 0, 0},
    {3, 10, 3}
};

/**
 * @brief Performs Canny edge detection on an input image.
 * 
 * @details Implements the full Canny edge detection pipeline:
 *          1. Converts image to grayscale
 *          2. Applies Gaussian blur for noise reduction
 *          3. Computes gradient magnitude and direction using Sobel operators
 *          4. Performs non-maximum suppression
 *          5. Applies hysteresis thresholding with edge tracking
 * 
 * @param[in,out] img           Pointer to the image pixel array in [R,G,B,...] or other format.
 *                              The data is overwritten with the filtered result.
 * @param[in] width             Image width in pixels.
 * @param[in] height            Image height in pixels.
 * @param[in] channels          Number of channels (e.g., 3 for RGB).
 * @param[in] low_thresh        Low threshold for hysteresis (ratio of max gradient, 0-1).
 * @param[in] high_thresh       High threshold for hysteresis (ratio of max gradient, 0-1).
 * 
 * @note 
 * - Uses modified Sobel kernels for better gradient estimation
 * - Implements 8-directional edge tracking during hysteresis
 * 
 * @warning
 * - Requires input image to have at least 3 channels
 * - Threshold values are relative to maximum gradient magnitude
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image
 * edge_detecting_Canny(image, 640, 480, 3, 0.1f, 0.3f);
 * @endcode
 */
void edge_detecting_Canny(unsigned char* img, int width, int height, int channels, float low_thresh, float high_thresh) {
    // Pre-processing
    toGrayscale(img, width, height, channels);
    gauss_blur(img, width, height, channels, 1.0f);

    float* gradient_magnitude = malloc(width * height * sizeof(float));
    float* gradient_direction = malloc(width * height * sizeof(float));
    unsigned char* suppressed = calloc(width * height, sizeof(unsigned char));

    if (!gradient_magnitude || !gradient_direction || !suppressed) {
        fprintf(stderr, "Error: Memory allocation failed in edge_detecting_Canny\n");
        free(gradient_magnitude);
        free(gradient_direction);
        free(suppressed);
        return;
    }

    // Gradient calculation using Sobel operators
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int gx = 0, gy = 0;

            // Convolution with Sobel kernels
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    unsigned char pixel = img[(y + ky) * width * channels + (x + kx) * channels];
                    gx += pixel * SOBEL_X[ky + 1][kx + 1];
                    gy += pixel * SOBEL_Y[ky + 1][kx + 1];
                }
            }

            gradient_magnitude[y * width + x] = sqrtf(gx * gx + gy * gy);
            gradient_direction[y * width + x] = roundf(atan2f(gy, gx) * (4.0 / M_PI)) * (M_PI / 4.0);
        }
    }

    // Non-maximum suppression
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int idx = y * width + x;
            float angle = gradient_direction[idx];
            float mag = gradient_magnitude[idx];

            // Determine neighbor pixels in gradient direction
            int dx = (int)round(cos(angle));
            int dy = -(int)round(sin(angle));

            // Compare with neighbors
            if ((x + dx >= 0 && x + dx < width && y + dy >= 0 && y + dy < height)) {
                if (mag < gradient_magnitude[(y + dy) * width + (x + dx)]) {
                    suppressed[idx] = 0;
                    continue;
                }
            }
            if ((x - dx >= 0 && x - dx < width && y - dy >= 0 && y - dy < height)) {
                if (mag < gradient_magnitude[(y - dy) * width + (x - dx)]) {
                    suppressed[idx] = 0;
                    continue;
                }
            }
            suppressed[idx] = (unsigned char)mag;
        }
    }

    // Hysteresis thresholding
    float max_mag = 0.0f;
    for (int i = 0; i < width * height; i++) {
        if (suppressed[i] > max_mag) max_mag = suppressed[i];
    }

    unsigned char low = (unsigned char)(max_mag * low_thresh);
    unsigned char high = (unsigned char)(max_mag * high_thresh);
    unsigned char *thresholded = malloc(width * height * sizeof(unsigned char));

    // Apply thresholds
    for (int i = 0; i < width * height; i++) {
        if (suppressed[i] >= high) {
            thresholded[i] = 255;       
        } else if (suppressed[i] >= low) {
            thresholded[i] = 127;       
        } else {
            thresholded[i] = 0;         
        }
    }

    // Edge tracking by hysteresis
    const int dx8[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    const int dy8[8] = {0, -1, -1, -1, 0, 1, 1, 1};
    
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int idx = y * width + x;
            if (thresholded[idx] == 127) {
                // Check 8-connected neighbors for strong edges
                for (int k = 0; k < 8; k++) {
                    int nx = x + dx8[k];
                    int ny = y + dy8[k];
                    if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                        if (thresholded[ny * width + nx] == 255) {
                            thresholded[idx] = 255;
                            break;
                        }
                    }
                }
                if (thresholded[idx] != 255) {
                    thresholded[idx] = 0;
                }
            }
        }
    }

    // Copy result back to original image
    for (int i = 0; i < width * height; i++) {
        unsigned char val = thresholded[i];
        for (int c = 0; c < channels; c++) {
            img[i * channels + c] = val;
        }
    }

    free(gradient_magnitude);
    free(gradient_direction);
    free(suppressed);
    free(thresholded);
}

// --- Sharpening Filter ---

/**
 * @brief Applies sharpening effect to an image using 3x3 convolution.
 * 
 * @details Enhances image edges by applying sharpening kernel
 *                  1. Creates a temporary copy of the original image
 *                  2. Applies convolution with the sharpening kernel
 *                  3. Clamps resulting values to [0, 255] range
 *                  4. Processes all channels independently
 * 
 * @param[in,out] img       Pointer to the image pixel array in [R,G,B,...] or other format.
 *                          The data is overwritten with the filtered result.
 * @param[in] width         Image width in pixels.
 * @param[in] height        Image height in pixels.
 * @param[in] channels      Number of channels (e.g., 3 for RGB).
 * 
 * @note 
 * - Uses zero-padding at image boundaries (skips edge pixels)
 * - Kernel weights sum to 1, preserving overall brightness
 * 
 * @warning
 * - Allocates temporary buffer of size width×height×channels
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image
 * sharpen(image, 640, 480, 3);
 * @endcode
 */
void sharpen(unsigned char* img, int width, int height, int channels) {
    // Sharpening kernel
    const int SHARPEN_KERNEL[3][3] = {
        {0, -1, 0},
        {-1, 5, -1},
        {0, -1, 0}
    };

    unsigned char* temp = malloc(width * height * channels * sizeof(unsigned char));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed in sharpen\n");
        return;
    }
    memcpy(temp, img, width * height * channels * sizeof(unsigned char));

    // Apply convolution with sharpening kernel
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            for (int c = 0; c < channels; c++) {
                int sum = 0;
                
                // 3x3 kernel convolution
                for (int ky = -1; ky <= 1; ky++) {
                    for (int kx = -1; kx <= 1; kx++) {
                        unsigned char pixel = temp[(y + ky) * width * channels + (x + kx) * channels + c];
                        sum += pixel * SHARPEN_KERNEL[ky + 1][kx + 1];
                    }
                }
                
                // Clamp values to [0, 255]
                sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum);
                img[y * width * channels + x * channels + c] = (unsigned char)sum;
            }
        }
    }

    free(temp);
}

// --- Ridge Detection ---

/**
 * @brief Performs ridge detection on an image using a 3x3 convolution kernel.
 * 
 * @details The function:         
 *                  1. Creates a temporary copy of the original image
 *                  2. Converts the image to grayscale
 *                  2. Applies convolution with the ridge kernel
 *                  3. Clamps resulting values to [0, 255] range
 *                  4. Processes all channels independently
 * 
 *                  This kernel enhances diagonal ridges and valleys in the image while
 *                  suppressing flat regions and orthogonal edges.
 * 
 * @param[in,out] img       Pointer to the image pixel array in [R,G,B,...] or other format.
 *                          The data is overwritten with the filtered result.
 * @param[in] width         Image width in pixels.
 * @param[in] height        Image height in pixels.
 * @param[in] channels      Number of channels (e.g., 3 for RGB).
 * 
 * @note
 * - Processes all channels but produces grayscale output (R=G=B)
 * - Kernel sum is zero (preserves uniform regions)
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3];  // RGB image
 * ridge(image, 640, 480, 3);      // Apply ridge detection
 * @endcode
 */
void ridge(unsigned char* img, int width, int height, int channels) {
    toGrayscale(img, width, height, channels);

    // Ridge detection kernel
    const int RIDGE_KERNEL[3][3] = {
        {-2, -1, 0},
        {-1, 0, 1},
        {0, 1, 2}
    };

    unsigned char* temp = malloc(width * height * channels * sizeof(unsigned char));
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed in ridge\n");
        return;
    }
    memcpy(temp, img, width * height * channels * sizeof(unsigned char));

    // Apply ridge detection kernel
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            for (int c = 0; c < channels; c++) {
                int sum = 0;
                
                for (int ky = -1; ky <= 1; ky++) {
                    for (int kx = -1; kx <= 1; kx++) {
                        unsigned char pixel = temp[(y + ky) * width * channels + (x + kx) * channels + c];
                        sum += pixel * RIDGE_KERNEL[ky + 1][kx + 1];
                    }
                }
                
                // Clamp values to [0, 255]
                sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum);
                img[y * width * channels + x * channels + c] = (unsigned char)sum;
            }
        }
    }

    free(temp);
}

// --- Histogram Equalization ---

/**
 * @brief Performs histogram equalization in HSV color space to enhance image contrast.
 * 
 * @details The function:
 *              1. Converts RGB image to HSV color space
 *              2. Equalizes the Value (V) channel using histogram stretching
 *              3. Preserves Hue (H) and Saturation (S) to maintain color integrity
 *              4. Converts back to RGB space
 * 
 * @param[in,out] img       Pointer to the image pixel array in [R,G,B,...] or other format.
 *                          The data is overwritten with the filtered result.
 * @param[in] width         Image width in pixels.
 * @param[in] height        Image height in pixels.
 * @param[in] channels      Number of channels (e.g., 3 for RGB).
 * 
 * @note
 * - Preserves color relationships by operating only on luminance (V channel)
 * - Handles empty histogram bins through minimum value detection
 * 
 * @warning
 * - Allocates temporary HSV buffer (3×width×height floats)
 * 
 * @code
 * // Usage example:
 * unsigned char image[640*480*3]; // RGB image
 * histogramEqualization(image, 640, 480, 3);
 * @endcode
 */
void histogramEqualization(unsigned char* img, int width, int height, int channels) { 
    // Temporary storage for HSV values
    float* hsv = malloc(width * height * 3 * sizeof(float));
    if (!hsv) {
        fprintf(stderr, "Error: Memory allocation failed in histogramEqualization\n");
        return;
    }

    // Convert RGB to HSV
    for (int i = 0; i < width * height; i++) {
        unsigned char r = img[i * channels + 0];
        unsigned char g = img[i * channels + 1];
        unsigned char b = img[i * channels + 2];
        
        float max_val = fmaxf(r, fmaxf(g, b)) / 255.0f;
        float min_val = fminf(r, fminf(g, b)) / 255.0f;
        float delta = max_val - min_val;
        
        // Calculate Hue
        float h = 0.0f;
        if (delta > 0.0f) {
            if (max_val == r / 255.0f) {
                h = 60.0f * fmodf((g - b) / 255.0f / delta, 6.0f);
            } else if (max_val == g / 255.0f) {
                h = 60.0f * (((b - r) / 255.0f) / delta + 2.0f);
            } else {
                h = 60.0f * (((r - g) / 255.0f) / delta + 4.0f);
            }
            if (h < 0.0f) h += 360.0f;
        }
        
        // Calculate Saturation
        float s = (max_val > 0.0f) ? (delta / max_val) : 0.0f;
        
        // Store HSV values
        hsv[i * 3 + 0] = h;
        hsv[i * 3 + 1] = s;
        hsv[i * 3 + 2] = max_val; // Value component
    }

    // Compute histogram for V channel
    const int HIST_SIZE = 256;
    int hist[HIST_SIZE] = {};
    for (int i = 0; i < width * height; i++) {
        int val = (int)(hsv[i * 3 + 2] * 255.0f);
        hist[val]++;
    }

    // Compute cumulative histogram
    int cdf[HIST_SIZE] = {};
    cdf[0] = hist[0];
    for (int i = 1; i < HIST_SIZE; i++) {
        cdf[i] = cdf[i - 1] + hist[i];
    }

    // Find minimum non-zero value in histogram
    int min_val = 0;
    while (min_val < HIST_SIZE && hist[min_val] == 0) min_val++;
    if (min_val == HIST_SIZE) min_val = 0;

    // Equalize V channel
    int total_pixels = width * height;
    for (int i = 0; i < width * height; i++) {
        float v = hsv[i * 3 + 2];
        int val = (int)(v * 255.0f);
        
        if (cdf[min_val] != total_pixels) {
            float scale = 1.0f / (total_pixels - cdf[min_val]);
            v = (cdf[val] - cdf[min_val]) * scale;
        }
        
        // Clamp and store back
        hsv[i * 3 + 2] = fmaxf(0.0f, fminf(1.0f, v));
    }

    // Convert HSV back to RGB
    for (int i = 0; i < width * height; i++) {
        float h = hsv[i * 3 + 0];
        float s = hsv[i * 3 + 1];
        float v = hsv[i * 3 + 2];
        
        float c = v * s;
        float x = c * (1.0f - fabsf(fmodf(h / 60.0f, 2.0f) - 1.0f));
        float m = v - c;
        
        float r, g, b;
        if (h < 60.0f) {
            r = c; g = x; b = 0.0f;
        } else if (h < 120.0f) {
            r = x; g = c; b = 0.0f;
        } else if (h < 180.0f) {
            r = 0.0f; g = c; b = x;
        } else if (h < 240.0f) {
            r = 0.0f; g = x; b = c;
        } else if (h < 300.0f) {
            r = x; g = 0.0f; b = c;
        } else {
            r = c; g = 0.0f; b = x;
        }
        
        // Store RGB values
        img[i * channels + 0] = (unsigned char)((r + m) * 255.0f);
        img[i * channels + 1] = (unsigned char)((g + m) * 255.0f);
        img[i * channels + 2] = (unsigned char)((b + m) * 255.0f);
    }

    free(hsv);
}