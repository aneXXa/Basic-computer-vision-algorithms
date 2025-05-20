#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <regex.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "filters.h"

int main(int argc, char *argv[]) {

    if (!strcmp(argv[1], "-help")) {
        printf("\nImage Processing Tool - Help\n");
        printf("==========================\n\n");
        
        printf("Basic Usage:\n");
        printf("  %s input_image <operation> [parameters] [-o output_image]\n\n", argv[0]);
        
        printf("Available Operations:\n");
        printf("  -gauss [sigma]       Apply Gaussian blur (default sigma=3.0)\n");
        printf("  -median [size]       Apply median filter (default size=3, must be odd)\n");
        printf("  -edge [low] [high]   Edge detection with Canny algorithm\n");
        printf("                       (default thresholds: low=0.5, high=0.6)\n");
        printf("  -grayscale           Convert image to grayscale\n");
        printf("  -sharpen             Apply sharpening filter\n");
        printf("  -ridge               Apply ridge detection\n");
        printf("  -histogram           Apply histogram equalization\n\n");
        
        printf("Examples:\n");
        printf("  %s photo.jpg -gauss 2.5 -o blurred.jpg\n", argv[0]);
        printf("  %s input.png -median 5\n", argv[0]);
        printf("  %s image.jpg -edge 0.2 0.4 -o edges.png\n\n", argv[0]);
        
        printf("Note:\n");
        printf("- For operations with parameters, if none are provided,\n");
        printf("  default values will be used\n");
        printf("- Supported input formats: JPG, PNG, BMP, TGA\n");
        printf("- Output format is determined by file extension (.jpg, .png)\n");
        return EXIT_SUCCESS;
    }

    // --- Check minimum argument requirements ---
    if (argc < 3) {
        fprintf(stderr, "Error: Incorrect syntax. Usage:\n");
        fprintf(stderr, "\t%s input_filename <operation> [operation_params] [-o output_filename]\n", argv[0]);
        fprintf(stderr, "\tAvailable operations:\n");
        fprintf(stderr, "\t-gauss [sigma] | -median [kernel_size] | -edge [low_thresh] [high_thresh]\n");
        fprintf(stderr, "\t-grayscale | -sharpen | -ridge | -histogram\n");
        return EXIT_FAILURE;
    }

    // --- Image Loading ---
    char *inputName = argv[1];
    int width, height, channels;
    unsigned char *img = stbi_load(inputName, &width, &height, &channels, 0);
    
    if (!img) {
        fprintf(stderr, "Error: Failed to load image '%s'. Supported formats: .png, .jpg, .jpeg\n", inputName);
        return EXIT_FAILURE;
    }
    printf("Successfully loaded image: %s (%dx%d pixels, %d channels)\n", inputName, width, height, channels);

    // --- Operation Selection and Parameter Handling ---
    char *operation = argv[2];
    
    if (strcmp(operation, "-gauss") == 0) {
        // Gaussian blur with optional sigma parameter
        float sigma = 3.0f; // default value
        
        if (argc > 3 && strcmp(argv[3], "-o") != 0) {
            sigma = atof(argv[3]);
            if (sigma == 0) {
                printf("Warning: Invalid sigma value (0). Using default value 3.0\n");
                sigma = 3.0f;
            } else if (sigma < 0) {
                printf("Warning: Negative sigma value (%.1f). Using absolute value\n", sigma);
                sigma = fabs(sigma);
            }
        }
        
        printf("Applying Gaussian blur with sigma=%.1f...\n", sigma);
        gauss_blur(img, width, height, channels, sigma);
        
    } else if (strcmp(operation, "-median") == 0) {
        // Median filter with optional kernel size
        int kernel = 3; // default value
        
        if (argc > 3 && strcmp(argv[3], "-o") != 0) {
            kernel = atoi(argv[3]);
            if (kernel == 0) {
                printf("Warning: Invalid kernel size (0). Using default 3x3 kernel\n");
                kernel = 3;
            } else if (kernel < 0) {
                printf("Warning: Negative kernel size (%d). Using absolute value\n", kernel);
                kernel = abs(kernel);
            } else if (kernel % 2 == 0) {
                printf("Warning: Even kernel size (%d). Using %dx%d kernel instead\n", kernel, kernel+1, kernel+1);
                kernel += 1;
            }
        }
        
        printf("Applying median filter with %dx%d kernel...\n", kernel, kernel);
        median_filter(img, width, height, channels, kernel);
        
    } else if (strcmp(operation, "-edge") == 0) {
        // Canny edge detection with optional thresholds
        float low_thresh = 0.5f;  // default
        float high_thresh = 0.6f; // default
        
        if (argc > 3 && strcmp(argv[3], "-o") != 0) {
            low_thresh = atof(argv[3]);
            if (argc > 4 && strcmp(argv[4], "-o") != 0) {
                high_thresh = atof(argv[4]);
            }
            
            // Validate thresholds
            if (low_thresh <= 0 || low_thresh > 1) {
                printf("Warning: Invalid low threshold (%.2f). Using default 0.1\n", low_thresh);
                low_thresh = 0.1f;
            }
            if (high_thresh <= 0 || high_thresh > 1) {
                printf("Warning: Invalid high threshold (%.2f). Using default 0.3\n", high_thresh);
                high_thresh = 0.3f;
            }
            if (low_thresh > high_thresh) {
                printf("Warning: Low threshold (%.2f) > high threshold (%.2f). Swapping values\n", low_thresh, high_thresh);
                float temp = high_thresh;
                high_thresh = low_thresh;
                low_thresh = temp;
            }
        }
        
        printf("Applying Canny edge detection with thresholds (%.2f, %.2f)...\n", low_thresh, high_thresh);
        edge_detecting_Canny(img, width, height, channels, low_thresh, high_thresh);
        
    } else if (strcmp(operation, "-grayscale") == 0) {
        // Simple grayscale conversion
        printf("Converting image to grayscale...\n");
        toGrayscale(img, width, height, channels);
        
    } else if (strcmp(operation, "-sharpen") == 0) {
        // Image sharpening
        printf("Applying sharpening filter...\n");
        sharpen(img, width, height, channels);
        
    } else if (strcmp(operation, "-ridge") == 0) {
        // Ridge detection
        printf("Applying ridge detection...\n");
        ridge(img, width, height, channels);
        
    } else if (strcmp(operation, "-histogram") == 0) {
        // Histogram equalization
        printf("Applying histogram equalization...\n");
        histogramEqualization(img, width, height, channels);
        
    } else {
        // Invalid operation specified
        fprintf(stderr, "Error: Unknown operation '%s'\n", operation);
        fprintf(stderr, "\tAvailable operations:\n");
        fprintf(stderr, "\t-gauss [sigma] | -median [kernel_size] | -edge [low_thresh] [high_thresh]\n");
        fprintf(stderr, "\t-grayscale | -sharpen | -ridge | -histogram\n");
        stbi_image_free(img);
        return EXIT_FAILURE;
    }

    // --- Image Saving ---
    char *outputName;
    char defaultName[] = "output.jpg"; // default output filename
    
    // Determine output filename
    if (argc >= 4 && strcmp(argv[argc - 2], "-o") == 0) {
        outputName = argv[argc - 1];
    } else {
        outputName = defaultName;
        printf("No output filename specified. Using default: %s\n", outputName);
    }

    // Handle file extension and save
    char *extension = strrchr(outputName, '.');
    if (extension) {
        if (strcmp(extension, ".png") == 0) {
            if (!stbi_write_png(outputName, width, height, channels, img, 0)) {
                fprintf(stderr, "Error: Failed to save PNG image\n");
            }
        } else if (strcmp(extension, ".jpg") == 0 || strcmp(extension, ".jpeg") == 0) {
            if (!stbi_write_jpg(outputName, width, height, channels, img, 100)) { // 100 = max quality
                fprintf(stderr, "Error: Failed to save JPG image\n");
            }
        } else {
            fprintf(stderr, "Error: Unsupported output format '%s'. Using JPG instead\n", extension);
            strcat(outputName, ".jpg");
            stbi_write_jpg(outputName, width, height, channels, img, 100);
        }
    } else {
        printf("No file extension specified. Using JPG format\n");
        strcat(outputName, ".jpg");
        stbi_write_jpg(outputName, width, height, channels, img, 100);
    }

    printf("Successfully processed and saved image to '%s'\n", outputName);
    stbi_image_free(img);
    return EXIT_SUCCESS;
}