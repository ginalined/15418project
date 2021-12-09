#include <algorithm>
#include "VCScene.h"
#include "image.h"
#include <iostream>
#include <string>

const int DATA_DUMP = 0;
const int RENDER_DUMP = 1;

void renderPicture();

static struct {
    int width;      // frame width
    int height;     // frame height

    VCScene* vs;
} gDisplay;

void
renderPicture() {

    // clear screen
    gDisplay.vs->clearImage();

    // render the particles< into the image
    gDisplay.vs->render();
}

void dumpData(const std::string& frameFilename, int frame, int num_tri) {
    char filename[1024];
    sprintf(filename, "%s_%04d.inp", frameFilename.c_str(), frame);
    
    FILE *fp = fopen(filename, "wb");

    if (!fp) {
        fprintf(stderr, "Error: could not open %s for write\n", filename);
        exit(1);
    }

    fprintf(fp, "%lf\n\n", (double)num_tri);
    gDisplay.vs->dumpTriangles(fp, num_tri);

    fclose(fp);
    // printf("Wrote image file %s\n", filename);
}

void dumpFrame(const std::string& frameFilename, int frame) {
    char filename[1024];
    sprintf(filename, "%s_%04d.ppm", frameFilename.c_str(), frame);
    
    Image *image = gDisplay.vs->getImage();
    if (!image) {
        std::cout <<"FATAL: gDisplay.vs->getImage() return null" << std::endl;
    }

    FILE *fp = fopen(filename, "wb");

    if (!fp) {
        fprintf(stderr, "Error: could not open %s for write\n", filename);
        exit(1);
    }

    // write ppm header
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", gDisplay.width, gDisplay.height);
    fprintf(fp, "255\n");

    // std::cout << "gDisplay.height = " << gDisplay.height << ", gDisplay.width = " << gDisplay.width << std::endl;

    for (int j=gDisplay.height-1; j>=0; j--) {
        for (int i=0; i<gDisplay.width; i++) {
            const float* ptr = &image->data[4 * (j*image->width + i)];

            char val[3];
            val[0] = static_cast<char>(255.f * CLAMP(ptr[0], 0.f, 1.f));
            val[1] = static_cast<char>(255.f * CLAMP(ptr[1], 0.f, 1.f));
            val[2] = static_cast<char>(255.f * CLAMP(ptr[2], 0.f, 1.f));

            fputc(val[0], fp);
            fputc(val[1], fp);
            fputc(val[2], fp);
        }
    }

    fclose(fp);
    // printf("Wrote image file %s\n", filename);
}

void
startRendererWithDisplay(VCScene* vs, int option, const std::string& frameFilename, int frame, int num_tri) {
    // setup the display
    vs->allocateImage(128, 128);
    const Image* img = vs->getImage();
    
    gDisplay.vs = vs;
    gDisplay.width = img->width;
    gDisplay.height = img->height;

    if (option == RENDER_DUMP) {
        renderPicture();
        dumpFrame(frameFilename, frame);
    } else if (option == DATA_DUMP) {
        dumpData(frameFilename, frame, num_tri);
    }
}
