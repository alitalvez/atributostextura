#ifndef READIMAGE_H
#define READIMAGE_H

#include "globalvar.h"

class ReadImage
{
private:
    char *fileName   = NULL;
    int coluna;
    int linha;
public:
    ReadImage(char *fileName, unsigned short int coluna, unsigned short int linha);
    tImage vectorImage();
    tImage vectorImage_();
};

#endif // READIMAGE_H

